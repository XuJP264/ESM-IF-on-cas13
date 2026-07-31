#!/usr/bin/env python
"""Exercise the complete Level-3 ingest/ranking path with labeled fixtures."""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np

from cas13_if.provenance import RunRecorder, atomic_write_text, sha256_file
from cas13_if.refold.level3 import (
    consistency_summary,
    contact_recovery,
    domain_rmsd,
    hepn_geometry,
    interface_confidence,
    interface_pae,
    load_pae,
    pareto_front,
    run_usalign,
    validate_level3_result,
)
from cas13_if.refold.providers import ManifestPredictionProvider, PredictionJob
from cas13_if.structures.contacts import annotate_rna_contacts
from cas13_if.structures.parser import (
    group_residues,
    parse_structure,
    residue_polymer_type,
)


def _direct_contacts(path: Path) -> set[tuple[int, str]]:
    atoms = [
        atom
        for atom in parse_structure(path)
        if atom.residue.chain_id == "A"
        or (
            atom.residue.chain_id in {"B", "C"}
            and residue_polymer_type(atom.residue.residue_name) == "rna"
        )
    ]
    return {
        (item.protein_residue.residue_number, item.protein_residue.insertion_code)
        for item in annotate_rna_contacts(atoms)
        if item.direct_rna_contact
    }


def _ca_coordinates(path: Path) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    residues = group_residues(parse_structure(path))
    ordered: list[np.ndarray] = []
    by_number: dict[int, np.ndarray] = {}
    for key, atoms in residues.items():
        if key.chain_id != "A" or residue_polymer_type(key.residue_name) != "protein":
            continue
        ca = next((atom.coordinate for atom in atoms if atom.name == "CA"), None)
        if ca is not None:
            coordinate = np.asarray(ca, dtype=float)
            ordered.append(coordinate)
            by_number[key.residue_number] = coordinate
    return np.asarray(ordered), by_number


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    config: dict[str, Any] = {
        "experiment": "stage-0003a-refold-mock-e2e",
        "is_mock": True,
        "source_structure": "data/experimental_structures/6e9f.pdb",
        "alignment_executable": ".tools/envs/bioinformatics/bin/TMalign",
    }
    artifact_root = repo / "artifacts/refold_mock_e2e/stage_0003a"
    report_root = repo / "reports/stage_0003a/refold_mock_e2e"
    if artifact_root.exists() or report_root.exists():
        raise FileExistsError("refusing to overwrite Stage-0003A fixture E2E")
    recorder = RunRecorder(
        root=repo / "results/runs",
        experiment=str(config["experiment"]),
        resolved_config=config,
        command=[sys.executable, str(Path(__file__).resolve())],
        repo_root=repo,
        is_mock=True,
    )
    try:
        source = repo / str(config["source_structure"])
        sequence = "ACDE"
        ingested_rows: list[dict[str, Any]] = []
        comparisons = []
        provider_qc: list[dict[str, Any]] = []
        for provider_name in ("alphafold3", "boltz"):
            provider = ManifestPredictionProvider(provider_name, is_mock=True)
            jobs = [
                PredictionJob(
                    candidate_id=f"fixture-{provider_name}-seed{seed}",
                    sequence=sequence,
                    provider=provider_name,
                    seed=seed,
                    shard=0,
                    is_mock=True,
                )
                for seed in (1, 2)
            ]
            result_root = artifact_root / provider_name
            for job in jobs:
                directory = result_root / job.candidate_id
                directory.mkdir(parents=True, exist_ok=False)
                prediction = directory / "prediction.pdb"
                shutil.copyfile(source, prediction)
                pae_path = directory / "pae.json"
                atomic_write_text(
                    pae_path,
                    json.dumps(
                        {
                            "predicted_aligned_error": [
                                [0.0, 3.0, 5.0, 6.0],
                                [4.0, 0.0, 4.0, 5.0],
                                [5.0, 4.0, 0.0, 2.0],
                                [6.0, 5.0, 2.0, 0.0],
                            ],
                            "is_mock": True,
                        },
                        sort_keys=True,
                    )
                    + "\n",
                )
                result = {
                    "candidate_id": job.candidate_id,
                    "provider": provider_name,
                    "mean_plddt": 80.0 + job.seed,
                    "per_residue_plddt": [80.0, 82.0, 78.0, 84.0],
                    "structure_path": "prediction.pdb",
                    "pae_path": "pae.json",
                    "seed": job.seed,
                    "is_mock": True,
                }
                validate_level3_result(result, expected_mock=True)
                atomic_write_text(
                    directory / "result.json",
                    json.dumps(result, indent=2, sort_keys=True) + "\n",
                )
            ingested = provider.ingest_outputs(jobs, result_root)
            provider_qc.append(provider.qc_outputs(ingested))
            for prediction in ingested:
                if prediction.status != "success" or not prediction.is_mock:
                    raise RuntimeError("fixture prediction failed strict ingest")
                comparison = run_usalign(
                    Path(str(prediction.structure_path)),
                    source,
                    executable=repo / str(config["alignment_executable"]),
                    is_mock=True,
                )
                if comparison.status != "success":
                    raise RuntimeError(
                        f"real TM-align fixture comparison failed: {comparison}"
                    )
                comparisons.append(comparison)
                ingested_rows.append(
                    {
                        **prediction.__dict__,
                        "tm_score": comparison.tm_score,
                        "rmsd": comparison.rmsd,
                        "is_mock": True,
                    }
                )
        pae = load_pae(artifact_root / "alphafold3/fixture-alphafold3-seed1/pae.json")
        coordinates, hepn_coordinates = _ca_coordinates(source)
        domains = {
            "N_terminal_fixture": list(range(0, min(100, len(coordinates)))),
            "C_terminal_fixture": list(
                range(max(0, len(coordinates) - 100), len(coordinates))
            ),
        }
        rmsds = domain_rmsd(coordinates.copy(), coordinates, domains)
        geometry = hepn_geometry(hepn_coordinates, [(295, 300), (849, 854)])
        reference_contacts = _direct_contacts(source)
        contacts = contact_recovery(reference_contacts, reference_contacts.copy())
        interface_pae_value = interface_pae(pae, {0, 1}, {2, 3})
        interface_plddt = interface_confidence([80.0, 82.0, 78.0, 84.0], {0, 1})
        multi_seed = consistency_summary(comparisons[:2])
        cross_model = consistency_summary([comparisons[0], comparisons[2]])
        pareto_rows = [
            {
                "candidate_id": "fixture-a",
                "sequence_novelty": 0.70,
                "monomer_structural_recovery": 0.95,
                "multi_state_compatibility": 0.80,
                "rna_interface_preservation": contacts["recall"],
                "model_agreement": 0.90,
                "candidate_diversity": 0.60,
                "is_mock": True,
            },
            {
                "candidate_id": "fixture-b",
                "sequence_novelty": 0.85,
                "monomer_structural_recovery": 0.90,
                "multi_state_compatibility": 0.82,
                "rna_interface_preservation": 0.90,
                "model_agreement": 0.88,
                "candidate_diversity": 0.75,
                "is_mock": True,
            },
        ]
        dimensions = [
            "sequence_novelty",
            "monomer_structural_recovery",
            "multi_state_compatibility",
            "rna_interface_preservation",
            "model_agreement",
            "candidate_diversity",
        ]
        pareto = pareto_front(pareto_rows, maximize=dimensions, minimize=[])
        summary = {
            "ingested_predictions": len(ingested_rows),
            "validated_outputs": len(ingested_rows),
            "provider_qc": provider_qc,
            "alignment": {
                "executable": str(config["alignment_executable"]),
                "comparison_count": len(comparisons),
                "minimum_tm_score": min(
                    float(item.tm_score)
                    for item in comparisons
                    if item.tm_score is not None
                ),
            },
            "domain_rmsd": rmsds,
            "hepn_geometry": geometry,
            "rna_contact_recovery": contacts,
            "interface_pae": interface_pae_value,
            "interface_confidence": interface_plddt,
            "multi_seed_consistency": multi_seed,
            "cross_model_consistency": cross_model,
            "pareto_front": pareto,
            "real_prediction_count": 0,
            "is_mock": True,
            "evidence_level": 0,
            "claim_scope": "fixture pipeline validation only",
        }
        if any(not row["is_mock"] for row in ingested_rows + pareto_rows):
            raise RuntimeError("unlabeled fixture row entered mock E2E")
        report_root.mkdir(parents=True, exist_ok=False)
        atomic_write_text(
            report_root / "ingested_predictions.jsonl",
            "".join(json.dumps(row, sort_keys=True) + "\n" for row in ingested_rows),
        )
        atomic_write_text(
            report_root / "pareto_fixture.jsonl",
            "".join(json.dumps(row, sort_keys=True) + "\n" for row in pareto_rows),
        )
        summary_path = report_root / "summary.json"
        atomic_write_text(
            summary_path, json.dumps(summary, indent=2, sort_keys=True) + "\n"
        )
        recorder.finish(
            success=True,
            metrics=summary,
            outputs=[
                {
                    "path": str(summary_path.relative_to(repo)),
                    "sha256": sha256_file(summary_path),
                }
            ],
        )
        print(json.dumps({**summary, "run_dir": str(recorder.run_dir)}, indent=2))
        return 0
    except Exception as error:
        recorder.record_failure("stage_0003a_refold_mock_e2e", str(error))
        recorder.finish(success=False)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
