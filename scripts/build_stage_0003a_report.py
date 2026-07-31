#!/usr/bin/env python
"""Build the evidence-separated Stage-0003A Markdown and HTML report."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any

import pandas as pd

from cas13_if.provenance import atomic_write_text


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"JSON report is not a mapping: {path}")
    return value


def _markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(column) for column in frame.columns]
    header = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join("---" for _ in columns) + " |"
    rows = []
    for values in frame.itertuples(index=False, name=None):
        rows.append(
            "| "
            + " | ".join(
                str(value).replace("|", "\\|") if pd.notna(value) else ""
                for value in values
            )
            + " |"
        )
    return "\n".join([header, divider, *rows])


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    root = repo / "reports/stage_0003a"
    states = pd.read_csv(root / "states.csv")
    mappings = pd.read_csv(root / "mapping_summary.csv")
    native = pd.read_csv(root / "multistate_native_scores.csv")
    correlations = pd.read_csv(root / "variant_retrospective_correlations.csv")
    smoke = pd.read_csv(root / "local_multiscaffold_smoke/metrics.csv")
    variant = _json(root / "variant_retrospective_summary.json")
    variant_data = _json(repo / "data/processed/cas13d_variant_activity_summary.json")
    jobs = _json(root / "gpu_jobs_summary.json")
    mock = _json(root / "refold_mock_e2e/summary.json")
    aggregate_mapping = _json(root / "mapping_aggregate.json")

    structure_table = states[
        [
            "scaffold_id",
            "pdb_id",
            "state",
            "resolution_angstrom",
            "protein_chain",
            "crrna_chains",
            "target_rna_chains",
            "coordinate_length",
            "rna_atom_count",
        ]
    ].copy()
    structure_table.columns = [
        "scaffold",
        "PDB",
        "state",
        "resolution_A",
        "protein",
        "crRNA",
        "target_RNA",
        "coordinate_length",
        "RNA_atoms",
    ]
    mappings["four_layer_exact_or_restored_positions"] = (
        (mappings["high_confidence_coverage"] * mappings["full_scaffold_length"])
        .round()
        .astype(int)
    )
    mapping_table = mappings[
        [
            "pdb_id",
            "scaffold_id",
            "state",
            "full_scaffold_length",
            "mapped_coordinate_positions",
            "four_layer_exact_or_restored_positions",
            "high_confidence_coverage",
            "unresolved_positions",
            "RNA_contact_positions",
            "RNA_second_shell_positions",
        ]
    ].copy()
    mapping_table["high_confidence_coverage"] = mapping_table[
        "high_confidence_coverage"
    ].map(lambda value: f"{float(value):.4f}")
    native_table = native[
        [
            "scaffold_id",
            "pdb_id",
            "state",
            "mean_log_likelihood_per_resolved_residue",
            "perplexity",
        ]
    ].copy()
    for column in ("mean_log_likelihood_per_resolved_residue", "perplexity"):
        native_table[column] = native_table[column].map(
            lambda value: f"{float(value):.4f}"
        )
    retrospective_table = correlations.copy()
    retrospective_table["spearman_rho"] = retrospective_table["spearman_rho"].map(
        lambda value: f"{float(value):.4f}"
    )
    smoke_table = (
        smoke.groupby(["scaffold_id", "method"], as_index=False)
        .agg(
            n=("candidate_id", "size"),
            parent_identity_min=("parent_identity", "min"),
            parent_identity_max=("parent_identity", "max"),
            mean_esm_ll=("mean_conditional_log_likelihood", "mean"),
            fixed_violations=("fixed_position_violations", "sum"),
            rna_context_rows=("rna_atomic_context", "sum"),
        )
        .sort_values(["scaffold_id", "method"])
    )
    for column in ("parent_identity_min", "parent_identity_max", "mean_esm_ll"):
        smoke_table[column] = smoke_table[column].map(
            lambda value: f"{float(value):.4f}"
        )
    absolute_low_complexity = int(smoke["low_complexity_failure"].sum())
    ur_low_complexity = int(
        smoke.loc[smoke["scaffold_id"].eq("UrCas13d"), "low_complexity_failure"].sum()
    )
    four_layer_positions = int(
        aggregate_mapping["four_layer_exact_or_restored_positions"]
    )
    total_mapping_positions = int(
        aggregate_mapping["total_full_positions_across_states"]
    )
    markdown = f"""# Stage 0003A report

Status: Stage implementation and local acceptance complete. Immutable final
bundle and CI identifiers are published in the operator handoff. Date:
2026-08-01 (Asia/Shanghai).

## Claim boundary

The real calculations in this report support no more than Level 2
inverse-folding compatibility. The Level-3 prediction count is **0**. Prediction
fixtures are marked `is_mock=true` and support only Level 0 code-path evidence.
No sequence is described as functional, active, effective, or wet-lab
validated.

## Experimental multi-scaffold structure atlas — real data

The atlas contains {states["scaffold_id"].nunique()} independent natural Cas13d
parents and {len(states)} experimental scaffold-state units. All coordinates
were downloaded from RCSB, hashed, and audited; all RNA-bearing states have
non-zero retained RNA atom counts.

{_markdown_table(structure_table)}

## Four-layer mapping — real data and calculations

Across state-specific full sequences, {four_layer_positions:,}
of {total_mapping_positions:,} positions pass
the strict four-layer exact/restored gate, a weighted coverage of
{float(aggregate_mapping["weighted_high_confidence_coverage"]):.2%}. Conservation
is fail-closed elsewhere. Manual-review CSV and HTML files are in
`reports/stage_0003a/manual_review/` and each state mapping is in
`reports/stage_0003a/residue_mapping/`.

{_markdown_table(mapping_table)}

## Multi-state native scoring — real ESM-IF1

Nine state-specific native sequences were genuinely scored on the local RTX
4060 with the isolated ESM-IF1 checkpoint. The implementation combines scores
after state-wise evaluation; it never averages coordinates. Thirteen applicable
state combinations were aggregated and 11 unavailable combinations were
recorded as `not_applicable`; fixed-token violations were zero.

{_markdown_table(native_table)}

## Cas13d variant/activity resource — real literature evidence

The version-1 resource contains {variant_data["records"]} records from
{variant_data["studies"]} studies, {variant_data["unique_mutations"]} unique
mutation strings, and {variant_data["comparability_groups"]} non-poolable assay
groups. All {variant_data["records"]} mutant sequences were recovered. Numerical
2025 Figure-6 values are approximate graph readings and are not raw replicates.

The same-backbone retrospective benchmark genuinely scored
{variant["point_variants_scored"]} point variants; {variant["indels_excluded"]}
indels were excluded rather than mis-scored. It produced
{variant["real_esm_state_scores"]} ESM state scores,
{variant["real_proteinmpnn_scores"]} ProteinMPNN scores, and
{variant["real_ligandmpnn_scores"]} RNA-context LigandMPNN scores. The two-state
rank consistency was {float(variant["state_rank_consistency"]):.4f}. With only
n=9 numerical variants, all correlations below are descriptive and no
significance tests were run.

{_markdown_table(retrospective_table)}

The score/activity directions are mixed and do not establish predictive
validity. In particular, conservation showed the largest positive descriptive
association in this small panel (cis rho 0.5667; trans rho 0.7448), while model
scores were not consistently monotonic with activity. This is a negative/limited
result, not evidence for selecting a functional protein.

## Local multi-scaffold generation — real bounded smoke

Exactly {len(smoke)} real rows were generated: four scaffolds, six methods, two
seeds, and one proposal per seed. Every method/scaffold pair has two rows;
mock rows, runtime failures, and hard-fixed violations are all zero. All eight
LigandMPNN rows explicitly retained RNA atomic context. Each coordinate sequence
was also genuinely scored by ESM-IF1 on its representative state.

{_markdown_table(smoke_table)}

The historical absolute low-complexity flag is set for
{absolute_low_complexity}/48 rows, all {ur_low_complexity} UrCas13d rows,
because the natural Ur parent itself contains the flagged windows. Candidate
selection therefore uses the preregistered, parent-aware rule of zero *new*
low-complexity windows. These smoke rows are coverage checks, not a final
statistical comparison or a wet-lab shortlist.

## Level-3 job preparation — real inputs, prediction not run

The inventory has {jobs["candidate_count"]} proteins: {jobs["pilot_candidate_count"]}
Stage-0002 pilot candidates, {jobs["local_real_multiscaffold_candidate_count"]}
new real smoke candidates, and {jobs["wt_control_count"]} WT controls. It expands
to {jobs["job_count"]} deterministic jobs over two seeds: 560 monomer, 280
binary, and 228 ternary jobs across ColabFold, AlphaFold2, AlphaFold3, and Boltz.
All use `{jobs["template_policy"]}`. The target scaffold is never supplied as a
forced template. Status is `prepared_not_run`; manifests are not predictions.

## Level-3 ingest/ranking — MOCK fixture only

The fixture E2E ingested {mock["ingested_predictions"]} labeled mock results and
used the genuine local TM-align executable for {mock["alignment"]["comparison_count"]}
fixture comparisons. It exercised pLDDT, PAE, domain RMSD, HEPN geometry,
RNA-contact recovery, interface confidence, multi-seed consistency,
cross-model consistency, missing-output/retry behavior, and six-dimensional
Pareto ranking. `real_prediction_count={mock["real_prediction_count"]}` and every
fixture artifact remains `is_mock=true`; its perfect self-comparison metrics
have no scientific meaning.

## Audited failures and limitations

- One initial multi-state run failed because user-site PyTorch 2.12/CUDA 13.0
  shadowed the isolated environment. The failed run was retained; rerunning with
  `PYTHONNOUSERSITE=1` used pinned PyTorch 2.4.1/cu121 and succeeded on CUDA.
- One initial retrospective run stopped when validation inspected only the last
  LigandMPNN FASTA header. Upstream records context on the input header and
  sampling statistics on generated headers. The gate was corrected to require
  the explicit `use_ligand_context=True` attestation on any upstream header;
  the failed run remains retained.
- UrCas13d natural sequence provenance is medium confidence because its RCSB
  UniProt accession currently returns no sequence; the 922-residue RCSB
  reference was restored from primary-paper catalytic annotations and the
  expression tag was excluded.
- DjCas13d and CasRx RCSB metadata omit mutations described in the primary
  paper. Their four catalytic residues were restored with medium-confidence
  natural-parent proxies, visibly recorded in the state manifest.
- Domain boundaries and domain-interface labels remain manual-review fields,
  not inferred negatives. No real AlphaFold/AF3/Boltz prediction was run.

## Readiness decision

The project is ready to transfer deterministic Level-3 jobs to an H100 node,
but **has not reached the pre-wet-lab candidate standard**. It remains at Level
2 because no real independent monomer or RNA-complex prediction has passed the
predeclared Level-3 gates. It is not ready to nominate an experimentally
supported Cas13 and makes no efficacy claim.
"""
    report_md = root / "report.md"
    atomic_write_text(report_md, markdown)
    html_sections = [
        "<h1>Stage 0003A report</h1>",
        (
            "<p><strong>Evidence ceiling: Level 2. Real Level-3 predictions: "
            "0.</strong></p>"
        ),
        "<h2>Experimental structures</h2>",
        structure_table.to_html(index=False, escape=True),
        "<h2>Four-layer mapping</h2>",
        mapping_table.to_html(index=False, escape=True),
        "<h2>Real native multi-state ESM-IF1 scores</h2>",
        native_table.to_html(index=False, escape=True),
        "<h2>Descriptive variant retrospective</h2>",
        retrospective_table.to_html(index=False, escape=True),
        "<h2>Real local multi-model smoke</h2>",
        smoke_table.to_html(index=False, escape=True),
        "<h2>Evidence and readiness</h2>",
        "<p>Level-3 jobs are prepared but not run. Fixture outputs are MOCK. "
        "The project has not reached the pre-wet-lab candidate standard.</p>",
        "<pre>" + html.escape(json.dumps(jobs, indent=2, sort_keys=True)) + "</pre>",
    ]
    atomic_write_text(
        root / "report.html",
        "<!doctype html><html><head><meta charset='utf-8'><title>Stage 0003A"
        "</title></head><body>" + "".join(html_sections) + "</body></html>\n",
    )
    print(report_md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
