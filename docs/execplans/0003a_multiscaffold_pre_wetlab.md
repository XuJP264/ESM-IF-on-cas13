# ExecPlan 0003A: Multi-scaffold pre-wet-lab preparation

## Purpose and evidence boundary

Advance the real VI-D platform from one 6E9F scaffold to an auditable
multi-scaffold, multi-state Cas13d preparation framework. The stage will build
and map experimental structures, score sequences without averaging state
coordinates, mine public Cas13d variant/activity evidence, run a deliberately
small real local multi-model smoke, preregister candidate selection, and export
deterministic H100 refold inputs plus a fully tested mock ingest path.

This stage can support Level 0 data/code validation, Level 1 sequence novelty,
and Level 2 inverse-folding compatibility only. It prepares Level-3 jobs but
does not produce Level-3 evidence unless independent real structure predictions
are later ingested. No result here establishes a functional, active, effective,
or wet-lab-validated Cas13.

## Current state

The plan starts on 2026-08-01 (Asia/Shanghai) from clean `main` at exact commit
`b5cd7949288e6b7ce20ed5e3a277d129d9c157ba`. Stage 0002 provides a strict
6E9F/EsCas13d four-layer map, a real two-seed nine-method matrix, local genuine
ESM-IF1/ProteinMPNN/LigandMPNN checkpoints, a 4,070-sequence Atlas novelty
resource, a length-gated VI-D MSA/conservation table, and provider-neutral
refold export/ingest primitives. Its 18 candidates remain computational only.

Known constraints are an 8,188 MiB RTX 4060 Laptop GPU, CPU-oriented local
execution, no trustworthy Atlas direct-repeat orientation, and no authorization
to launch large AlphaFold 2/3 or Boltz jobs. Public structure and literature
sources must be official or lawfully open and must retain URLs, access dates,
licenses, and hashes. The requested EsCas13d, RfxCas13d/CasRx, DjCas13d, and
UrCas13d structure availability must be discovered rather than assumed.

## Milestones

### A0 — Source discovery and frozen scope

Search primary RCSB records and original/open publications for Cas13d
experimental structures and variant/activity evidence. Record all hits and
exclusions. Freeze at least three independent Cas13d scaffolds and six
scaffold-state units if the public record supports them; otherwise retain an
explicit evidence-backed blocker and maximize real independent work.

### A1 — Experimental Cas13d structure atlas

Download PDB/mmCIF and official metadata atomically, hash every file, and build
one scaffold table plus one state table. Parse protein/RNA chains, full and
coordinate sequences, missing residues, alternate locations, noncanonical
residues, chain breaks, RNA atoms, resolution, state, publication, provenance,
license, and pairwise parent identities. Exclusions remain visible.

### A2 — Multi-scaffold four-layer mapping and annotations

Generalize the Stage-0002 mapper to configurable PDB/scaffold/state inputs and
map each included protein chain through coordinate, full natural sequence, and
VI-D MSA columns. Add per-position domain/state, coordinate availability,
RNA-contact/second-shell, core, domain-interface, HEPN, conservation, entropy,
gap fraction, status, and confidence. Conservation is fail-closed unless the
mapping gate passes. Emit machine tables plus CSV/HTML manual review.

### A3 — Multi-state scoring and masks

Implement state-group validation, normalized weights, per-state score,
weighted mean, minimum, variance, rank consistency, and state-combination
selection. Build intersection hard, union risk, and variable/hinge masks from
mapped positions without averaging coordinates. Stop on missing states,
inconsistent residue mappings, inconsistent fixed tokens, or invalid weights.
Run genuine ESM-IF1 multi-state scoring where compatible state groups exist.

### A4 — Cas13d variant/activity evidence

Create a source-level manifest and a non-pooled record table from public
primary/supplementary evidence. Preserve assay, guide/target, cis/trans,
knockdown, expression, solubility, normalization, numbering, extraction method,
replicate evidence, comparability group, and missingness. Recover full mutant
sequences only when unambiguous. Deduplicate within source/assay/variant, map
numbering, and perform explicitly versioned label-threshold sensitivity.

### A5 — Retrospective descriptive benchmark

For recoverable variants, compute genuine or clearly labelled unavailable
ESM-IF1 delta score, ProteinMPNN/LigandMPNN scores where structure/context are
valid, conservation, contact class, and multi-state score. Keep assays separate.
Report effect estimates and intervals only when supported; do not force a
classifier or significance claim from insufficient or incomparable data.

### A6 — Real local multi-scaffold smoke

For every locally runnable scaffold, execute two seed blocks and a few genuine
proposals for common-mask ESM-IF1, conservation ESM-IF1, conservation-plus-RNA
ESM-IF1, ProteinMPNN, LigandMPNN, and the two-source consensus. Require one
common safety mask per scaffold, zero fixed violations, actual RNA context for
LigandMPNN when present, no mock formal rows, and explicit CPU/GPU provenance.
This is engineering/Level-2 coverage, not final method inference.

### A7 — Candidate preregistration and Level-3 job export

Freeze `pre_wetlab_candidate_protocol.yaml` before viewing future refold
outputs. Require Level-1 QC, Atlas threshold, zero hard violations, multiple
inverse-folding support, no severe state failure, composition QC, sequence
cluster diversity, and conservative/moderate/aggressive novelty strata. Export
the current 18 candidates, WT controls, and new real smoke candidates into
deterministic monomer/binary/ternary manifests, shards, expected outputs, and
retry manifests for ColabFold/AF2, AF3, and Boltz. Target-scaffold templates are
forbidden as forced templates.

### A8 — Level-3 ingest and ranking fixture E2E

Using only `is_mock=true` prediction fixtures, validate structure/PAE/plDDT
ingest, US-align/TM-align invocation, domain RMSD, HEPN geometry, RNA-contact
recovery, interface confidence, multi-seed and cross-model consistency, missing
output audits, and Pareto ranking. Formal ranking must not collapse to a single
ESM score. No fixture metric supports a scientific claim.

### A9 — Clean H100 bundle, acceptance, and publication record

Commit all executable inputs before bundle export. Export a clean bundle from
the final execution commit with Stage-0003 manifests included, empty missing
assets, passing internal and large-asset hashes, and the single H100 tmux entry.
Run all required local gates and real smokes, publish the report, push without
force, and verify GitHub Actions.

## Progress

- [x] 2026-08-01: Read `AGENTS.md`, `.agent/PLANS.md`, Stage-0002 summary,
  experiment protocol, claims/evidence policy, and GPU migration guide in full;
  confirmed clean synchronized starting HEAD `b5cd794`.
- [x] 2026-08-01: Created this living ExecPlan and the Stage-0003A summary
  skeleton before source discovery, implementation, or new downloads.
- [x] 2026-08-01: Completed A0 source discovery using official RCSB records,
  open primary articles, and lawful public supplements; froze four independent
  Cas13d parents and nine experimental state units.
- [x] 2026-08-01: Completed A1 structure atlas/QC with official URLs, access
  date, licenses, hashes, full/coordinate sequences, chain assignments, RNA
  atoms, construct restorations, pairwise parent identity, and inclusion audit.
- [x] 2026-08-01: Completed A2 for all nine states. Strict weighted four-layer
  high-confidence coverage is 7,507/8,272 (90.75%); conservation is fail-closed
  outside the mapping gate and state-specific manual-review CSV/HTML exists.
- [x] 2026-08-01: Completed A3 with genuine CUDA ESM-IF1 scores for nine states,
  13 applicable state combinations, 11 explicit non-applicable combinations,
  and intersection/union/hinge masks with zero fixed-token violations.
- [x] 2026-08-01: Completed A4–A5. The real variant table has 22 records in
  four assay groups; the descriptive same-backbone benchmark scored 10 point
  variants with real ESM-IF1, ProteinMPNN, and RNA-context LigandMPNN and
  excluded six indels.
- [x] 2026-08-01: Completed A6 with 48 genuine local candidates over four
  scaffolds, six methods, and two seeds; runtime failures, mocks, and fixed
  violations are zero and all eight LigandMPNN rows retain RNA context.
- [x] 2026-08-01: Completed A7–A8. Candidate thresholds were frozen before
  Level-3 outputs; 70 proteins expand to 1,068 deterministic jobs, and four
  labeled mock predictions passed the complete ingest/TM-align/ranking E2E.
- [ ] Complete A9 final bundle, acceptance, push, and CI.

## Decisions

- “Independent scaffold” means a distinct natural Cas13d parent sequence, not
  another state or construct of the same parent. States remain separate
  statistical/structural units.
- Experimental states are never coordinate-averaged. Scores are evaluated on
  individual backbones and combined only after per-state scoring.
- A public paper naming a protein is not evidence that its atomic structure is
  deposited. Missing Rfx/Dj/Ur structures will be reported rather than replaced
  by predicted structures under an experimental label.
- VI-D conservation remains a proposal annotation only after a high-confidence
  scaffold-specific full-to-MSA map. It never automatically establishes a HEPN
  catalytic residue.
- Variant assay values are stratified by comparability group. Cross-assay
  normalization does not create a pooled biological scale.
- The local 8 GB GPU may be used only for bounded real smoke when safe; large
  AF2/AF3/Boltz inference is deferred as `not_run`.

## Discoveries

- RCSB supports four independent requested parent sequences and nine useful
  states: EsCas13d 6E9E/6E9F, UrCas13d 6IV9, DjCas13d
  9M38/9M30/9M33/9M34, and CasRx 9M31/9M8Q.
- The 2025 primary paper states that the DjCas13d and CasRx structures use
  four-residue nuclease-inactive constructs even though the corresponding RCSB
  entity metadata omits that mutation declaration. Natural proxies therefore
  restore those tokens and remain medium confidence.
- UrCas13d's RCSB-linked UniProt accession currently returns no sequence. The
  RCSB 922-residue reference plus primary-paper catalytic restoration is used;
  the eight-residue expression tag is excluded and provenance remains medium.
- An initial real multi-state run imported user-site PyTorch 2.12/cu130 and
  failed despite selecting the project analysis interpreter. Enforcing
  `PYTHONNOUSERSITE=1` allowed the isolated ESM environment (PyTorch 2.4.1,
  cu121) to use the RTX 4060 successfully. Both failed and repaired runs remain.
- LigandMPNN writes `use_ligand_context=True` on its input FASTA header, while
  generated headers carry sampling statistics. A first retrospective run
  correctly stopped because it inspected only the final header; the repaired
  gate requires the explicit context attestation anywhere in the upstream
  record set and then completed with zero forced-sequence violations.
- The numerical retrospective panel contains only nine assayed point variants;
  directions are mixed and class counts are insufficient for significance.
  Conservation has descriptive rho 0.5667 (cis) and 0.7448 (trans), but no
  model score establishes predictive validity.
- The Ur parent contains all absolute low-complexity windows flagged in its 12
  smoke rows. The selection protocol therefore fails on newly introduced
  windows relative to the parent baseline, rather than rejecting every design
  for inheriting a natural parent segment.

## Execution details

Work from `/home/junpeng/ESM-IF`. Use project-local analysis, ESM-IF1,
LigandMPNN, and bioinformatics environments. New remote assets are downloaded
to temporary files, validated, then atomically renamed. Raw structures and
large supplementary files remain ignored; manifests and processed tables are
versioned.

Authoritative planned outputs are:

```text
reports/stage_0003a/
data/processed/cas13d_variant_activity.parquet
data/processed/cas13d_variant_activity_sources.csv
docs/CAS13D_VARIANT_DATA_CARD.md
experiments/preregistered/pre_wetlab_candidate_protocol.yaml
docs/PRE_WETLAB_SELECTION_PROTOCOL.md
artifacts/gpu_jobs/stage_0003/
configs/stage_0003_refold.yaml
```

Formal generation uses seed blocks `20260801` and `20260802`, genuine local
checkpoints, and offline inference. Historical Atlas and Stage-0002 pilot/run
directories are read-only inputs and are not rebuilt.

## Validation and acceptance

Unit/fixture coverage must include state weight normalization, missing-state
failure, residue-map inconsistency, fixed-token inconsistency, identical-state
aggregation, state ranking, intersection/union/hinge masks, variant numbering,
deduplication, assay separation, refold job determinism, output missingness,
mock labelling, and Pareto inputs.

Required final commands are:

```bash
make lint
make typecheck
make test
make smoke-cpu
make smoke-esm-if1
make smoke-proteinmpnn
make smoke-ligandmpnn
make verify-reproducibility
```

Acceptance additionally requires three real independent Cas13d parent
scaffolds and six experimental state units if primary public structures exist;
all included mappings/contact annotations; genuine per-scaffold local model
smokes where backend inputs are valid; a real variant table and descriptive
benchmark or an exact source-data blocker; deterministic Level-3 inputs; mock
ingest E2E; zero mock formal rows; zero fixed-position violations; and no Level
3/4 claim without corresponding real evidence.

## Artifacts and provenance

Every download records official URL, retrieval UTC date, byte size, SHA256,
license/access status, and inclusion decision. Every real model run uses an
immutable `results/runs/<run_id>/` record with resolved config, command,
environment, hardware, git state, input/output manifests, logs, failure file,
exit code, and `SUCCESS`/`FAILED`. Canonical reports summarize but do not
replace run provenance.

## Blockers and recovery

Public structure or assay availability may prevent the requested scaffold or
variant coverage. After exhausting RCSB primary metadata, original papers, and
lawful supplements, record the exact missing evidence and continue with all
independent real scaffolds/interfaces. GPU scale is not a blocker for job
export/ingest fixtures. Atlas direct-repeat orientation remains irrelevant to
this stage unless a coevolution constraint is explicitly attempted.

## Outcomes and retrospective

A0–A8 are implemented and locally executed. The real ledger contains four
parents, nine structures/mappings/native state scores, a 22-record variant
resource, a descriptive 10-point-variant multi-model benchmark, and a 48-row
multi-scaffold real generation smoke. The mock ledger contains four explicitly
labeled prediction fixtures used only for Level-3 I/O and metric validation.
The 1,068 production refold jobs are `prepared_not_run`; real Level-3 count is
zero. A9 remains active until full acceptance, final clean bundle, push, and CI
pass. Without real Level-3 refolds, the wet-lab readiness answer remains “not
yet”.
