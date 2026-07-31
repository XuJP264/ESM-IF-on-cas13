# ExecPlan 0002: VI-D Mapping and Matched Baselines

## Purpose and evidence boundary

Advance the existing platform from a runnable real pilot to an auditable,
fair-method-comparison platform on one experimental VI-D scaffold: EsCas13d
chain A from PDB 6E9F. This plan establishes the strict residue-to-MSA mapping,
one common design mask, a preregistered small real CPU matrix, matched novelty
selection, hierarchical descriptive statistics, and a GPU/HPC continuation
manifest. It does not repeat the completed Atlas production parse or the
historical 72-candidate pilot, and it does not launch a large GPU task.

Mapping, code, and I/O verification support Level 0. Statistical novelty can
support Level 1 only when every registered gate passes. Genuine local
inverse-folding scores or samples can support Level 2 compatibility with the
input backbone. No output here is Level 3 without real refolding/multimodel
structural support, and no computational output is Level 4 or an effective or
validated Cas13.

## Current state

The plan starts on 2026-07-31 (Asia/Shanghai) from clean `main` at exact commit
`a9a530d14434e74dc0cfc47896847e201431c1c2`. The complete Atlas parse, exact
deduplication, six MMseqs2 cluster thresholds, strict 40% split, VI-D MSA, and
historical ESM-IF1 pilot already exist and will be consumed read-only.

6E9F chain A has a 954-residue deposited/full scaffold sequence and 864
resolved coordinate residues. The canonical VI-D MSA contains 182
cluster-representative sequences and 3,524 columns; 724 columns have at least
0.8 sequence coverage. Conservation is not yet permitted to enter any design
mask because the scaffold-to-MSA mapping and its confidence audit do not yet
exist.

Real CPU ESM-IF1, 6E9F ProteinMPNN, and 6E9F LigandMPNN ligand-checkpoint
smokes pass. The 5XWP MPNN smokes and LigandMPNN protein/soluble checkpoint
smokes have not yet run. The local WSL GPU device is unavailable, so this plan
uses preregistered small CPU execution and exports larger work without claiming
it ran.

Two GitHub Actions runs have different evidentiary roles and must never be
conflated: `30633893318` is the primary fixture/code CI run; `30635299503` is
the final handoff/documentation CI run for the preceding stage.

## Milestones

### P0 — Immutable handoff checkpoint

From exact clean HEAD `a9a530d14434e74dc0cfc47896847e201431c1c2`, export a
new GPU bundle. Acceptance requires `git_worktree_dirty_at_export=false`, an
empty `missing_assets`, exact manifest commit equality, all internal hashes,
and all 14 registered large-asset hashes. Create annotated tag
`v0.1.0-data-pipeline` only for this verified data-pipeline checkpoint and
record the two preceding CI roles in maintained documentation.

### P1 — Strict 6E9F/VI-D coordinate system

Build a deterministic four-layer mapping:

```text
PDB residue ID (chain, residue number, insertion code)
<-> resolved coordinate-sequence index
<-> full deposited scaffold-sequence index
<-> original VI-D MSA column
```

The mapper must preserve and classify terminal truncation, internal missing
coordinates, insertion codes, unresolved residues, query insertion columns,
internal MSA gaps, and residue disagreements. Each full-scaffold position must
carry mapping status, confidence, MSA coverage, weighted conservation,
weighted entropy, and gap fraction. It must validate that adding the scaffold
does not alter the original MSA columns. A machine table plus human-review CSV
and HTML are required. No conservation-derived hard or soft constraint is
released unless the registered mapping-confidence gate passes.

### P2 — One fair mask and small real baseline matrix

Use only resolved 6E9F chain-A positions with complete N/CA/C coordinates as
the canonical comparison universe. Translate this universe explicitly to every
backend's numbering. Every method receives the same common hard-fixed set, the
same free set, the same selected seeds, the same final candidate count, and the
same Atlas novelty search parameters. Any method-specific conservation or RNA
information may change proposal probabilities or allowed residues only within
the common free set; it may not silently change which positions are evaluated.

The common hard set comprises structurally non-designable positions plus
manually reviewed catalytic residues that map unambiguously. To satisfy both
the requested method list and the invariant mask, `unconstrained_esm_if1`
means no evolutionary/RNA proposal bias beyond this common safety mask, while
`catalytic_only_fixed_esm_if1` is a same-seed deterministic technical control
using the same common mask. Their equality is a reproducibility/control
contrast, not evidence for a distinct biological treatment. This limitation
must be visible in the methods table and report.

The preregistered local matrix uses two independent seeds and selects one real
candidate per seed per method after deterministic matching. Generation may
oversample within each seed, but final selected rows are balanced. The nine
requested methods are matched random mutation, MSA-profile sampling,
unconstrained ESM-IF1, catalytic-only fixed ESM-IF1,
conservation-constrained ESM-IF1, conservation-plus-RNA-contact ESM-IF1,
ProteinMPNN, LigandMPNN, and ESM-IF1/LigandMPNN consensus. Mock outputs are
prohibited from the formal matrix.

Parent identity and designed-position identity are matched by a preregistered
common binning/assignment procedure, not by inspecting performance endpoints.
If no nonempty common bin supports all methods, the run fails the formal
comparison and exports a GPU/HPC continuation manifest rather than silently
changing masks or counts. Atlas novelty uses the exact-unique 4,070-sequence
resource, MMseqs2 sensitivity 7.5, minimum query coverage 0.8, and the same
failure-closed no-hit semantics used previously.

### P3 — Metrics, statistics, failures, and report

For each selected candidate record conditional log-likelihood, perplexity,
parent and designed-position identity, maximum Atlas identity or explicit
no-coverage-hit status, fixed-position violations, low-complexity and
homopolymer gates, buried/core, RNA-interface, RNA-second-shell, and HEPN
recovery, model agreement, and diversity. All failed proposals and selection
failures remain in the funnel.

Candidates from one scaffold are not independent biological replicates.
Statistics therefore use seed as the paired resampling block and describe
candidate distributions without candidate-level biological p-values. Produce
paired bootstrap confidence intervals, paired effect sizes, matched-novelty
summaries, and Benjamini-Hochberg adjusted exploratory comparisons. With only
two local seeds, inferential results must be labelled low-power descriptive
evidence; the exported GPU manifest increases seeds without changing the
registered analysis.

### P4 — Missing real CPU model checks

Run and retain genuine CPU evidence for 5XWP ProteinMPNN, 5XWP LigandMPNN with
RNA atomic context, and the LigandMPNN upstream `protein_mpnn` and
`soluble_mpnn` checkpoints. A real failure is retained with command, exit code,
log, and exact recovery action; it is never replaced with mock output.

### P5 — Acceptance, handoff, and publication record

Generate all required files under `reports/matched_baselines/`, run the full
local quality gate, execute the formal small real matrix, verify identical
positions/zero fixed violations/matched identities/no mocks/evidence wording,
export the deterministic larger job manifest, update maintained status and
decision records, commit coherent verified changes, push without force, and
verify GitHub Actions.

## Progress

- [x] 2026-07-31: User explicitly resumed the long-term goal after the prior
  blocked state; inspected repository-wide instructions and confirmed a clean
  exact starting HEAD at `a9a530d14434e74dc0cfc47896847e201431c1c2`.
- [x] 2026-07-31: Created this ExecPlan before implementing the new mapping or
  matched-matrix features.
- [x] 2026-07-31: Exported
  `gpu-bundle-a9a530d14434-7540febfb2` from the exact clean starting HEAD;
  manifest commit equality, `dirty=false`, `missing_assets=[]`, all internal
  hashes, and 14/14 source-asset hashes passed. Created annotated tag
  `v0.1.0-data-pipeline` at that commit and documented the two prior CI roles.
- [x] 2026-07-31: Completed the real strict VI-D mapping. All 3,524 original
  MSA columns were preserved; 864/954 full-scaffold positions have exact
  four-layer mappings, 90 are coordinate-unresolved but full-to-MSA mapped,
  and no query-only MSA insertion occurred. Wrote the manual CSV/HTML audit.
- [x] 2026-07-31: Implemented the method-invariant mask, deterministic identity
  matching, seed-level statistics, environment-separated genuine model
  proposal/scoring runners, and MPNN adapters. Fixture tests cover the actual
  mapping writer, region builder, temporary biological-residue restoration,
  missing-slot translation, and fixed-position behavior.
- [x] 2026-07-31: Completed all four missing genuine CPU checks: 5XWP
  ProteinMPNN, 5XWP RNA-context LigandMPNN, and LigandMPNN's protein and
  soluble checkpoints. Every accepted artifact records CUDA unavailable for
  that forced-CPU process.
- [x] 2026-07-31: Completed the final preregistered real small matrix at
  `results/runs/20260731-vi-d-matched-baselines-3e2655746e-f5932f7/`.
  All nine methods contributed two selected candidates, all 18 received real
  ESM-IF1 scores, 14 passed Level-1 novelty/QC, mock count was zero, and the
  common four-fixed/860-free mask had zero violations.
- [x] 2026-07-31: Passed final local acceptance at `acd1137`: Ruff/format,
  strict mypy (47 source files), 61/61 tests with 72.17% branch coverage,
  fixture CPU smoke, matched-report invariants, bundle internal hashes, and
  14/14 source-asset hashes.
- [x] 2026-07-31: Pushed `main` and annotated tag without force. Initial
  stage-0002 CI run `30642823247` failed strict mypy in a fresh pip environment;
  reproduced the eight errors locally, added the missing pandas stubs and
  explicit NumPy dtype annotations, reran every local gate plus a fresh Python
  3.11 pip-only mypy check, and obtained successful code/results CI run
  `30643470631` at exact commit
  `b1a33fffd1d7a9254a55247c449017f240f02683`.

## Decisions

- The original VI-D MSA column numbering remains authoritative. Query-only
  insertion columns are represented explicitly with a null original-MSA
  column; the mapper must not renumber conservation data.
- The scaffold is added to the existing alignment with MAFFT's add-to-existing
  workflow. Original aligned rows are checked byte-for-byte after removing
  query insertion columns; a changed original alignment is a hard failure.
- Mapping confidence is rule-based and auditable: exact one-to-one links across
  all available layers are high; missing coordinates with an exact full-to-MSA
  mapping are medium; substitutions, ambiguity, or query-only insertions are
  low/failed. Only high-confidence, coverage-qualified positions are eligible
  for automatic conservation proposal bias, never automatic catalytic truth.
- The formal local scale is intentionally small and forced to CPU. GPU access
  recovered during the stage, but the local 8,188 MiB device was not used to
  expand the preregistered task. The result demonstrates a real, fair pipeline
  and provides descriptive estimates; it is not powered to declare a winner.
- All nine methods share the same hard/free position manifest. Requested labels
  do not override that comparability invariant.

## Discoveries

- The final-HEAD bundle contains five checkpoints, the 5.27 GB Atlas JSON, and
  eight experimental PDB/mmCIF assets in its external checksum manifest. All
  14 were present and matched; no large byte was embedded in the bundle.
- The real 6E9F mapping contains 61 terminal and 29 internal unresolved full
  positions. The 864 resolved positions all map exactly to original VI-D MSA
  columns; 712 additionally have MSA coverage at least 0.8 and are eligible
  for conservation proposal logic.
- The four deposited 6E9F catalytic sites are alanine at PDB positions
  295/300/849/854, but all map with high confidence. Every fair baseline
  restores and fixes their literature-supported biological R/H identities.
- GPU access recovered during this stage. An initial diagnostic MPNN process
  auto-selected the RTX 4060, so the required CPU checks were rerun with
  `CUDA_VISIBLE_DEVICES=''`; no large GPU task was launched.
- The first full test audit had 55/55 passing tests but only 65.90% coverage,
  below the 70% gate. Adding production-path fixture tests (without excluding
  modules or lowering the gate) yielded 59/59 tests and 72.03% coverage.
- The first two formal-run attempts stopped before proposal generation because
  the legacy ESM environment intentionally lacks pandas and nullable mapping
  indices were serialized as decimal integer strings. Both I/O defects were
  corrected without changing any sampling or endpoint rule, and their FAILED
  run records remain in `results/runs/`.
- The next endpoint-free identity audit exposed duplicate actual seeds across
  adjacent seed blocks and a consensus proposal just above the common identity
  interval. Before any conditional likelihood, novelty, or recovery endpoint
  was computed, the preregistration was amended to use disjoint proposal seeds
  and an identity-matched two-source consensus that never invents a residue
  outside the genuine ESM-IF1 and LigandMPNN source tokens.
- The first pushed stage-0002 CI exposed an environment parity gap hidden by
  the local Conda environment: pip `[dev]` omitted `pandas-stubs`, and current
  NumPy typing required dtype-parameterized arrays. This was reproduced in two
  clean pip environments, including Python 3.11, then corrected without any
  runtime algorithm or report change. The repaired GitHub run passed.

## Execution details

Work from `/home/junpeng/ESM-IF`. The immutable input sources are:

```text
data/experimental_structures/6e9f.cif
data/manifests/experimental_structures.yaml
data/manifests/cas13_functional_residues.yaml
data/processed/atlas/v1.0/msa/vi-d/alignment.fasta
data/processed/atlas/v1.0/conservation/vi-d.parquet
data/processed/atlas/v1.0/clusters/cas13_exact_unique.fasta
```

Use analysis, ESM-IF1, LigandMPNN, and bioinformatics tools only from their
project environments. Formal sampling uses local checkpoints with offline
loading. Initial CPU seeds are `20260731` and `20260732`; the larger exported
manifest uses a deterministic registered extension and is not marked run.

Historical production Atlas and pilot run directories are read-only inputs.
Commands that would rebuild them are out of scope for this plan.

## Validation and acceptance

Unit and integration tests must cover mapping gaps/substitutions/insertion
codes/terminal truncation, original-MSA preservation, confidence gating,
backend index translation, common mask equality, deterministic matching,
fixed-residue preservation, no-mock filtering, failure funnel completeness,
paired seed statistics, and report evidence labels.

Required final commands are:

```bash
make lint
make typecheck
make test
make smoke-cpu
make verify-reproducibility
```

The formal real run additionally fails unless:

- every selected row has `is_mock=false`;
- every method has the same selected count and seed set;
- every method's hard-fixed and free-position hashes are identical;
- every selected candidate has zero hard-fixed violations;
- the registered identity matching succeeds;
- all Atlas novelty parameters have one shared hash;
- reports distinguish Level 1 novelty from Level 2 model compatibility and
  contain no functional-validity claim.

## Artifacts and provenance

Authoritative outputs are planned at:

```text
reports/vi_d_mapping/
reports/matched_baselines/methods_table.csv
reports/matched_baselines/candidates.jsonl
reports/matched_baselines/candidate_funnel.csv
reports/matched_baselines/matched_statistics.csv
reports/matched_baselines/per_region_metrics.csv
reports/matched_baselines/failure_analysis.csv
reports/matched_baselines/report.md
reports/matched_baselines/report.html
reports/matched_baselines/figures/
results/runs/<matched-run-id>/
artifacts/bundles/<final-head-bundle>/
```

Every real run records resolved config, command, environment, hardware, git,
input/output manifests, seeds, stdout/stderr, failures, exit code, and SUCCESS
or FAILED. Generated reports are canonical summaries; immutable detailed run
data remain under `results/runs/`.

## Blockers and recovery

Local GPU access recovered, but the 8,188 MiB device is not accepted as the
large-scale target and no large GPU task was launched. Paired-repeat DCA remains
blocked by absent trustworthy Atlas repeat orientation and is not part of this
stage. If a larger model matrix is needed, execute only the exported one-command
tmux job on a suitable GPU node; do not insert mocks or loosen the registered
comparison post hoc.

## Outcomes and retrospective

The strict mapping gate passed for 864/954 (90.566%) four-layer positions; the
remaining 90 full-sequence positions are coordinate-unresolved (61 terminal,
29 internal) but retain full-to-MSA mappings. Conservation proposal bias was
released only for 712 resolved, high-confidence, coverage-qualified positions.

The final CPU matrix selected 18 genuine candidates (nine methods, two seed
blocks), with one common hard-position hash and one common free-position hash.
All 18 had zero fixed-position violations and real ESM-IF1 compatibility
scores; 14 passed every Level-1 novelty/QC gate. Three failures were
fail-closed missing Atlas coverage hits and one was a low-complexity failure.
All 18 unselected but identity-eligible oversampled proposals remain in the
failure table. With only two seed blocks, all statistical results are labelled
low-power descriptive; no adjusted exploratory p-value was below 1.0.

Two initial runs failed before sampling because the legacy ESM environment did
not include pandas and because nullable integer CSV fields used decimal text.
A later run failed before endpoint scoring when the consensus missed the common
identity interval and exposed overlapping actual seeds. All failures are
retained. A successful pre-handoff run and the final run reproduced all six
core scientific tables byte-for-byte; the final rerun corrected only GPU
handoff metadata. The larger ten-seed GPU extension remains explicitly
`not_run` and is not part of the reported result.
