# ExecPlan 0001: Bootstrap and Real Baselines

## Purpose and evidence boundary

Build a publication-grade, reproducible ESM-IF-on-Cas13 repository and execute
the first real data/model baselines. The plan covers Milestones 0–3 and the
locally feasible portions of Milestones 4–5. Repository and fixture validation
can establish Level 0. Novelty analyses can establish Level 1 for generated
candidates. Inverse-folding scores can establish Level 2. No output from this
plan is Level 4, and no computational candidate will be called an effective or
validated Cas13.

## Current state

At plan creation on 2026-07-31 (Asia/Shanghai), the repository has an empty
`main` branch, no commits, and an `origin` remote pointing to
`https://github.com/XuJP264/ESM-IF-on-cas13.git`. Machine capacity, network
access, GitHub authentication, model availability, and upstream data sizes
have not yet been measured.

## Milestones

### M0 — Repository and research governance

Create the maintained project tree, research documents, configuration schemas,
CLI, fixture workflow, CI, and packaging. Acceptance requires installable
package metadata, working command help, lint, type checking, unit tests, mock
integration, and explicit evidence-level language.

### M1 — Machine, environments, dependencies, and model assets

Capture the requested hardware/software audit; establish four isolated
environment specifications and reproducible bootstrap paths; pin and license
third-party dependencies; download and hash locally feasible model assets;
execute genuine ESM-IF1, ProteinMPNN, and LigandMPNN smokes where resources
permit. Missing tools or assets are accepted only as explicit reproducible
blockers, not as successful results.

### M2 — Experimental-structure benchmark

Use official RCSB data for 6E9F, 5XWP, and discovered related Cas13 structures.
Build a hashed manifest, protein/RNA chain classification, coordinate/sequence
mapping, strict QC, contacts, functional-region annotations, unified backend
schemas, real inverse-folding baselines, and a report. Acceptance requires at
least one real toy and one real Cas13 ESM-IF1 inference, unless a documented
external blocker persists after reproducible installation attempts.

### M3 — CRISPR-Cas Atlas dataset

Inspect and download the official v1.0 source after size and disk checks. Stream
records into auditable tables, retain ambiguous and failed pairings, deduplicate
sequences, cluster at all registered thresholds, and build leakage-safe split
audits and a data card. If the official asset is unavailable, preserve URL,
HTTP/error evidence, disk evidence, tested fixture paths, and the exact resume
command while continuing independent work.

### M4/M5 — Locally feasible evolutionary constraints and generation

Build subtype-specific MSA validation, weighted conservation/entropy,
protein–repeat MI/APC with bootstrap and permutation support, constraint
manifests, a correctly constrained autoregressive ESM-IF1 decoder, matched
baselines, sequence QC, novelty scoring, and candidate reports. Formal DCA and
large-scale generation remain `not_run` when hardware or data are insufficient;
their export/ingest contracts must still be exercised on fixtures.

### M6 — Migration and reproducibility handoff

Implement deterministic refold job export/ingest, manifest-based GPU bundles,
tmux launch scripts, mock prediction E2E, full acceptance command audit, staged
commits, and verified push.

## Progress

- [x] 2026-07-31: Inspected initial git state and remote.
- [x] 2026-07-31: Created repository instructions, ExecPlan standard, and this
  initial plan.
- [x] 2026-07-31: Ran and persisted the requested machine, GPU, software,
  container, authentication-tool, Git LFS, network, disk, memory, CPU, and
  bioinformatics-tool audit.
- [x] 2026-07-31: Completed and validated M0 with Ruff, strict mypy, 31 fixture
  tests, 79.67% branch-aware coverage, and a provenance-recorded CPU mock smoke.
- [x] 2026-07-31: Fetched and hashed four pinned upstream repositories, ESM-IF1
  and MPNN checkpoints, and four RCSB Cas13 structures.
- [x] 2026-07-31: Completed real Level 0 structure/QC manifests for 6E9F,
  5XWP, 6E9E, and 5XWY; separately audited deposited catalytic mutations.
- [x] 2026-07-31: Implemented and fixture-tested the offline ESM-IF1 backend,
  causal constrained decoder, Atlas Parquet pipeline, six-threshold clustering
  CLI, subtype MSA/conservation pipeline, SASA regions, benchmark runner, and
  evidence-aware project report.
- [x] 2026-07-31: Re-ran the current CPU quality gate after network
  interruption: Ruff/format pass for package, tests, and scripts; strict mypy
  pass for 42 source files; 42/42 tests; and 70.08% branch coverage.
- [ ] Restore current OS-level CUDA access and execute genuine GPU smoke tests.
  This requires recovery of WSL `/dev/dxg`; CPU work continues independently.
- [x] 2026-07-31: Repaired ESM user-site isolation and passed genuine CPU
  checkpoint load, toy score/sample, all-fixed zero-violation sampling, 6E9F
  score, and 5XWP score. GPU validation remains blocked by the current WSL
  `/dev/dxg` failure.
- [x] 2026-07-31: Passed a genuine 6E9F catalytic hard-fixed ESM-IF1 CPU
  sample with zero violations and a genuine ProteinMPNN CPU smoke with
  missing-coordinate slots explicitly masked and audited.
- [x] 2026-07-31: Repaired the LigandMPNN environment using its pinned upstream
  NumPy 1.23.5 and `dm-tree` requirements, then passed a genuine 6E9F CPU smoke
  with RNA B/C atomic context, fixed-residue preservation, statistics, and
  backbone outputs validated.
- [x] 2026-07-31: Built the fourth isolated bioinformatics environment, exported
  exact/conda/pip locks, and audited MMseqs2, MAFFT, HMMER, Infernal, seqkit,
  Foldseek, TM-align, and Git LFS with no missing executable.
- [x] 2026-07-31: Passed a genuine PyTorch Geometric GCNConv CPU smoke.
- [x] 2026-07-31: Completed the corrected 72-candidate ESM-IF1 pilot benchmark
  for 6E9F/5XWP across three temperatures and three constraint conditions;
  all candidate IDs are unique and all fixed-position violation counts are
  zero. The report remains Level 2 and does not claim a matched method win.
- [x] 2026-07-31: Exported a clean-commit GPU bundle and verified all internal
  files plus all currently available model/experimental-structure assets.
  The incomplete Atlas JSON is explicitly listed in `missing_assets`; target
  GPU-node validation remains pending.
- [x] 2026-07-31: Completed the resumed official Atlas download, verified the
  exact 5,267,508,328-byte size and SHA256
  `5b4ba2fb99638d279e0c126100e19a4b77aba487b37b7df118e4bf4acd494720`,
  atomically finalized the manifest, and passed an offline idempotence check.
- [ ] Complete GPU validation for M1; the local CPU/dependency portion passes.
- [ ] Complete the full ProteinMPNN/LigandMPNN/matched-novelty matrix for M2;
  the ESM pilot and individual real smokes pass.
- [x] 2026-07-31: Completed the core M3 production path: clean deterministic
  streaming parse, real funnel, exact dedup, all six registered MMseqs2
  thresholds, and strict 40% cluster leakage audit. Auxiliary subtype-held-out
  and scaffold-held-out splits remain a later M3 extension.
- [x] 2026-07-31: Completed feasible real M4 work on this node: length-gated
  subtype VI-B/D/F/I MAFFT alignments and coverage-gated weighted conservation.
  Scaffold mapping remains; paired MSA/MI/DCA are source-data blocked because
  Atlas exposes no trustworthy repeat orientation.
- [x] 2026-07-31: Completed an M5 real Level 1 candidate audit against the
  4,070-sequence Atlas resource. Fourteen of 72 pilot candidates passed all
  registered novelty/QC gates; low-complexity and no-hit failures remain in the
  audit and no candidate receives a functional claim.
- [x] 2026-07-31: Exported Atlas-complete clean source bundle
  `gpu-bundle-3e53026923aa-7540febfb2` from commit `3e53026923aa`; internal
  hashes plus five checkpoints, the 5.27 GB Atlas JSON, and eight experimental
  structure assets passed, with no missing assets.
- [ ] Complete M6 and the final acceptance audit.

## Decisions

- The primary package is `cas13_if` with a `cas13-if` Typer CLI.
- Snakemake is the formal research workflow; direct CLI and Make targets remain
  usable for diagnosis and small runs.
- Parquet is the canonical tabular output when `pyarrow` is available.
  Explicit fixture-only fallback formats may be used in CPU CI but may not be
  presented as the real Atlas result.
- All historical runs are immutable under `results/runs/<run_id>/`.
- Experimental structures are processed before waiting on predicted Atlas
  structures.

## Discoveries

- Initial repository: no commits and no tracked project files.
- Initial branch: `main`.
- Initial remote: target GitHub repository already configured.
- Measured GPU: NVIDIA GeForce RTX 4060 Laptop GPU with 8188 MiB VRAM;
  driver 560.94 and `nvidia-smi` driver-compatibility label CUDA 12.6. `nvcc`
  is absent.
- Measured host: 32 logical CPUs, 15 GiB RAM, about 737 GiB available on `/`.
- Docker 28.5.1 and Apptainer 1.4.5 are present. Conda 26.3.2 is present;
  micromamba, `gh`, functional Git LFS, and most requested bioinformatics tools
  were absent at initial audit.
- HTTPS checks for GitHub and RCSB succeeded. Git operations needed an
  out-of-sandbox network allowance because the default sandbox could not resolve
  GitHub.
- The first Conda attempt exposed two real bootstrap defects: unwritable global
  package caches and an editable path interpreted relative to `envs/`. The
  bootstrap now uses project-local caches and installs the project by absolute
  path after environment creation.
- Analysis environment validation: Ruff passed, strict mypy passed for 36 source
  files, and 31/31 tests passed with 79.67% coverage.
- The ESM checkpoint and all requested MPNN checkpoints were fetched and hashed.
- The first ESM Conda transaction reused `fair-esm` from the user site and
  failed the isolation gate. It was repaired by installing the pinned local ESM
  source with `PYTHONNOUSERSITE=1`; locks were regenerated and genuine CPU
  inference passed.
- PyTorch 2.4.1 with CUDA 12.1 was installed, but the current host reports both
  `torch.cuda.is_available()=false` and an OS-level `nvidia-smi` access block.
- Atlas v1.0 HEAD and a 2 MiB schema probe succeeded before the interruption.
  The subsequent full fetch failed with DNS resolution error before a
  production `.part` file was created.
- Current CPU validation after the new implementation work is 42/42 tests with
  70.08% branch-aware coverage.
- ProteinMPNN's 6E9F parser represents 864 resolved residues as 893 numbering
  slots, with 29 internal missing-coordinate slots masked as `X`. The initial
  smoke validator incorrectly required a length of 864; the corrected
  validation passed and preserves this distinction.
- All four project-local environments now exist with exact/conda/pip locks.
  The bioinformatics audit reports no missing executable and includes local
  Git LFS 3.7.1.
- The first full ESM pilot computed genuine outputs but generated ambiguous
  candidate IDs across temperature/constraint conditions. Historical outputs
  remain immutable; the corrected run uses condition-digest IDs and validates
  72/72 uniqueness before downstream use.
- A read-only scan of the first 100,000 production Atlas records showed that
  `summary.subtype` is often generic `VI` or empty while the Cas HMM contains
  the precise subtype. The parser now preserves raw/source/conflict fields and
  resolves explicit HMM subtypes before the full parse. Conflicting non-VI
  summaries are retained but cannot enter high-confidence pairing.
- The complete real Atlas parse yielded 1,246,088 operons, 12,353 Cas13
  records, 4,070 exact-unique sequences, 3,500 evolution-eligible exact-unique
  sequences, zero processing failures, zero high-confidence oriented pairs,
  and 11,727 ambiguous pairs.
- Six-threshold MMseqs2 clustering yielded 3,877/1,797/1,323/1,003/783/516
  clusters at 100/90/70/50/40/30%. The strict 40% split has
  3,335/160/575 train/validation/test sequences and passed the leakage gate.
- Atlas `truncated=00` did not establish full-length effectors. The inclusive
  MSA had 48–80 aa representatives and no 90%-coverage column. A preregistered
  700–1600 aa screen restored hundreds of high-coverage columns per subtype,
  but those columns still require explicit scaffold mapping.
- The real candidate novelty audit found 14/72 Level 1 pass rows. Fifty-three
  candidates had no Atlas hit at the required 80% query coverage and were
  failed closed; 41 failed low-complexity and 18 failed homopolymer gates.

## Execution details

Work from `/home/junpeng/ESM-IF`. Persist raw audit output in
`artifacts/system/software_initial.txt` and normalized values in
`artifacts/system/hardware.json`. Commands that may be unavailable are run
individually so one missing executable does not hide later results.

Use fixed fixture seeds (default `20260731`) and store production seeds in
resolved configs. Before every external download, record the source, expected
destination, available disk, HTTP metadata where possible, and checksum after
completion. Download into a temporary `.part` file and atomically rename only
after verification.

Do not let inference code silently access the network. Fetch scripts populate
local model/data manifests; runtime backends fail clearly if assets are absent.

## Validation and acceptance

The final audit runs, records, and truthfully reports:

```text
make lint
make typecheck
make test
make smoke-cpu
make smoke-esm-if1
make smoke-proteinmpnn
make smoke-ligandmpnn
make process-atlas
make cluster
make benchmark-experimental
make report
make export-gpu-bundle
make verify-reproducibility
```

Required behavioral checks include deterministic manifests; Atlas streaming and
ambiguous-pair routing; valid MSA enforcement; insertion-code and missing-atom
handling; RNA contact annotation; hard-fixed preservation with zero violations;
seed reproducibility; mock refold `is_mock=true`; leakage-stop behavior; and
immutable run directories.

## Artifacts and provenance

Planned authoritative paths:

- system audit: `artifacts/system/`
- data/model/third-party manifests: their respective `manifest.yaml` files
- run records: `results/runs/`
- benchmark outputs: `reports/runs/` and `reports/latest`
- status and decisions: `docs/STATUS.md`, `docs/DECISIONS.md`
- GPU bundle manifests: `artifacts/bundles/`

Paths become evidence only after the corresponding command and validation are
recorded here and in the run manifest.

M0 validation commands actually run from `/home/junpeng/ESM-IF`:

```text
bash scripts/bootstrap_local.sh
make lint
make typecheck
make test
make smoke-cpu
```

The smoke run is explicitly mock and Level 0:
`results/runs/20260731-benchmark-experimental-d2e8c34d52-nogit/`.

## Blockers and recovery

Initial M1 gaps included micromamba, `gh`, Git LFS, MAFFT, HMMER, Infernal,
seqkit, Foldseek, and US-align/TM-align. The isolated bioinformatics environment
now supplies every required bioinformatics executable, TM-align, and Git LFS.
Conda remains the selected environment manager, so micromamba is not required.
`gh` is absent; direct Git remote access works, and write authentication will be
tested at push time.

## Outcomes and retrospective

M0 outcome remains complete. Since the original retrospective was written,
genuine ESM-IF1 toy/6E9F/5XWP CPU inference, genuine constrained 6E9F sampling,
a genuine ProteinMPNN 6E9F CPU smoke, a genuine RNA-context LigandMPNN smoke,
and a corrected 72-candidate ESM pilot benchmark have been produced. The real
Atlas parse, exact dedup, six-threshold clustering, strict split, subtype MSA,
coverage-gated conservation, and candidate novelty audit are now also complete.
The full matched method matrix, scaffold-to-MSA mapping, real refold and target
GPU validation remain pending; paired-repeat analysis is blocked by the source
orientation field rather than compute. No result constitutes Level 4
functional validation.
