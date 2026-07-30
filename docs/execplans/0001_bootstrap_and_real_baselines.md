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
- [ ] Complete and validate M1.
- [ ] Complete and validate M2.
- [ ] Complete and validate M3.
- [ ] Complete feasible M4/M5 work.
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

Current M1 tool gaps are recorded rather than treated as blockers for M0:
micromamba, `gh`, Git LFS, MAFFT, HMMER, Infernal, seqkit, Foldseek, and
US-align/TM-align were unavailable at initial audit. Conda/container bootstrap
routes exist, and independent work continues. GitHub push authentication remains
unverified because `gh` is absent; it will be tested at push time.

## Outcomes and retrospective

M0 outcome: the repository, CLI surface, immutable provenance, fixture
workflow, core Atlas/alignment/structure/constraint/refold algorithms, tests,
research governance, manifests, and CI are implemented and locally validated.
The machine audit is real. The CPU preflight is mock and supports Level 0 only.
No genuine model inference, experimental-structure benchmark, real Atlas parse,
or scientific candidate result has yet been produced.
