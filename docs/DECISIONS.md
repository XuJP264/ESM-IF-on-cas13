# Decision Log

## 2026-07-31 — Evidence terminology

Adopt Levels 0–4 globally. Levels 1–3 remain computational candidates; Level 4
requires wet-lab evidence outside the present scope.

## 2026-07-31 — Experimental structures first

Benchmark 6E9F and 5XWP before large-scale structure prediction. This supplies
real backbones and RNA-context annotations without waiting for AlphaFold.

## 2026-07-31 — Local asset loading

All formal model inference uses explicitly fetched local checkpoints. Runtime
network downloads are errors because they prevent reliable hashing and audit.

## 2026-07-31 — Split isolation

Use sequence-cluster-level splits with a strict 40% identity generalization
test. Auxiliary 50%/70% splits may be added if subtype or structure coverage is
insufficient, but cannot replace the strict test.

## 2026-07-31 — Measured local resource limits

Plan for a measured RTX 4060 Laptop GPU with 8188 MiB VRAM, not 32 GB. Use
single-structure/small batches, CPU fallback where algorithmically equivalent,
and GPU/HPC exports for jobs that do not fit.

## 2026-07-31 — Project-local Conda state

Use `.tools/pkgs`, `.tools/cache`, and `.tools/envs/*` so environment creation
does not depend on writable global Miniforge or user cache directories. Install
the project from its absolute repository path after Conda resolves the
environment; relative `-e .` entries inside `envs/*.yml` are ambiguous because
Conda resolves them relative to the environment-file directory.

## 2026-07-31 — M0 quality scope

Apply strict mypy to the `cas13_if` package and require at least 70% branch-aware
coverage for CPU fixture tests. M0 passed at 79.67%. Genuine model/network tests
remain explicitly marked and excluded from CPU CI.
