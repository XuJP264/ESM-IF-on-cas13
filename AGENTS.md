# ESM-IF-on-Cas13 Repository Instructions

These instructions apply to the entire repository.

## Research scope and claims

This project develops a Cas13 inverse-folding framework combining structural,
evolutionary, and RNA-interface constraints. It preserves experimentally or
computationally supplied Cas13 backbones and does not claim to create new
backbones.

Every result and document must use the following evidence levels:

- Level 0: code paths and input/output behavior are verified.
- Level 1: generated sequences have measured statistical novelty.
- Level 2: an inverse-folding model supports compatibility with the target
  backbone.
- Level 3: multiple models and structure prediction support functional
  plausibility.
- Level 4: wet-lab validation supports calling a sequence an effective Cas13.

Levels 1–3 are computational candidates, never "validated" or "effective"
Cas13 proteins. This repository does not design or provide wet-lab protocols
for cloning, culture, expression, delivery, or experimental operation.

## Planning and milestone discipline

- Use an ExecPlan for every major feature or research milestone.
- Follow `.agent/PLANS.md`; plans are living records and must be updated while
  work proceeds.
- After every milestone, update the relevant ExecPlan, `docs/STATUS.md`, and
  `docs/DECISIONS.md`.
- Prefer autonomous problem solving and continue independent work when a
  recoverable problem blocks one path.
- Keep `README.md` and documented commands runnable throughout development.
- Do not stop at describing what could be done next when safe, in-scope
  implementation or verification remains possible.

## Evidence, provenance, and failures

- Never fabricate results, metrics, downloads, tool availability, or
  successful runs.
- Mark mock data and mock results explicitly with `is_mock: true` in
  machine-readable output and with a visible MOCK label in reports.
- Mark unavailable or deferred analyses as `not_run` and record the exact
  blocker. Never substitute an approximate method under the name of a formal
  algorithm.
- Every real experiment records resolved configuration, random seed, software
  versions, git commit, hardware, input hashes, output hashes, logs, exit code,
  and failures.
- A failed data-leakage check stops the affected experiment. Test-set data may
  not be used for threshold choice, masking, model selection, or temperature
  tuning.
- Failed tests mean the milestone is not complete. All failed runs remain in
  the audit trail.

## Data, models, references, and third-party code

- Do not commit raw large datasets, model weights, paper PDFs, credentials,
  secrets, tokens, or private data.
- Keep manifests, hashes, licenses, download sources, and `.gitkeep`/README
  placeholders under version control.
- Project-authored code is MIT licensed. Third-party code, data, papers, and
  weights retain their original licenses and must not be relicensed as MIT.
- Pin third-party repositories to commits or tags. Do not edit vendored or
  submodule source directly. Put adapters in `src/cas13_if/backends/`; put any
  unavoidable upstream changes in standalone patch files and record them.
- Only download papers from lawful, clearly open sources. Record metadata for
  paywalled papers without bypassing access controls.
- Formal inference uses local checkpoints. Missing checkpoints must cause an
  explicit failure; runtime code must not silently download weights.

## Engineering quality

- Keep analysis, ESM-IF1, LigandMPNN, and bioinformatics environments isolated.
  Do not install into system Python.
- CPU CI must not download large models or datasets.
- Use deterministic identifiers and manifests. Do not overwrite historical
  runs.
- Treat indexing conventions, chain mappings, insertion codes, missing atoms,
  RNA/protein classification, and fixed-residue enforcement as tested
  correctness requirements.
- Do not treat RNA as an ESM-IF1 protein input; RNA supplies contact annotation
  or atomic context for compatible models.
- Statistical inference uses scaffold, subtype, sequence cluster, or structure
  state as the independent unit. Candidate-level distributions are descriptive
  and must not create pseudoreplication.

## Git workflow

- Make a small, explicit commit after each verified milestone or coherent
  verified feature.
- Do not commit a milestone with failing required tests.
- Never force-push.
- Push only after the relevant local validation passes.
- Preserve unrelated user changes and inspect the worktree before committing.

