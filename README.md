# ESM-IF-on-Cas13

ESM-IF-on-Cas13 is a reproducible research framework for Cas13 sequence design
on existing backbones using structural, evolutionary, and RNA-interface
constraints. It benchmarks ESM-IF1, constrained ESM-IF1, ProteinMPNN,
LigandMPNN, MSA-profile sampling, and matched random mutation under comparable
design masks and novelty ranges.

## Evidence boundary

The repository labels evidence from Level 0 (verified code and I/O) through
Level 4 (wet-lab validation). Generated sequences supported only by statistical,
inverse-folding, or structure-prediction evidence are computational candidates,
not validated or effective Cas13 proteins. This project does not contain
wet-lab protocols.

## Quick start

```bash
make bootstrap
conda run -p .tools/envs/analysis cas13-if --help
make lint
make typecheck
make test
make smoke-cpu
```

Create and lock the optional model/bioinformatics environments, then run the
genuine local smokes:

```bash
make bootstrap-specialized
make smoke-pyg
make smoke-esm-if1
make smoke-proteinmpnn
make smoke-ligandmpnn
```

Large data and model assets are explicit fetch steps and are never downloaded
by import or inference:

```bash
make fetch-references
make fetch-third-party
make fetch-models
make fetch-structures
make fetch-atlas
```

After the Atlas fetch has completed its size/hash checks:

```bash
make process-atlas
make cluster
make msa
make conservation
make coevolution-smoke
make candidate-novelty
make report
```

`coevolution-smoke` is a fixture-only MI/APC implementation check. The real
Atlas v1.0 paired analysis is currently data-blocked because repeat orientation
is unavailable; it is not reported as DCA. `candidate-novelty` uses the
authoritative clean 72-candidate pilot path declared in
`configs/candidate_filtering.yaml` and fails closed when no Atlas hit reaches
the registered query-coverage threshold.

Production workflows are defined in Snakemake and exposed through `make`.
Machine-specific paths belong in `configs/local.yaml`, which is ignored by Git;
copy `configs/paths.example.yaml` when adapting the project elsewhere.

The completed VI-D mapping and preregistered real small matched matrix are
reproduced without rebuilding Atlas or the historical pilot:

```bash
make map-vi-d
make matched-baselines CONFIG=configs/matched_baselines.yaml
```

The second command requires a clean worktree and local genuine checkpoints; it
rejects mock candidates and refuses to overwrite an existing canonical report.

## Reproducible runs

Every CLI execution writes an immutable run under
`results/runs/<date>-<experiment>-<config-hash>-<git-sha>/` with resolved
configuration, command, environment, hardware, git metadata, input/output
manifests, metrics, failures, logs, exit code, and a `SUCCESS` or `FAILED`
marker. See `docs/REPRODUCIBILITY.md`.

## Project status

Current real, fixture, mock, failed, and not-run work is reported in
`docs/STATUS.md`. Research decisions are recorded in `docs/DECISIONS.md`, and
the completed current-stage living plan is
`docs/execplans/0002_vi_d_mapping_and_matched_baselines.md`.
The network-interruption recovery audit and GPU/local work split are in
`docs/PHASE_SUMMARY_2026-07-31.md`; executable migration instructions are in
`docs/GPU_MIGRATION.md`. Stage-0002 results and failures are summarized in
`docs/STAGE_0002_SUMMARY_2026-07-31.md`.

## License

Project-authored source code is MIT licensed. Data, papers, model weights, and
third-party projects retain their own licenses; consult their manifests before
use or redistribution.
