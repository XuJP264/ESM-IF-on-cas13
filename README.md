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

Large data and model assets are explicit fetch steps and are never downloaded
by import or inference:

```bash
make fetch-references
make fetch-third-party
make fetch-models
make fetch-structures
make fetch-atlas
```

Production workflows are defined in Snakemake and exposed through `make`.
Machine-specific paths belong in `configs/local.yaml`, which is ignored by Git;
copy `configs/paths.example.yaml` when adapting the project elsewhere.

## Reproducible runs

Every CLI execution writes an immutable run under
`results/runs/<date>-<experiment>-<config-hash>-<git-sha>/` with resolved
configuration, command, environment, hardware, git metadata, input/output
manifests, metrics, failures, logs, exit code, and a `SUCCESS` or `FAILED`
marker. See `docs/REPRODUCIBILITY.md`.

## Project status

Current real, fixture, mock, failed, and not-run work is reported in
`docs/STATUS.md`. Research decisions are recorded in `docs/DECISIONS.md`, and
the active living plan is
`docs/execplans/0001_bootstrap_and_real_baselines.md`.

## License

Project-authored source code is MIT licensed. Data, papers, model weights, and
third-party projects retain their own licenses; consult their manifests before
use or redistribution.

