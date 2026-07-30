# Reproducibility

Each run is immutable and records configuration, command, environment, hardware,
git state, inputs, outputs, metrics, failures, standard streams, exit code, and
success/failure marker. Run IDs contain local date, experiment, canonical
configuration hash, and git short SHA. Existing run directories are never
overwritten.

Real experiments additionally record seed, checkpoint/file SHA256, package
lists, software versions, and timing. Mock and fixture runs are visibly marked.
Inference never silently downloads weights.

The project uses isolated Conda prefixes under `.tools/envs/`, versioned
environment specifications under `envs/`, and exact exported locks after
successful environment validation. `make verify-reproducibility` audits
manifests, required files, fixture tests, and the current GPU bundle.

