# Model Card

This project orchestrates independently licensed upstream models: ESM-IF1,
ProteinMPNN, and LigandMPNN. Checkpoints are fetched explicitly into `models/`,
hashed, and never committed. Runtime loading is offline and fails when a local
checkpoint is missing.

ESM-IF1 conditions on protein N/CA/C backbone coordinates and does not consume
RNA as a protein chain. RNA is used for contact annotations. LigandMPNN is the
RNA atomic-context baseline; ProteinMPNN is the matched protein-only graph
baseline. The local constrained ESM-IF1 decoder is left-to-right causal:
earlier positions cannot observe future fixed tokens. Post-hoc replacement is
not conditional generation.

Scientific limitations and genuine smoke-test status are tracked per checkpoint
in `models/manifest.yaml` and `docs/STATUS.md`.

