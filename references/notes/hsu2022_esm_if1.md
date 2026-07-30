# Hsu et al. 2022 — ESM-IF1

- **Research question:** Can predicted structures expand inverse-folding
  training and improve recovery on structurally held-out backbones?
- **Data:** Experimental structures plus approximately 12 million
  AlphaFold2-predicted structures with structural holdouts.
- **Method:** Invariant geometric input processing followed by an
  autoregressive sequence model conditioned on protein backbone coordinates.
- **Loss:** Autoregressive token cross-entropy / conditional log likelihood.
- **Evaluation:** Native recovery, perplexity, buried/surface recovery,
  structural holdouts, complexes, interfaces, and multi-state examples.
- **Reproducible resources:** PMLR paper and `facebookresearch/esm`.
- **Direct use here:** Genuine scoring/sampling, temperature baselines,
  per-position probabilities, and protein multi-chain conditioning.
- **Not directly transferable:** RNA is not a standard protein-chain input;
  generic training does not establish Cas13 functional competence.
- **Key risks:** Silent checkpoint download, chain/coordinate mapping,
  autoregressive fixed-token semantics, and training-data overlap.

