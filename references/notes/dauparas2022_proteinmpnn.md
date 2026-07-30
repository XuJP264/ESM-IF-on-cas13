# Dauparas et al. 2022 — ProteinMPNN

- **Research question:** Can message passing over protein backbones improve
  sequence design robustness and recovery?
- **Data:** Experimentally determined protein structures with sequence/structure
  cluster-aware splits.
- **Method:** Order-agnostic autoregressive message-passing network on protein
  backbone graphs.
- **Loss:** Categorical cross-entropy on amino-acid identity.
- **Evaluation:** Native recovery, soluble expression, structural accuracy, and
  experimentally tested designs.
- **Reproducible resources:** Public paper/manuscript, MIT GitHub code, weights.
- **Direct use here:** Protein-only graph baseline using the same backbone and
  design positions as LigandMPNN.
- **Not directly transferable:** Does not observe RNA atoms; experimental
  successes on other proteins do not validate Cas13 candidates.
- **Key risks:** Chain JSON semantics, fixed-position enforcement, checkpoint
  selection, and incomparable raw score scales.

