# Høie et al. 2025 — AntiFold

- **Research question:** Does antibody-specific ESM-IF1 fine-tuning improve
  inverse folding and mutation scoring in antibody regions?
- **Data:** Solved antibody complexes and large predicted antibody-structure
  collections with identity-aware splits.
- **Method:** ESM-IF1 fine-tuning with masking strategies, region weighting, and
  layer-wise learning-rate choices.
- **Loss:** Amino-acid reconstruction cross-entropy under the selected masking
  scheme.
- **Evaluation:** Region recovery, refolded structural agreement, and zero-shot
  binding-affinity ranking.
- **Reproducible resources:** Open article, BSD-3-Clause package and code.
- **Direct use here:** Domain-adaptation precedent and leakage-aware evaluation.
- **Not directly transferable:** Antibodies have different architecture,
  alignments, functions, data volume, and test tasks than Cas13.
- **Key risks:** Fine-tuning without a preregistered zero-shot defect, predicted
  structure leakage, and optimizing recovery rather than functional regions.

