# Ruffolo et al. 2025 — OpenCRISPR / CRISPR-Cas Atlas

- **Research question:** Can large-scale CRISPR-Cas mining and language-model
  generation yield diverse genome-editor candidates?
- **Data:** Approximately 1.25 million mined CRISPR-Cas operons from assembled
  genomes and metagenomes; official counts differ slightly by release text.
- **Method:** Atlas mining/annotation followed by protein-language-model
  training and candidate prioritization.
- **Loss:** Language-model next-token objective plus downstream model-specific
  objectives.
- **Evaluation:** Sequence diversity, biochemical/cellular assays and editor
  characterization reported by the authors.
- **Reproducible resources:** Atlas GitHub, v1.0 JSON, OpenCRISPR repository.
- **Direct use here:** Natural Cas13 sequences, subtypes, operons, arrays, and
  direct repeats.
- **Not directly transferable:** A record is not an independent sequence;
  spacers are not verified targets; ambiguous arrays/effectors need routing.
- **Key risks:** Very large JSON, CC BY-NC terms, orientation ambiguity,
  annotation error, duplication, and random-record leakage.

