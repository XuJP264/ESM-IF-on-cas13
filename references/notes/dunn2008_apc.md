# Dunn et al. 2008 — APC-corrected mutual information

- **Research question:** Can background mutual information caused by entropy and
  phylogeny be reduced for residue-contact prediction?
- **Data:** Protein-family multiple sequence alignments and structural contacts.
- **Method:** Mutual information followed by average product correction (APC).
- **Loss:** Not applicable; this is a statistical estimator.
- **Evaluation:** Residue-contact prediction relative to uncorrected MI.
- **Reproducible resources:** Publisher metadata and formulas.
- **Direct use here:** Required MI/APC baseline on paired protein–repeat MSAs.
- **Not directly transferable:** APC does not separate direct from indirect
  coupling and does not remove all phylogenetic bias.
- **Key risks:** Gap treatment, undersampling, inconsistent sequence weighting,
  and calling APC-MI a DCA result.

