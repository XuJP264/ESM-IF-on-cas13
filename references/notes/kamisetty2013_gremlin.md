# Kamisetty et al. 2013 — GREMLIN

- **Research question:** When do global coevolution models yield useful
  residue-contact predictions?
- **Data:** Deep protein-family MSAs and known structures.
- **Method:** Pseudolikelihood-based global Potts/Markov random-field model with
  corrected coupling scores and optional structural priors.
- **Loss:** Regularized pseudolikelihood for sequence-family observations.
- **Evaluation:** Top-L/L/2 contact precision and dependence on alignment depth.
- **Reproducible resources:** Public author copy and historical GREMLIN service.
- **Direct use here:** One acceptable formal direct-coupling family and
  structural validation design.
- **Not directly transferable:** Long protein–repeat concatenations may lack
  effective depth; protein–RNA alphabet/regularization need explicit handling.
- **Key risks:** Computational scale, paralog/pairing errors, phylogeny,
  overfitting, and substituting MI when formal DCA is not run.

