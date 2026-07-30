# Skopintsev et al. 2026 — SynTnpB

- **Research question:** Can structure- and evolution-constrained generation
  redesign minimal RNA-guided nucleases beyond natural sequence variation?
- **Data:** TnpB homologs, structures, evolutionary annotations, generated
  candidates, and experimental characterization.
- **Method:** ESM-IF1 generation with evolutionary/functional position
  constraints; implementation details are inspected in the pinned repository.
- **Loss:** Uses the pretrained inverse-folding objective; this project does not
  infer an additional loss without code/paper evidence.
- **Evaluation:** Computational selection and the authors' experimental and
  structural validation.
- **Reproducible resources:** MIT `pyskop/SynTnpBs` repository; publisher paper
  metadata. Local access to the full paper is not assumed.
- **Direct use here:** Reference for `partial_seq` semantics and regression
  behavior on a small fixture.
- **Not directly transferable:** TnpB architecture/function differs from Cas13,
  and reported wet-lab success does not transfer to our candidates.
- **Key risks:** Mistaking post-hoc replacement for causal conditioning,
  future-fixed-token visibility, notebook/environment drift, and data leakage.

