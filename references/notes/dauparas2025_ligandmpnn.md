# Dauparas et al. 2025 — LigandMPNN

- **Research question:** Can explicit nonprotein atomic context improve sequence
  design at ligand, nucleotide, and metal interfaces?
- **Data:** PDB assemblies containing protein and nonprotein atomic contexts,
  clustered for held-out tests.
- **Method:** Protein graph, intraligand graph, protein–ligand graph, and an
  autoregressive sequence/side-chain decoder.
- **Loss:** Categorical sequence cross-entropy; side-chain objectives are
  described in the upstream work.
- **Evaluation:** Native recovery near small molecules, nucleotides and metals,
  side-chain conformations, and experimental designs.
- **Reproducible resources:** Public article, MIT GitHub code and weights.
- **Direct use here:** RNA-context baseline for Cas13–crRNA–target complexes.
- **Not directly transferable:** Reported nucleotide averages are not Cas13
  functionality and raw scores are not calibrated against ESM-IF1.
- **Key risks:** RNA atoms being deleted/misclassified, wrong design chains,
  fixed residue leakage, and packing mode changes.

