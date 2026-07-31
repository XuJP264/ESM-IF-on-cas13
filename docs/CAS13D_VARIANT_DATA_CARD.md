# Cas13d variant/activity data card

## Scope and evidence boundary

This dataset is a curated retrospective resource for testing whether structure-
and sequence-based scores track previously reported Cas13d phenotypes. It is
not a training set, does not establish causality, and does not convert a
computational candidate into an active Cas13. The curated observations are
experimental facts from their source studies; the processing and model-score
links in this repository are Level 0–2 evidence only.

Canonical outputs are
`data/processed/cas13d_variant_activity.parquet`,
`data/processed/cas13d_variant_activity_sources.csv`, and
`data/processed/cas13d_variant_label_sensitivity.csv`. Raw papers and public
supplements remain outside Git. The machine-readable curation source is
`data/manifests/cas13d_variant_activity_curated.yaml`.

## Version 1 contents

- 22 records from four studies and four non-poolable comparability groups.
- 20 distinct mutation strings; all 22 full mutant sequences were recovered
  after validating the one-based CasRx natural-sequence residue numbering.
- 16 records are the 2025 purified in-vitro cis/trans panel. Fifteen have
  numerical graph-read values and one purification failure is retained as a
  failure, not imputed as zero activity.
- Four 2018 human-cell deletion variants, one 2023 high-fidelity variant, and
  one 2024 high-precision variant are retained as qualitative, assay-specific
  observations.
- Primary labels at the preregistered thresholds are 8 active, 7 partial,
  5 inactive, 1 cis-retained/trans-reduced, and 1 not assayed.

## Sources and extraction

The sources table records DOI, official URL, access status, asset hash where a
public supplement was downloaded, and the extraction route. The 2025 Figure 6
values are manual approximate readings to 0.01 from the public figure; they are
marked `numeric_is_approximate=true`, with `n=3` metadata. The supplement ZIP
has SHA256
`bc5e73b1f9c4e864511c03b4b9ba84abc43360b3f6b66f5d700a16fffcf070a1`.
The 2023 public supplement has SHA256
`3fee0eb1071c91902c894e3767ad5076b1b86af38d35741e5b31ac00cc12fb7a`.
Copyrighted source assets are not redistributed.

## Harmonization

Mutation expressions support substitutions, inclusive deletions, compound
deletions, and explicit insertions. Every stated wild-type residue is checked
against the natural CasRx sequence before applying a substitution. Exact
records are deduplicated by study, scaffold, mutation, assay, guide, target,
and comparability group. Activities are WT-normalized only within the source
group; values from purified biochemical assays, reporter knockdown, collateral
assays, cell assays, and mouse studies are never pooled on a common scale.

The sensitivity table recomputes active/partial/inactive calls across several
documented cutoffs. Threshold sensitivity is descriptive and must not be used
to tune a retrospective classifier on the same observations.

## Missingness and limitations

- The dataset is small, deliberately heterogeneous, and biased toward variants
  selected by individual studies rather than systematic saturation mutagenesis.
- Several studies provide qualitative categories rather than recoverable
  numerical endpoints.
- Figure-derived numerical values have digitization uncertainty and are not raw
  replicate data.
- Deletions and insertions change sequence/backbone correspondence. They are
  excluded from same-backbone inverse-folding comparisons unless a matching
  experimental mutant structure exists.
- The CasRx structures used for score comparisons are nuclease-inactive
  experimental constructs; catalytic residues are restored in the biological
  sequence mapping and this distinction is audited.
- Scaffold, assay, guide, target, expression system, and phenotype definitions
  can confound score/activity associations. Any benchmark therefore reports
  descriptive estimates and confidence intervals without claiming biological
  validation.

## Intended retrospective use

Within a single comparability group, the benchmark may compare ESM-IF1 delta
log-likelihood, ProteinMPNN and LigandMPNN model scores, conservation, RNA
contact class, and multi-state ESM-IF1 score against active, partial, inactive,
or cis-retained/trans-reduced observations. The statistical unit is the
experimentally reported variant; repeated structure states are not independent
biological replicates. Results must remain explicitly retrospective and must
not be used as wet-lab readiness evidence.
