# Preregistered Experiment Protocol

## Scope and evidence

Experiments use fixed Cas13 backbones and compare computational sequence-design
methods. Levels 1–3 describe novelty, inverse-folding compatibility, and
multi-model structural plausibility. No computational result is Level 4.

## Hypotheses

**H1.** At equal parent identity or equal designed-position fraction,
evolution-constrained ESM-IF1 will retain HEPN critical sites, RNA interfaces,
highly conserved natural positions, and multi-state compatibility better than
unconstrained ESM-IF1.

**H2.** LigandMPNN will outperform ProteinMPNN and ESM-IF1 at RNA-contact and
RNA second-shell residue recovery and conditional score because it observes RNA
atomic context.

**H3.** Combining conservation, direct RNA contact, and validated
direct-repeat coevolution will outperform any single constraint without wholly
sacrificing sequence novelty.

**H4.** Multi-state constrained designs will yield fewer candidates compatible
with only one conformation than designs based solely on an apo or bound state.

**Future H5.** Cas13 domain adaptation of ESM-IF1 should outperform the original
model only when zero-shot benchmarks show a clear systematic defect. Domain
adaptation will not be performed merely to add a fine-tuning experiment.

## Units and splits

Independent units for inference are scaffolds, subtypes, sequence clusters, and
structural states. Thousands of candidates from one scaffold are not thousands
of independent protein experiments. The strict sequence test split is isolated
at 40% MMseqs2 identity; 50% and 70% splits may be auxiliary, never replacements.
Subtype-held-out and scaffold-held-out analyses are separate. Thresholds and
temperature choices use training/validation only.

## Comparability

Methods share scaffolds, designable positions, fixed residues, random-seed
registry, and matched novelty bins. ProteinMPNN and LigandMPNN receive the same
protein backbone; only the latter receives retained RNA atomic context.

## Statistics

Report scaffold-aware bootstrap confidence intervals, paired effects, effect
sizes, multiplicity-adjusted comparisons, permutation nulls, matched-novelty
analyses, sensitivity analyses, and seed variability. Keep all failed runs in
the audit.

