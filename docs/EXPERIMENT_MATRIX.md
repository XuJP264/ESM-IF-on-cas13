# Experiment Matrix

| Family | Method or ablation | Protein backbone | RNA atoms | Constraints |
|---|---|---:|---:|---|
| Reference | WT/native | yes | annotation | native |
| Baseline | matched random mutation | no model | no | matched mask/rate |
| Baseline | MSA profile sampling | no structure model | no | subtype profile |
| Baseline | ESM-IF1 | yes | no | none |
| Baseline | catalytic-only ESM-IF1 | yes | no | catalytic hard-fixed |
| Baseline | conservation ESM-IF1 | yes | no | conservation |
| Baseline | conservation + RNA ESM-IF1 | yes | contact annotation | conservation/contact |
| Baseline | conservation + RNA + coevolution ESM-IF1 | yes | annotation | combined |
| Baseline | ProteinMPNN | yes | no | matched fixed positions |
| Baseline | LigandMPNN | yes | yes | matched fixed positions |
| Consensus | ESM-IF1 / LigandMPNN | yes | model-specific | agreement filter |

Ablations remove conservation, direct RNA contacts, RNA second shell,
coevolution, or multi-state constraints; compare single/multiple protein chains;
and vary temperature, mask threshold, sequence weighting, subtype-specific
versus mixed MSA, and single versus multiple structural states.

