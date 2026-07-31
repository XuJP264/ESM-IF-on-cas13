# Project Status

Last updated: 2026-07-31 (Asia/Shanghai)

| Milestone | State | Evidence |
|---|---|---|
| M0 repository/governance | complete | Ruff/format pass for package, tests and scripts; strict mypy pass; pytest 42/42; 70.08% branch coverage; CPU mock smoke pass |
| M1 environments/assets | local CPU complete; GPU pending | four isolated environments and three lock forms each; all requested bioinformatics tools plus Git LFS available; genuine PyG, ESM-IF1, ProteinMPNN and LigandMPNN CPU smokes passed; GPU blocked by current WSL `/dev/dxg` failure |
| M2 experimental benchmark | ESM pilot complete; full matrix pending | four RCSB structures downloaded/QC passed; genuine 72-candidate ESM-IF1 temperature/constraint benchmark has 72 unique IDs and zero fixed violations; genuine ProteinMPNN and RNA-context LigandMPNN smokes passed; matched full method matrix pending |
| M3 Atlas | official download in progress | official 5,267,508,328-byte source HEAD/schema verified; initial DNS failure preserved; resumable production `.part` is now downloading |
| M4 evolution | implementation/fixture only | subtype MSA/conservation and MI/APC code tested; real Atlas-derived analysis pending |
| M5 constrained generation | real pilot, matrix pending | genuine causal decoder implemented; 6E9F catalytic hard-fixed CPU sample completed with zero fixed-position violations; baseline matrix/novelty report pending |
| M6 migration/refold | interface/fixture partial | provider-neutral mock E2E, bundle/sync/verify/tmux scripts and migration guide exist; bundle export and target-node validation pending |

Genuine ESM-IF1 CPU scoring now supplies Level 2 inverse-folding compatibility
for the toy, 6E9F, and 5XWP inputs. A genuine 6E9F catalytic hard-fixed sample
has zero fixed-position violations; ProteinMPNN and RNA-context LigandMPNN have
completed genuine 6E9F CPU smokes. These are not Level 4 functional validation.
The corrected ESM pilot benchmark is
`results/runs/20260731-benchmark-experimental-a998ff40aa-7599be0-r003`;
its raw recovery values are not a matched-design-position or matched-novelty
method comparison. No candidate yet has formal Level 1 Atlas novelty or Level 3
refold evidence.
Real Level 0 outputs also include the machine audit, model/third-party hashes,
RCSB downloads, structure QC, and functional-residue manifest. Fixture
candidate/refold/MSA paths remain mock or test-only. See
`docs/PHASE_SUMMARY_2026-07-31.md`.
