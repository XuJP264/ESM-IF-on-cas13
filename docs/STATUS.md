# Project Status

Last updated: 2026-07-31 (Asia/Shanghai)

| Milestone | State | Evidence and remaining gate |
|---|---|---|
| M0 repository/governance | complete | ExecPlan/governance/tree/CI present; Ruff, format, strict mypy and 47 tests pass |
| M1 environments/assets | local CPU complete; GPU pending | four isolated environments and three lock forms each; requested bioinformatics tools plus Git LFS available; genuine PyG/ESM-IF1/ProteinMPNN/LigandMPNN CPU smokes pass; current WSL lacks `/dev/dxg` |
| M2 experimental benchmark | ESM pilot complete; full matrix pending | four RCSB structures downloaded/QC; clean 72-candidate ESM-IF1 pilot has unique IDs and zero fixed violations; genuine MPNN CPU smokes pass; matched design-position/novelty matrix not run |
| M3 Atlas | core production path complete; auxiliary splits pending | official 5,267,508,328-byte source verified; 1,246,088 operons streamed with zero processing failures; 4,070 exact-unique Cas13; six MMseqs thresholds and strict 40% leakage gate pass; subtype/scaffold-held-out auxiliary splits remain |
| M4 evolution | real subtype MSA/conservation partial; paired analysis data-blocked | length-gated VI-B/D/F/I MSAs and coverage-gated conservation are real; scaffold mapping remains; Atlas repeat orientation is unavailable, yielding 0 high-confidence and 11,727 ambiguous pairs, so real MI/APC/DCA are blocked |
| M5 constrained generation | real pilot and Level 1 audit partial | causal constrained decoder and 6E9F real sample have zero violations; 14/72 pilot candidates pass registered Level 1 novelty/QC; full baselines, matched novelty and ablations remain |
| M6 migration/refold | Atlas-complete source bundle verified; target pending | provider-neutral mock E2E and migration scripts exist; clean bundle `gpu-bundle-3e53026923aa-7540febfb2` passed internal plus 14 large-asset hashes with no missing assets; target transfer/bootstrap/GPU/refold remain |

## Current scientific evidence

- Level 0: repository/I/O, structures, Atlas production parse, clustering, MSA
  and conservation pipelines have real or fixture evidence as individually
  labeled.
- Level 1: 14 of 72 ESM-IF1 pilot candidates pass all registered sequence
  novelty/QC gates. The other candidates are not promoted.
- Level 2: genuine ESM-IF1 score/sample and experimental-structure pilot;
  genuine ProteinMPNN and RNA-context LigandMPNN CPU smokes.
- Level 3: not available; no real candidate refold/multimodel structural
  validation.
- Level 4: not available and outside the current no-wet-lab scope.

No candidate is described as a validated or effective Cas13.

## Authoritative paths

- phase audit: `docs/PHASE_SUMMARY_2026-07-31.md`;
- clean Atlas parse:
  `results/runs/20260731-atlas-processing-e8356ef7b5-eebc1a5-r001/`;
- Atlas funnel: `data/processed/atlas/v1.0/data_funnel.json`;
- cluster summary:
  `data/processed/atlas/v1.0/clusters/clustering_summary.json`;
- MSA manifest: `data/processed/atlas/v1.0/msa/msa_manifest.json`;
- conservation manifest:
  `data/processed/atlas/v1.0/conservation/conservation_manifest.json`;
- clean ESM pilot:
  `results/runs/20260731-benchmark-experimental-a998ff40aa-ab6c9c5/`;
- real novelty audit:
  `results/runs/20260731-candidate-filtering-b14455d461-6d258de/`;
- latest real report at this update:
  `results/runs/20260731-benchmark-experimental-bcfd0be469-3ebd1c9/report/`.

## Current quality gate

- Ruff lint: pass;
- Ruff format: pass;
- strict mypy: pass for 43 source files;
- pytest: 47/47 pass;
- branch-aware coverage: 70.89% (required 70%);
- clean Atlas parse and cluster leakage gate: pass;
- real candidate novelty audit: pass with fail-closed missing-hit semantics;
- target GPU validation: not run.
