# Project Status

Last updated: 2026-07-31 (Asia/Shanghai)

| Milestone | State | Evidence and remaining gate |
|---|---|---|
| M0 repository/governance | complete | ExecPlan/governance/tree/CI present; Ruff, format, strict mypy and 47 tests pass |
| M1 environments/assets | local CPU complete; GPU acceptance pending | four isolated environments and three lock forms each; requested bioinformatics tools plus Git LFS available; genuine PyG/ESM-IF1/ProteinMPNN/LigandMPNN CPU smokes pass; GPU access recovered during stage 0002 but the requested large GPU work remains intentionally not run |
| M2 experimental benchmark | real small fair matrix complete; scale pending | four RCSB structures downloaded/QC; clean 72-candidate ESM-IF1 pilot remains characterization only; genuine 6E9F/5XWP ProteinMPNN and RNA-context LigandMPNN smokes plus protein/soluble checkpoint smokes pass; the 9-method × 2-seed real VI-D matrix completed with 18 Level-2-scored candidates, no mocks, and zero fixed violations |
| M3 Atlas | core production path complete; auxiliary splits pending | official 5,267,508,328-byte source verified; 1,246,088 operons streamed with zero processing failures; 4,070 exact-unique Cas13; six MMseqs thresholds and strict 40% leakage gate pass; subtype/scaffold-held-out auxiliary splits remain |
| M4 evolution | VI-D mapping complete; paired analysis data-blocked | length-gated VI-B/D/F/I MSAs and coverage-gated conservation are real; 6E9F has 864/954 exact four-layer mapped positions and 90 coordinate-unresolved positions, with 712 resolved coverage-qualified conservation positions; Atlas repeat orientation is unavailable, yielding 0 high-confidence and 11,727 ambiguous pairs, so real MI/APC/DCA are blocked |
| M5 constrained generation | preregistered small fair matrix complete; scale/ablations pending | causal constrained decoder and all 18 selected VI-D matrix candidates have zero fixed violations; 14/18 matrix candidates pass registered Level 1 novelty/QC; larger seed count, multistate/scaffold replication, coevolution arm, and full ablations remain |
| M6 migration/refold | final-HEAD source bundle verified; target pending | provider-neutral mock E2E and migration scripts exist; clean bundle `gpu-bundle-a9a530d14434-7540febfb2` pins final handoff HEAD and passed internal plus 14 large-asset hashes with no missing assets; target transfer/bootstrap/GPU/refold remain |

## Current scientific evidence

- Level 0: repository/I/O, structures, Atlas production parse, clustering, MSA
  and conservation pipelines have real or fixture evidence as individually
  labeled.
- Level 1: 14 of 18 selected matched-matrix candidates pass every registered
  sequence novelty/QC gate. Historical pilot results remain separately labelled.
- Level 2: all 18 selected matched-matrix candidates received genuine ESM-IF1
  compatibility scores on the 6E9F backbone; genuine ProteinMPNN and
  RNA-context LigandMPNN CPU sampling also completed.
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
- strict VI-D mapping: `reports/vi_d_mapping/`;
- final matched run:
  `results/runs/20260731-vi-d-matched-baselines-3e2655746e-f5932f7/`;
- matched report: `reports/matched_baselines/report.md` and `report.html`;
- stage-0002 audit: `docs/STAGE_0002_SUMMARY_2026-07-31.md`.

## Current quality gate

- Ruff lint: pass;
- Ruff format: pass;
- strict mypy: pass for 47 source files;
- pytest: 61/61 pass;
- branch-aware coverage: 72.17% (required 70%);
- matched-report acceptance: pass for 18 candidates, nine methods, two common
  seed blocks, one fixed/free mask, no mocks, and zero fixed violations;
- final local commands `make lint`, `make typecheck`, `make test`,
  `make smoke-cpu`, and `make verify-reproducibility`: pass at `acd1137`;
- clean Atlas parse and cluster leakage gate: pass;
- real candidate novelty audit: pass with fail-closed missing-hit semantics;
- GitHub Actions run `30633893318`: primary fixture/code validation success;
- GitHub Actions run `30635299503`: preceding-stage final handoff/documentation
  validation success; it is not the primary code-run identifier;
- tag `v0.1.0-data-pipeline` points to final data-pipeline HEAD
  `a9a530d14434e74dc0cfc47896847e201431c1c2`;
- final-HEAD GPU bundle: `gpu-bundle-a9a530d14434-7540febfb2`, export
  `dirty=false`, `missing_assets=[]`, internal hashes and 14/14 source-asset
  hashes passed;
- target GPU validation: not run.
