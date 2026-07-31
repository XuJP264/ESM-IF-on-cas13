# Project Status

Last updated: 2026-08-01 (Asia/Shanghai)

| Milestone | State | Evidence and remaining gate |
|---|---|---|
| M0 repository/governance | complete | ExecPlan/governance/tree/CI present; Ruff, format, strict mypy and 47 tests pass |
| M1 environments/assets | local CPU complete; GPU acceptance pending | four isolated environments and three lock forms each; requested bioinformatics tools plus Git LFS available; genuine PyG/ESM-IF1/ProteinMPNN/LigandMPNN CPU smokes pass; GPU access recovered during stage 0002 but the requested large GPU work remains intentionally not run |
| M2 experimental benchmark | real multi-scaffold preparation complete; Level-3 pending | Stage 0003A adds four independent Cas13d parents and nine experimental state units, all with strict four-layer mapping/RNA contacts; genuine ESM-IF1 scored all states and a 48-row four-scaffold × six-method × two-seed real smoke completed with no runtime failures, mocks, or fixed violations |
| M3 Atlas | core production path complete; auxiliary splits pending | official 5,267,508,328-byte source verified; 1,246,088 operons streamed with zero processing failures; 4,070 exact-unique Cas13; six MMseqs thresholds and strict 40% leakage gate pass; subtype/scaffold-held-out auxiliary splits remain |
| M4 evolution | VI-D mapping complete; paired analysis data-blocked | length-gated VI-B/D/F/I MSAs and coverage-gated conservation are real; 6E9F has 864/954 exact four-layer mapped positions and 90 coordinate-unresolved positions, with 712 resolved coverage-qualified conservation positions; Atlas repeat orientation is unavailable, yielding 0 high-confidence and 11,727 ambiguous pairs, so real MI/APC/DCA are blocked |
| M5 constrained generation | multi-scaffold smoke and pre-wet-lab protocol complete; formal selection pending | causal masks/generation now cover Es, Ur, Dj, and CasRx representatives; the frozen selection protocol requires parent-aware Level-1 QC, two real inverse-folding model families, multi-state gates, novelty strata, and cluster diversity; no final shortlist exists |
| M6 migration/refold | Stage-0003 job/ingest interface and clean-bundle procedure complete; real target run pending | 70 proteins expand to 1,068 deterministic monomer/binary/ternary jobs; four labeled mock outputs pass ingest, genuine TM-align invocation, structural metrics, consistency, retry, and Pareto E2E; a clean execution-HEAD bundle passed 1,204 internal-file and 28 large-asset hashes; real prediction count is zero |

## Current scientific evidence

- Level 0: repository/I/O, structures, Atlas production parse, clustering, MSA
  and conservation pipelines have real or fixture evidence as individually
  labeled.
- Level 1: 14 of 18 selected matched-matrix candidates pass every registered
  sequence novelty/QC gate. Historical pilot results remain separately labelled.
- Level 2: all 18 selected matched-matrix candidates and 48 Stage-0003A
  multi-scaffold smoke candidates received genuine ESM-IF1 compatibility
  scores; genuine ProteinMPNN and RNA-context LigandMPNN sampling covers all
  four Stage-0003A representatives.
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
- stage-0002 audit: `docs/STAGE_0002_SUMMARY_2026-07-31.md`;
- stage-0003A audit: `docs/STAGE_0003A_SUMMARY.md`;
- multi-scaffold report: `reports/stage_0003a/report.md` and `report.html`;
- variant data card: `docs/CAS13D_VARIANT_DATA_CARD.md`;
- frozen candidate protocol:
  `experiments/preregistered/pre_wetlab_candidate_protocol.yaml`;
- Stage-0003 GPU input root: `artifacts/gpu_jobs/stage_0003/`.

## Current quality gate

- Ruff lint: pass;
- Ruff format: pass;
- strict mypy: pass for 52 source files in both the isolated conda environment
  and a fresh CI-equivalent Python 3.11 pip environment;
- pytest: 79/79 pass;
- branch-aware coverage: 70.93% (required 70%);
- matched-report acceptance: pass for 18 candidates, nine methods, two common
  seed blocks, one fixed/free mask, no mocks, and zero fixed violations;
- Stage-0003A local commands `make lint`, `make typecheck`, `make test`,
  `make smoke-cpu`, `make smoke-esm-if1`, `make smoke-proteinmpnn`,
  `make smoke-ligandmpnn`, and `make verify-reproducibility`: pass at execution
  commit `437ca988ff92bc4ab0e728df5c706fbf817f0754`;
- clean Atlas parse and cluster leakage gate: pass;
- real candidate novelty audit: pass with fail-closed missing-hit semantics;
- GitHub Actions run `30633893318`: primary fixture/code validation success;
- GitHub Actions run `30635299503`: preceding-stage final handoff/documentation
  validation success; it is not the primary code-run identifier;
- GitHub Actions run `30642823247`: stage-0002 initial push failed strict mypy
  because the fresh pip environment lacked pandas stubs and required explicit
  NumPy dtype annotations; the failure remains visible;
- GitHub Actions run `30643470631`: stage-0002 repaired code/results validation
  success at exact commit `b1a33fffd1d7a9254a55247c449017f240f02683`;
- GitHub Actions run `30662641614`: Stage-0003A first push failed strict mypy
  because current pip `pandas-stubs` widened `itertuples()` field types; the
  failure remains visible and no scientific output was affected;
- GitHub Actions run `30662821387`: Stage-0003A repaired code/results validation
  success at exact commit `b251adf` after an explicit type-only narrowing;
- tag `v0.1.0-data-pipeline` points to final data-pipeline HEAD
  `a9a530d14434e74dc0cfc47896847e201431c1c2`;
- final-HEAD GPU bundle: `gpu-bundle-a9a530d14434-7540febfb2`, export
  `dirty=false`, `missing_assets=[]`, internal hashes and 14/14 source-asset
  hashes passed;
- Stage-0003A execution bundle:
  `gpu-bundle-437ca988ff92-0a8201ce8f`, export `dirty=false`,
  `missing_assets=[]`, all 1,204 embedded input/internal hashes and 28/28
  source-asset hashes passed;
- target GPU validation: not run.

The final documentation-HEAD bundle and its post-commit CI are replayed without
tracked edits; their exact immutable identifiers are reported in the operator
handoff. No large refold task has been started.
