# Stage 0003A summary

Status: Stage implementation and local acceptance are complete as of
2026-08-01 (Asia/Shanghai). Execution commit `b251adf` passed the clean pip-CI
run `30662821387`; the exact final documentation-HEAD bundle and CI identifiers
are published in the operator handoff. The detailed generated report is
`reports/stage_0003a/report.md` and `report.html`.

## Evidence boundary

The highest real evidence in this stage is Level 2. No real AlphaFold2,
ColabFold, AlphaFold3, or Boltz result has been produced. The prediction fixture
suite is `is_mock=true` and supports Level 0 only. No candidate is called
functional, active, effective, or wet-lab validated.

## Real completed work

- Four independent VI-D parents and nine experimental state units were
  downloaded from official RCSB endpoints: EsCas13d (6E9E, 6E9F), UrCas13d
  (6IV9), DjCas13d (9M38, 9M30, 9M33, 9M34), and CasRx/RfxCas13d (9M31,
  9M8Q). This covers one apo, four binary, and four ternary variants/states.
- All nine states have strict PDB-residue ↔ coordinate-sequence ↔ full-natural-
  sequence ↔ VI-D-MSA mappings and RNA contact/second-shell annotations.
  Weighted high-confidence four-layer coverage is 7,507/8,272 (90.75%).
- Nine native state sequences were genuinely scored by ESM-IF1 on the local
  RTX 4060. Thirteen valid state combinations were aggregated without averaging
  coordinates; fixed-token violations are zero.
- The Cas13d variant dataset contains 22 records, 20 mutation strings, four
  studies, and four non-pooled assay groups. All full mutant sequences are
  recoverable; 15 cis/trans values are approximate readings from a public
  figure and are labeled accordingly.
- The retrospective benchmark genuinely scored 10 point variants (plus WT
  reference) with 22 ESM state scores, 11 ProteinMPNN scores, and 11
  RNA-context LigandMPNN scores. Six indels were explicitly excluded from an
  invalid unchanged-backbone comparison. Only nine point variants have numeric
  activity, so results are descriptive; significance tests were not run.
- The local multi-scaffold smoke generated 48 real rows: four parents × six
  methods × two seeds. It has mock=0, runtime failures=0, fixed violations=0,
  and eight of eight LigandMPNN rows with verified RNA atomic context.
- Candidate-selection thresholds were preregistered before real Level-3
  outputs. The production input inventory contains 18 Stage-0002 pilots, 48 new
  real smoke sequences, and four WT controls.

## Fixture/mock completed work

Four `is_mock=true` AlphaFold3/Boltz-style results passed strict ingest and QC.
The E2E invoked the genuine local TM-align binary and exercised pLDDT, PAE,
domain RMSD, HEPN geometry, RNA-contact recovery, interface confidence,
multi-seed and cross-model consistency, retries, and Pareto ranking. Perfect
self-comparison values are fixture checks and have no scientific meaning.

## H100 preparation

The 70 proteins expand to 1,068 deterministic, two-seed jobs: 560 monomer, 280
binary, and 228 ternary jobs for ColabFold/AlphaFold2, AlphaFold3, and Boltz.
Every job forbids forcing its target scaffold as a template. The single tmux
entry is:

```bash
bash scripts/launch_gpu_tmux.sh \
  configs/stage_0003_refold.yaml \
  stage-0003-refold
```

Predictor installation, licensed databases/weights, and the site adapter are
node-specific prerequisites. The dispatcher refuses GPUs below 40,000 MiB,
validates job hashes, and refuses to start without an explicit executable site
adapter. This local 8 GB node was not used for large prediction.

## Audited failures

- The first multi-state attempt failed because user-site PyTorch shadowed the
  project environment. `PYTHONNOUSERSITE=1` plus the pinned ESM environment
  repaired it; both runs remain in `results/runs/`.
- The first retrospective attempt stopped at the LigandMPNN RNA-context gate
  because context and sampling metadata occupy different FASTA headers. The
  validator now requires the upstream context attestation across all headers;
  the repaired real run passed.
- UrCas13d natural-sequence confidence is medium because its RCSB-linked
  UniProt accession currently has no sequence response. DjCas13d and CasRx
  natural-parent proxies are also medium because the primary paper declares
  inactive catalytic constructs that RCSB entity metadata does not.

## Readiness

The platform is ready to transfer Level-3 jobs, but it has **not reached the
pre-wet-lab candidate standard**. Real Level-3 monomer and RNA-complex outputs
must pass the frozen gates before any shortlist can advance. Level 4 remains
unavailable and outside this stage.
