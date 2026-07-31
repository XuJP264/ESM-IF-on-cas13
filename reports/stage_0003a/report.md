# Stage 0003A report

Status: Stage implementation and local acceptance complete. Immutable final
bundle and CI identifiers are published in the operator handoff. Date:
2026-08-01 (Asia/Shanghai).

## Claim boundary

The real calculations in this report support no more than Level 2
inverse-folding compatibility. The Level-3 prediction count is **0**. Prediction
fixtures are marked `is_mock=true` and support only Level 0 code-path evidence.
No sequence is described as functional, active, effective, or wet-lab
validated.

## Experimental multi-scaffold structure atlas — real data

The atlas contains 4 independent natural Cas13d
parents and 9 experimental scaffold-state units. All coordinates
were downloaded from RCSB, hashed, and audited; all RNA-bearing states have
non-zero retained RNA atom counts.

| scaffold | PDB | state | resolution_A | protein | crRNA | target_RNA | coordinate_length | RNA_atoms |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| EsCas13d | 6E9E | binary | 3.4 | A | B |  | 863 | 1088 |
| EsCas13d | 6E9F | ternary | 3.3 | A | B | C | 864 | 1679 |
| UrCas13d | 6IV9 | binary | 1.86 | A | B |  | 874 | 1054 |
| DjCas13d | 9M38 | apo | 3.47 | A |  |  | 716 | 0 |
| DjCas13d | 9M30 | binary | 3.05 | A | B |  | 844 | 1085 |
| DjCas13d | 9M33 | ternary_tr1 | 3.27 | A | B | C | 750 | 1444 |
| DjCas13d | 9M34 | ternary_tr2 | 3.46 | A | B | C | 750 | 1444 |
| CasRx | 9M31 | binary | 2.86 | A | B |  | 945 | 1082 |
| CasRx | 9M8Q | ternary | 3.46 | A | B | C | 904 | 1696 |

## Four-layer mapping — real data and calculations

Across state-specific full sequences, 7,507
of 8,272 positions pass
the strict four-layer exact/restored gate, a weighted coverage of
90.75%. Conservation
is fail-closed elsewhere. Manual-review CSV and HTML files are in
`reports/stage_0003a/manual_review/` and each state mapping is in
`reports/stage_0003a/residue_mapping/`.

| pdb_id | scaffold_id | state | full_scaffold_length | mapped_coordinate_positions | four_layer_exact_or_restored_positions | high_confidence_coverage | unresolved_positions | RNA_contact_positions | RNA_second_shell_positions |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 6E9E | EsCas13d | binary | 954 | 863 | 863 | 0.9046 | 91 | 137 | 132 |
| 6E9F | EsCas13d | ternary | 954 | 864 | 864 | 0.9057 | 90 | 128 | 135 |
| 6IV9 | UrCas13d | binary | 922 | 874 | 871 | 0.9447 | 48 | 127 | 145 |
| 9M38 | DjCas13d | apo | 877 | 716 | 716 | 0.8164 | 161 | 0 | 0 |
| 9M30 | DjCas13d | binary | 877 | 844 | 844 | 0.9624 | 33 | 134 | 130 |
| 9M33 | DjCas13d | ternary_tr1 | 877 | 750 | 750 | 0.8552 | 127 | 67 | 85 |
| 9M34 | DjCas13d | ternary_tr2 | 877 | 750 | 750 | 0.8552 | 127 | 59 | 85 |
| 9M31 | CasRx | binary | 967 | 945 | 945 | 0.9772 | 22 | 152 | 164 |
| 9M8Q | CasRx | ternary | 967 | 904 | 904 | 0.9349 | 63 | 121 | 144 |

## Multi-state native scoring — real ESM-IF1

Nine state-specific native sequences were genuinely scored on the local RTX
4060 with the isolated ESM-IF1 checkpoint. The implementation combines scores
after state-wise evaluation; it never averages coordinates. Thirteen applicable
state combinations were aggregated and 11 unavailable combinations were
recorded as `not_applicable`; fixed-token violations were zero.

| scaffold_id | pdb_id | state | mean_log_likelihood_per_resolved_residue | perplexity |
| --- | --- | --- | --- | --- |
| EsCas13d | 6E9E | binary | -2.2303 | 9.3024 |
| EsCas13d | 6E9F | ternary | -2.1019 | 8.1816 |
| UrCas13d | 6IV9 | binary | -2.3860 | 10.8703 |
| DjCas13d | 9M38 | apo | -2.2447 | 9.4379 |
| DjCas13d | 9M30 | binary | -2.3806 | 10.8116 |
| DjCas13d | 9M33 | ternary_tr1 | -2.2397 | 9.3902 |
| DjCas13d | 9M34 | ternary_tr2 | -2.2448 | 9.4382 |
| CasRx | 9M31 | binary | -2.5221 | 12.4543 |
| CasRx | 9M8Q | ternary | -2.3898 | 10.9116 |

## Cas13d variant/activity resource — real literature evidence

The version-1 resource contains 22 records from
4 studies, 20 unique
mutation strings, and 4 non-poolable assay
groups. All 22 mutant sequences were recovered. Numerical
2025 Figure-6 values are approximate graph readings and are not raw replicates.

The same-backbone retrospective benchmark genuinely scored
10 point variants; 6
indels were excluded rather than mis-scored. It produced
22 ESM state scores,
11 ProteinMPNN scores, and
11 RNA-context LigandMPNN scores. The two-state
rank consistency was 0.5152. With only
n=9 numerical variants, all correlations below are descriptive and no
significance tests were run.

| endpoint | metric | spearman_rho | n | p_value | inference |
| --- | --- | --- | --- | --- | --- |
| cis_activity | esm_binary_delta_vs_wt | -0.0667 | 9 |  | descriptive_only |
| cis_activity | esm_ternary_delta_vs_wt | -0.4167 | 9 |  | descriptive_only |
| cis_activity | multi_state_min_score | -0.0667 | 9 |  | descriptive_only |
| cis_activity | proteinmpnn_delta_nll_vs_wt | 0.1500 | 9 |  | descriptive_only |
| cis_activity | ligandmpnn_delta_log_probability_vs_wt | -0.3167 | 9 |  | descriptive_only |
| cis_activity | mean_mutated_position_conservation | 0.5667 | 9 |  | descriptive_only |
| trans_activity | esm_binary_delta_vs_wt | -0.4770 | 9 |  | descriptive_only |
| trans_activity | esm_ternary_delta_vs_wt | -0.2176 | 9 |  | descriptive_only |
| trans_activity | multi_state_min_score | -0.4770 | 9 |  | descriptive_only |
| trans_activity | proteinmpnn_delta_nll_vs_wt | 0.1255 | 9 |  | descriptive_only |
| trans_activity | ligandmpnn_delta_log_probability_vs_wt | -0.4519 | 9 |  | descriptive_only |
| trans_activity | mean_mutated_position_conservation | 0.7448 | 9 |  | descriptive_only |

The score/activity directions are mixed and do not establish predictive
validity. In particular, conservation showed the largest positive descriptive
association in this small panel (cis rho 0.5667; trans rho 0.7448), while model
scores were not consistently monotonic with activity. This is a negative/limited
result, not evidence for selecting a functional protein.

## Local multi-scaffold generation — real bounded smoke

Exactly 48 real rows were generated: four scaffolds, six methods, two
seeds, and one proposal per seed. Every method/scaffold pair has two rows;
mock rows, runtime failures, and hard-fixed violations are all zero. All eight
LigandMPNN rows explicitly retained RNA atomic context. Each coordinate sequence
was also genuinely scored by ESM-IF1 on its representative state.

| scaffold_id | method | n | parent_identity_min | parent_identity_max | mean_esm_ll | fixed_violations | rna_context_rows |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CasRx | common_safety_mask_esm_if1 | 2 | 0.3113 | 0.3289 | -2.1913 | 0 | 0 |
| CasRx | conservation_esm_if1 | 2 | 0.3382 | 0.3464 | -2.1964 | 0 | 0 |
| CasRx | conservation_rna_esm_if1 | 2 | 0.3547 | 0.3733 | -2.1887 | 0 | 0 |
| CasRx | esm_ligand_consensus | 2 | 0.3919 | 0.3919 | -2.1880 | 0 | 2 |
| CasRx | ligandmpnn | 2 | 0.3444 | 0.3506 | -2.7585 | 0 | 2 |
| CasRx | proteinmpnn | 2 | 0.3195 | 0.3237 | -2.8491 | 0 | 0 |
| DjCas13d | common_safety_mask_esm_if1 | 2 | 0.3090 | 0.3101 | -2.0659 | 0 | 0 |
| DjCas13d | conservation_esm_if1 | 2 | 0.3330 | 0.3387 | -2.0744 | 0 | 0 |
| DjCas13d | conservation_rna_esm_if1 | 2 | 0.3535 | 0.3660 | -2.1242 | 0 | 0 |
| DjCas13d | esm_ligand_consensus | 2 | 0.4436 | 0.4436 | -2.0216 | 0 | 2 |
| DjCas13d | ligandmpnn | 2 | 0.3592 | 0.3660 | -2.7757 | 0 | 2 |
| DjCas13d | proteinmpnn | 2 | 0.3318 | 0.3330 | -2.7941 | 0 | 0 |
| EsCas13d | common_safety_mask_esm_if1 | 2 | 0.3690 | 0.3784 | -2.0183 | 0 | 0 |
| EsCas13d | conservation_esm_if1 | 2 | 0.3899 | 0.4004 | -2.0014 | 0 | 0 |
| EsCas13d | conservation_rna_esm_if1 | 2 | 0.4088 | 0.4193 | -2.0152 | 0 | 0 |
| EsCas13d | esm_ligand_consensus | 2 | 0.4109 | 0.4109 | -2.0162 | 0 | 2 |
| EsCas13d | ligandmpnn | 2 | 0.3920 | 0.4067 | -2.5199 | 0 | 2 |
| EsCas13d | proteinmpnn | 2 | 0.3784 | 0.3795 | -2.5909 | 0 | 0 |
| UrCas13d | common_safety_mask_esm_if1 | 2 | 0.4208 | 0.4219 | -1.6984 | 0 | 0 |
| UrCas13d | conservation_esm_if1 | 2 | 0.4436 | 0.4458 | -1.7696 | 0 | 0 |
| UrCas13d | conservation_rna_esm_if1 | 2 | 0.4620 | 0.4620 | -1.7887 | 0 | 0 |
| UrCas13d | esm_ligand_consensus | 2 | 0.3839 | 0.3839 | -2.1085 | 0 | 2 |
| UrCas13d | ligandmpnn | 2 | 0.4892 | 0.4935 | -2.9462 | 0 | 2 |
| UrCas13d | proteinmpnn | 2 | 0.4501 | 0.4664 | -3.0446 | 0 | 0 |

The historical absolute low-complexity flag is set for
12/48 rows, all 12 UrCas13d rows,
because the natural Ur parent itself contains the flagged windows. Candidate
selection therefore uses the preregistered, parent-aware rule of zero *new*
low-complexity windows. These smoke rows are coverage checks, not a final
statistical comparison or a wet-lab shortlist.

## Level-3 job preparation — real inputs, prediction not run

The inventory has 70 proteins: 18
Stage-0002 pilot candidates, 48
new real smoke candidates, and 4 WT controls. It expands
to 1068 deterministic jobs over two seeds: 560 monomer, 280
binary, and 228 ternary jobs across ColabFold, AlphaFold2, AlphaFold3, and Boltz.
All use `no_target_scaffold_as_forced_template`. The target scaffold is never supplied as a
forced template. Status is `prepared_not_run`; manifests are not predictions.

## Level-3 ingest/ranking — MOCK fixture only

The fixture E2E ingested 4 labeled mock results and
used the genuine local TM-align executable for 4
fixture comparisons. It exercised pLDDT, PAE, domain RMSD, HEPN geometry,
RNA-contact recovery, interface confidence, multi-seed consistency,
cross-model consistency, missing-output/retry behavior, and six-dimensional
Pareto ranking. `real_prediction_count=0` and every
fixture artifact remains `is_mock=true`; its perfect self-comparison metrics
have no scientific meaning.

## Audited failures and limitations

- One initial multi-state run failed because user-site PyTorch 2.12/CUDA 13.0
  shadowed the isolated environment. The failed run was retained; rerunning with
  `PYTHONNOUSERSITE=1` used pinned PyTorch 2.4.1/cu121 and succeeded on CUDA.
- One initial retrospective run stopped when validation inspected only the last
  LigandMPNN FASTA header. Upstream records context on the input header and
  sampling statistics on generated headers. The gate was corrected to require
  the explicit `use_ligand_context=True` attestation on any upstream header;
  the failed run remains retained.
- UrCas13d natural sequence provenance is medium confidence because its RCSB
  UniProt accession currently returns no sequence; the 922-residue RCSB
  reference was restored from primary-paper catalytic annotations and the
  expression tag was excluded.
- DjCas13d and CasRx RCSB metadata omit mutations described in the primary
  paper. Their four catalytic residues were restored with medium-confidence
  natural-parent proxies, visibly recorded in the state manifest.
- Domain boundaries and domain-interface labels remain manual-review fields,
  not inferred negatives. No real AlphaFold/AF3/Boltz prediction was run.

## Readiness decision

The project is ready to transfer deterministic Level-3 jobs to an H100 node,
but **has not reached the pre-wet-lab candidate standard**. It remains at Level
2 because no real independent monomer or RNA-complex prediction has passed the
predeclared Level-3 gates. It is not ready to nominate an experimentally
supported Cas13 and makes no efficacy claim.
