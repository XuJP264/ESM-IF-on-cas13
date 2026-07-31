# CRISPR-Cas Atlas / experimental-structure data card

Last updated: 2026-07-31

## Scope and licenses

The primary source is official CRISPR-Cas Atlas v1.0, licensed CC BY-NC 4.0.
The project MIT license applies only to project-authored code and does not
relicense Atlas records, RCSB structures, papers, or model weights. Raw JSON,
large processed tables and structure files remain local and are not committed.

Official Atlas asset:

- local path: `data/raw/atlas/v1.0/crispr-cas-atlas-v1.0.json`;
- size: 5,267,508,328 bytes;
- SHA256:
  `5b4ba2fb99638d279e0c126100e19a4b77aba487b37b7df118e4bf4acd494720`;
- manifest: `data/manifests/atlas_v1.0.yaml`;
- manifest status: `downloaded_verified`.

Experimental structures 6E9F, 5XWP, 6E9E and 5XWY were obtained through RCSB
and retain their upstream terms and citations.

## Processing and record semantics

The top-level Atlas JSON array is decoded incrementally. Operon, CRISPR array
and effector tables are written in batches; exact Cas13 dedup uses a local
SQLite aggregation. The pipeline does not assume that:

- record count equals independent sequence count;
- every Type VI operon contains one Cas13;
- every operon contains one CRISPR array;
- repeat orientation is known;
- a spacer is a verified natural target RNA.

The clean production run is:

`results/runs/20260731-atlas-processing-e8356ef7b5-eebc1a5-r001/`

It records `dirty=false`, input/output hashes, exit code 0 and SUCCESS.

## Production funnel

| Item | Count |
|---|---:|
| Atlas operons | 1,246,088 |
| Type VI operons | 11,707 |
| Cas effector annotations | 6,174,375 |
| Cas13 records | 12,353 |
| Cas13 exact-unique sequences | 4,070 |
| Evolution-eligible exact-unique sequences | 3,500 |
| High-confidence Cas13–direct-repeat pairs | 0 |
| Ambiguous pairs | 11,727 |
| Processing failures | 0 |

Resolved Cas13 record subtypes:

| Subtype | Records |
|---|---:|
| VI-B | 5,163 |
| VI-D | 6,857 |
| VI-F | 166 |
| VI-I | 167 |

Forty records have a conflict between the operon summary and explicit effector
annotation. Raw subtype, resolution source and conflict state are retained.

## Direct-repeat limitation

The inspected Atlas v1.0 records do not provide a recoverable repeat orientation.
The high-confidence policy requires exactly one unambiguous Cas13, one array, a
nonempty repeat, an explicit subtype and declared/reliably recovered orientation.
Consequently, all otherwise relevant records remain in
`ambiguous_pairs.parquet`; none are silently strand-flipped or used for
coevolution. Real paired MSA, MI/APC and DCA are data-blocked until a trustworthy
orientation source is available.

## Deduplication, clustering and split

Exact uniqueness uses the full amino-acid sequence SHA256. MMseqs2 clustering
uses coverage 0.8, coverage mode 0, cluster mode 2 and 16 threads.

| Minimum identity | Clusters |
|---:|---:|
| 100% | 3,877 |
| 90% | 1,797 |
| 70% | 1,323 |
| 50% | 1,003 |
| 40% | 783 |
| 30% | 516 |

The strict 40% cluster split contains 3,335 train, 160 validation and 575 test
sequences; the cluster leakage gate passed. A 100% MMseqs cluster is not the
same as exact dedup because the registered 80% coverage permits fragment/full
length relationships.

Subtype-held-out and scaffold-held-out auxiliary splits remain to be
materialized. Test clusters may not be used for threshold, mask, temperature
or model selection.

## Evolutionary-input quality

Atlas `truncated=00` does not guarantee a full-length Cas13: the inclusive MSA
contained representatives as short as 48–80 aa and no subtype had a 90%
coverage column. That inclusive output is retained for audit only.

Before candidate test metrics, the project preregistered a broad 700–1600 aa
screen, required canonical amino acids and an unconflicted complete occurrence,
and selected the longest eligible member of each 70% cluster/subtype.

The resulting real MAFFT sets are:

- VI-B: 489 sequences;
- VI-D: 182;
- VI-F: 45;
- VI-I: 50.

Conservation tables retain weighted/unweighted frequencies, entropy, gap
fraction, coverage and effective sequence count. A column requires coverage
of at least 0.8 plus a successful scaffold mapping and mapping-confidence gate
before it can inform constraints. Conservation never automatically proves a
catalytic or functional position.

## Known biases and nonclaims

- Atlas discovery/HMM/taxonomic sampling biases are not corrected by exact
  deduplication alone.
- The observed Atlas Cas13 table has no VI-A records, limiting interpretation
  of 5XWP/Cas13a candidate similarity.
- Low-count VI-F and VI-I alignments have less statistical support than VI-B/D.
- MSA length/gap sensitivity and scaffold mapping remain pending.
- Neither Atlas annotation, sequence novelty nor inverse-folding score proves
  Cas13 activity.

