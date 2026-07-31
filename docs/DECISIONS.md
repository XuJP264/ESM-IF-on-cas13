# Decision Log

## 2026-07-31 — Evidence terminology

Adopt Levels 0–4 globally. Levels 1–3 remain computational candidates; Level 4
requires wet-lab evidence outside the present scope.

## 2026-07-31 — Experimental structures first

Benchmark 6E9F and 5XWP before large-scale structure prediction. This supplies
real backbones and RNA-context annotations without waiting for AlphaFold.

## 2026-07-31 — Local asset loading

All formal model inference uses explicitly fetched local checkpoints. Runtime
network downloads are errors because they prevent reliable hashing and audit.

## 2026-07-31 — Split isolation

Use sequence-cluster-level splits with a strict 40% identity generalization
test. Auxiliary 50%/70% splits may be added if subtype or structure coverage is
insufficient, but cannot replace the strict test.

## 2026-07-31 — Measured local resource limits

Plan for a measured RTX 4060 Laptop GPU with 8188 MiB VRAM, not 32 GB. Use
single-structure/small batches, CPU fallback where algorithmically equivalent,
and GPU/HPC exports for jobs that do not fit.

## 2026-07-31 — Project-local Conda state

Use `.tools/pkgs`, `.tools/cache`, and `.tools/envs/*` so environment creation
does not depend on writable global Miniforge or user cache directories. Install
the project from its absolute repository path after Conda resolves the
environment; relative `-e .` entries inside `envs/*.yml` are ambiguous because
Conda resolves them relative to the environment-file directory.

## 2026-07-31 — M0 quality scope

Apply strict mypy to the `cas13_if` package and require at least 70% branch-aware
coverage for CPU fixture tests. M0 passed at 79.67%. Genuine model/network tests
remain explicitly marked and excluded from CPU CI.

## 2026-07-31 — Catalytic construct mutations

Treat deposited experimental construct sequences separately from literature-
supported biological sequences. 6E9F carries alanine at all four annotated
HEPN catalytic positions, and 5XWP carries alanine at its second R/H pair.
Design masks restore/protect the literature-supported R/H tokens; reports show
construct recovery and biological recovery separately.

## 2026-07-31 — Atlas repeat orientation

The inspected Atlas v1.0 schema does not provide repeat orientation. Unknown
orientation remains ambiguous and is excluded from high-confidence paired
protein/direct-repeat analyses. Do not silently select a strand or use spacer
content as proof of orientation.

## 2026-07-31 — Atlas subtype resolution preserves raw disagreements

Atlas v1.0 frequently stores only generic `VI` (or an empty subtype) while the
Cas HMM contains an explicit family such as `CAS-VI-D`. Resolve Cas13 subtype
from an explicit Cas HMM first and, where the HMM lacks a subtype, from
canonical Cas13a–j nomenclature. Preserve `subtype_raw`, `subtype_source`, and
`subtype_conflict` on every Cas13 record. A clear non-Type-VI summary that
conflicts with a Cas13 HMM is retained but routed to ambiguous pairing; it is
never silently promoted to a high-confidence pair.

## 2026-07-31 — Evolutionary inputs require auditable complete records

Retain all Cas13 HMM hits in `cas13_records` and exact-unique tables, including
short, truncated, or subtype-conflicting annotations. Preserve HMM e-value,
score, source length field, and two-end truncation flag. For subtype MSA
eligibility, require at least one occurrence with no subtype conflict and the
explicit Atlas completeness flag `truncated=00`. Select one eligible member per
70% cluster and subtype even if MMseqs chose an ineligible member as the raw
cluster representative. This is a data-QC rule, not evidence of Cas13 function.

## 2026-07-31 — Environment isolation gate

An environment is not accepted merely because its Conda transaction completes.
Imports must pass with `PYTHONNOUSERSITE=1`, and dependency paths must resolve
inside the declared prefix. The first ESM environment failed this gate because
pip reused `fair-esm` from the user site. The repair installs ESM from the
pinned local upstream checkout into the environment prefix and regenerates the
lock.

## 2026-07-31 — Missing-coordinate semantics across MPNN backends

Do not require ProteinMPNN's upstream residue-number tensor to have the same
length as the strict observed-coordinate sequence. For 6E9F chain A, the
strict parser sees 864 resolved residues while ProteinMPNN builds 893 residue
number slots and masks 29 missing-coordinate positions. Validation must prove
that removing the upstream `X` slots recovers the strict sequence and must
report both lengths; it must not silently delete or impute missing coordinates.

## 2026-07-31 — LigandMPNN upstream-compatible environment

Pin LigandMPNN to NumPy 1.23.5 and include `dm-tree`, matching the pinned
upstream requirements. The first isolated smoke exposed missing `tree`; after
adding it, the pinned upstream OpenFold code exposed its use of removed
`np.int` under NumPy 1.26. Resolve these by environment pinning, not by editing
third-party source.

## 2026-07-31 — Candidate identifiers encode sampling conditions

Candidate IDs must distinguish backend, scaffold, parent sequence, temperature,
seed, fixed positions, allowed-residue filters, conditioning chains, and sample
index. Do not use only scaffold plus sample index: it collides across
temperatures and ablations and makes refold manifests ambiguous. Paths are
excluded from the digest so the same declared request has the same identifier
after migration to another node.

## 2026-07-31 — Pilot recovery is not a matched method comparison

The first real ESM benchmark varies the number of hard-fixed positions across
conditions. Report its recovery values as pipeline/pilot characterization only.
Do not interpret higher raw recovery in a more heavily fixed condition as a
method improvement. Formal H1–H4 comparisons must match scaffold, design
positions, parent identity/novelty range, and statistical unit.

## 2026-07-31 — Bioinformatics and Git LFS remain project-local

Use the locked `.tools/envs/bioinformatics` environment for MMseqs2, MAFFT,
HMMER, Infernal, seqkit, Foldseek, TM-align, and Git LFS. Git LFS is available
for repository tooling but model weights, raw Atlas data, PDFs, and prediction
outputs remain ignored rather than added to LFS.

## 2026-07-31 — Atlas completeness flags do not establish full-length Cas13

The first inclusive subtype MSA retained records marked `truncated=00`, but its
shortest representatives were 48–80 aa and no subtype had an alignment column
with 90% sequence coverage. Treat that alignment as an audit result, not a
design constraint source. Before any scaffold mapping or candidate test metric,
preregister a broad 700–1600 aa full-length screen, retain canonical sequences
with a nonconflicting complete occurrence, and select the longest eligible
member of each 70% MMseqs2 cluster with SHA256 as the deterministic tie-break.
Column conservation still requires an explicit scaffold mapping, mapping
confidence, and at least 80% column coverage; conservation alone never
automatically creates a hard-fixed position. Preserve the inclusive alignment
outside the canonical result path and report length-threshold sensitivity.
