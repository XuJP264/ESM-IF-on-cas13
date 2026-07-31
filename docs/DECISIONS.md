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

## 2026-07-31 — Atlas novelty searches fail closed on missing coverage hits

For candidate novelty, search the full 4,070-sequence exact-unique Atlas FASTA
with MMseqs2 and require at least 80% query coverage. Report the maximum
identity only among returned alignments and preserve candidates with no such
hit as `no_atlas_hit_at_required_query_coverage`; do not reinterpret missing
hits as proof of extreme novelty. A Level 1 row must also pass the registered
parent-identity, Atlas-identity, homopolymer, low-complexity, entropy, and
fixed-position gates. The absence of VI-A in the current Atlas Cas13 HMM table
is a database-coverage limitation for 5XWP, not favorable evidence for its
candidates.

## 2026-07-31 — Final data-pipeline checkpoint and CI roles

Tag `v0.1.0-data-pipeline` identifies exact commit
`a9a530d14434e74dc0cfc47896847e201431c1c2`. Its authoritative GPU source
bundle is `gpu-bundle-a9a530d14434-7540febfb2`; the export was clean, had no
missing assets, and passed internal hashes plus all 14 source-asset hashes.
GitHub Actions run `30633893318` remains the primary fixture/code CI evidence;
run `30635299503` is the final handoff/documentation CI evidence. Do not cite
the latter as if it were the primary code-validation run.

## 2026-07-31 — Fair baseline masks are method-invariant

For the 6E9F VI-D matrix, all methods use one identical hard-fixed set and one
identical free set. Conservation or RNA-interface information may bias allowed
residues or proposal probabilities only inside the common free set; it may not
change the evaluated design positions. `unconstrained_esm_if1` means no such
proposal bias beyond the common safety mask. The requested
`catalytic_only_fixed_esm_if1` condition is a same-seed deterministic control
of that same mask and is not interpreted as a separate biological treatment;
its equality with `unconstrained_esm_if1` is expected and must not be counted as
independent evidence. Conservation remains disabled until a high-confidence full
scaffold-to-VI-D-MSA mapping passes its audit.

## 2026-07-31 — Identity matching uses a preregistered common interval

Before the formal matched run, calibrate only on historical/diagnostic identity
values, never endpoint performance. Select one proposal per method and seed in
the common parent- and designed-position identity interval 0.18–0.32, nearest
to target 0.25 with candidate ID as deterministic tie-break. Use random
mutation probability 0.78, profile temperature 2.0, and model temperature 1.0.
Limit the local conservation and RNA-contact proposal filters to the top 48
mapped positions each so a real small CPU matrix can meet the same novelty
interval without post-hoc sequence editing. The larger GPU extension may add
seeds but may not change this matching rule after seeing performance metrics.

## 2026-07-31 — Consensus matching may choose only genuine source tokens

The ESM-IF1/LigandMPNN consensus preserves exact source agreement. Where the
sources disagree and exactly one retains the parent token, it ranks those
parent choices by the source-confidence advantage until the preregistered
parent-identity target is reached as closely as feasible. Other disagreements
use the higher-confidence source token. It never invents a third amino acid.
This identity-only amendment was registered after a fail-closed proposal audit
and before conditional likelihood, novelty, recovery, or other endpoint
calculation. Oversampling seeds are `seed_block + 1,000,000 * proposal_index`
to prevent adjacent blocks from sharing an actual model seed.

## 2026-07-31 — GPU asset and execution commits are distinct

The verified 14-asset bundle intentionally pins the requested immutable
data-pipeline commit `a9a530d14434e74dc0cfc47896847e201431c1c2`. The matched
GPU extension requires the later execution-code commit recorded in
`reports/matched_baselines/gpu_hpc_job_manifest.jsonl`. Verify/synchronize the
asset bundle at its pinned commit first, then checkout the recorded execution
commit before bootstrapping and launching. Do not treat the asset-bundle commit
as proof that later experiment code was included in that bundle.

## 2026-07-31 — CI type dependencies must be explicit in pip extras

Passing mypy in the Conda analysis environment is insufficient if that
environment contains type packages absent from `.[dev]`. Keep `pandas-stubs`
in the pip development extra and use `numpy.typing.NDArray` with explicit
dtypes. Validate dependency-sensitive type changes once in a fresh Python 3.11
pip-only environment so GitHub Actions does not depend on hidden Conda state.

## 2026-08-01 — Independent scaffold and state semantics

Treat EsCas13d, UrCas13d, DjCas13d, and CasRx/RfxCas13d as four independent
natural parents. Different apo, binary, ternary, mismatch, or trans-RNA states
of one parent are state units, not additional scaffold replicates. Never average
their atomic coordinates into a synthetic backbone; score each state and
combine normalized log-likelihoods only after state-wise evaluation.

## 2026-08-01 — Natural-parent restoration is evidence-qualified

Separate deposited inactive constructs from natural-parent proxies. Restore
paper-supported catalytic residues for EsCas13d, UrCas13d, DjCas13d, and CasRx,
but assign medium sequence-provenance confidence when RCSB lacks the mutation
declaration or its UniProt cross-reference no longer supplies a sequence. Show
every substitution and source reason in mapping/manual-review output. Do not
treat a restored proxy as an independently sequenced construct.

## 2026-08-01 — Multi-state safety masks have distinct meanings

The intersection hard mask contains literature-supported HEPN residues and
high-confidence direct RNA contacts present in every observed state. The union
risk mask records any mapping failure, core, direct/second-shell RNA contact,
or chain-break-adjacent position in any state. The variable hinge mask is the
union-minus-intersection difference. Core annotations do not automatically
enter the hard mask; otherwise an experimental single-state core call would
freeze most of a scaffold and obscure state variability.

## 2026-08-01 — Variant assays remain non-poolable

Keep purified cis/trans cleavage, reporter knockdown, collateral activity,
cellular allele-specific suppression, and mouse phenotypes in separate
comparability groups. Figure-read numerical values remain approximate. Indels
are excluded from unchanged-backbone inverse-folding scoring. With only nine
numeric point variants, correlations are descriptive and no significance or
predictive-validity claim is allowed.

## 2026-08-01 — Low-complexity QC is parent-aware

Fail a designed sequence for low complexity when it introduces a flagged
window beyond the natural-parent baseline. Preserve and report the absolute
window count, but do not reject all variants of a parent solely because an
unresolved or inherited natural segment already triggers the same heuristic.
Homopolymer and composition gates remain absolute as preregistered.

## 2026-08-01 — Level-3 preparation is not Level-3 evidence

Stage-0003 manifests, shards, backend-formatted inputs, expected outputs, and
retry manifests are Level 0 preparation. Mock prediction fixtures remain
`is_mock=true`, even when the metric implementation invokes genuine TM-align.
Only real independent monomer and RNA-complex outputs passing the frozen gates
may support Level 3. Conditional ESM-IF1 likelihood alone may never determine a
pre-wet-lab ranking.

## 2026-08-01 — Large refolds are H100-only and site-explicit

The Stage-0003 dispatcher requires at least 40,000 MiB visible GPU memory,
verifies input hashes, and requires an explicit executable site adapter for
the installed predictor/database layout. The local 8 GB RTX 4060 must fail
before prediction starts. Every job uses
`no_target_scaffold_as_forced_template`; ordinary backend template search may
be used only without forcing the target experimental scaffold.
