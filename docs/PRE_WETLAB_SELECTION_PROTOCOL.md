# Pre-wet-lab candidate selection protocol

Version 1 of this protocol was frozen on 2026-08-01, before any real Level 3
prediction was ingested. The authoritative thresholds are machine-readable in
`experiments/preregistered/pre_wetlab_candidate_protocol.yaml`.

## Evidence boundary

Stage 0003A can nominate computational candidates at Level 1 or Level 2. It
cannot call a sequence functional, active, effective, or ready for biological
use. Real, independent monomer and RNA-complex prediction is needed for Level
3; experimental function is needed for Level 4. Fixture predictions exercise
only code paths and remain `is_mock=true` at Level 0.

## Sequential gates

1. Level 1 requires a standard amino-acid sequence, zero hard-fixed
   violations, Atlas nearest-neighbour identity at or below 0.80 using the
   registered MMseqs2 parameters, no homopolymer longer than seven residues,
   no newly introduced registered low-complexity window beyond the parent
   baseline, and composition deviation no greater than 0.25 from its parent.
2. Level 2 requires genuine outputs from at least two of the ESM-IF1,
   ProteinMPNN, and LigandMPNN model families. A LigandMPNN support claim is
   allowed only when its output audit proves that RNA atom context was retained.
3. A candidate cannot have a multi-state per-residue log-likelihood more than
   1.0 below its parent's worst registered native state, or state-score variance
   above 0.25. No severe single-state failure is permitted.
4. Candidates are separated into conservative (0.50–0.70 parent identity),
   moderate (0.30–0.50), and aggressive (0.18–0.30) novelty tiers. Boundaries
   are lower-inclusive and upper-exclusive.
5. Candidates are clustered at 70% identity with at least 80% coverage. At most
   one primary candidate is selected from a cluster within each novelty tier,
   preventing many near-duplicates from dominating a shortlist.

## Level 3 and ranking

The predeclared Level 3 gates are monomer TM-score at least 0.70, mean pLDDT at
least 70, interface PAE no greater than 12 Å, no severe state failure, and both
multi-seed and cross-model consistency. TM-score must come from US-align or
TM-align; no home-built proxy may be labeled TM-score.

Ranking uses Pareto dominance across sequence novelty, monomer structural
recovery, multi-state compatibility, RNA-interface preservation, model
agreement, and candidate diversity. ESM-IF1 conditional likelihood alone can
never determine the final ranking.

## Change control

Primary thresholds may not be changed after viewing real Level 3 results. A
necessary amendment must create a new versioned protocol with timestamp,
rationale, affected candidates, and explicit sensitivity-analysis status. The
original version and results under it remain in the audit trail.
