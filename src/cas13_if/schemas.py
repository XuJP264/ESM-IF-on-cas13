"""Shared validated schemas for backends, candidates, scores, and traces."""

from __future__ import annotations

from enum import IntEnum
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


def normalize_protein_sequence(value: str | None) -> str | None:
    if value is None:
        return value
    sequence = "".join(value.split()).upper()
    invalid = sorted(set(sequence).difference(STANDARD_AA))
    if invalid:
        raise ValueError(f"sequence contains non-standard residues: {invalid}")
    return sequence


class EvidenceLevel(IntEnum):
    IO_VALIDATED = 0
    STATISTICAL_NOVELTY = 1
    INVERSE_FOLDING_COMPATIBILITY = 2
    MULTI_MODEL_PLAUSIBILITY = 3
    WET_LAB_VALIDATED = 4


class BackendCapabilities(BaseModel):
    model_config = ConfigDict(extra="forbid")

    scoring: bool
    sampling: bool
    protein_multichain: bool = False
    rna_atomic_context: bool = False
    hard_fixed: bool = False
    allowed_residue_filter: bool = False
    per_residue_probabilities: bool = False


class PositionTrace(BaseModel):
    model_config = ConfigDict(extra="forbid")

    index: int = Field(ge=0)
    logits: list[float]
    probabilities: list[float]
    selected_token: str = Field(min_length=1, max_length=1)
    fixed: bool
    temperature: float = Field(gt=0)
    seed: int

    @field_validator("selected_token")
    @classmethod
    def validate_token(cls, value: str) -> str:
        value = value.upper()
        if value not in STANDARD_AA:
            raise ValueError(f"non-standard amino acid token: {value}")
        return value


class Candidate(BaseModel):
    model_config = ConfigDict(extra="forbid")

    candidate_id: str = Field(min_length=1)
    scaffold_id: str = Field(min_length=1)
    backend: str = Field(min_length=1)
    sequence: str = Field(min_length=1)
    parent_sequence: str | None = None
    seed: int
    temperature: float = Field(gt=0)
    is_mock: bool
    evidence_level: EvidenceLevel
    fixed_positions: dict[int, str] = Field(default_factory=dict)
    traces: list[PositionTrace] = Field(default_factory=list)
    metadata: dict[str, Any] = Field(default_factory=dict)

    @field_validator("sequence", "parent_sequence")
    @classmethod
    def validate_sequence(cls, value: str | None) -> str | None:
        return normalize_protein_sequence(value)

    @field_validator("fixed_positions")
    @classmethod
    def validate_fixed_positions(cls, value: dict[int, str]) -> dict[int, str]:
        normalized: dict[int, str] = {}
        for index, amino_acid in value.items():
            token = amino_acid.upper()
            if index < 0 or token not in STANDARD_AA:
                raise ValueError(f"invalid fixed position {index}: {amino_acid}")
            normalized[index] = token
        return normalized

    @model_validator(mode="after")
    def fixed_positions_are_preserved(self) -> Candidate:
        for index, amino_acid in self.fixed_positions.items():
            if index >= len(self.sequence):
                raise ValueError(f"fixed position {index} exceeds sequence length")
            if self.sequence[index] != amino_acid:
                raise ValueError(
                    f"fixed-position violation at {index}: "
                    f"{self.sequence[index]} != {amino_acid}"
                )
        if self.is_mock and self.evidence_level > EvidenceLevel.IO_VALIDATED:
            raise ValueError("mock candidates cannot support evidence above Level 0")
        return self


class ScoreRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    scaffold_id: str
    structure_path: str
    sequence: str
    protein_chains: list[str] = Field(default_factory=list)
    seed: int = 20260731

    _validate_sequence = field_validator("sequence")(normalize_protein_sequence)


class ScoreResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    scaffold_id: str
    backend: str
    sequence: str
    conditional_log_likelihood: float
    perplexity: float = Field(gt=0)
    per_residue_log_probabilities: list[float]
    is_mock: bool
    evidence_level: EvidenceLevel
    metadata: dict[str, Any] = Field(default_factory=dict)

    @model_validator(mode="after")
    def lengths_and_evidence_are_consistent(self) -> ScoreResult:
        if len(self.sequence) != len(self.per_residue_log_probabilities):
            raise ValueError("per-residue scores must match sequence length")
        if self.is_mock and self.evidence_level > EvidenceLevel.IO_VALIDATED:
            raise ValueError("mock scores cannot support evidence above Level 0")
        return self


class SampleRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    scaffold_id: str
    structure_path: str
    parent_sequence: str
    count: int = Field(default=1, ge=1)
    temperature: float = Field(default=1.0, gt=0)
    seed: int = 20260731
    fixed_positions: dict[int, str] = Field(default_factory=dict)
    allowed_residues: dict[int, set[str]] = Field(default_factory=dict)
    protein_chains: list[str] = Field(default_factory=list)

    _validate_parent = field_validator("parent_sequence")(normalize_protein_sequence)

    @model_validator(mode="after")
    def positions_are_valid(self) -> SampleRequest:
        length = len(self.parent_sequence)
        for index, amino_acid in self.fixed_positions.items():
            if index < 0 or index >= length:
                raise ValueError(f"fixed position {index} outside [0, {length})")
            if amino_acid.upper() not in STANDARD_AA:
                raise ValueError(f"invalid fixed amino acid: {amino_acid}")
        for index, allowed in self.allowed_residues.items():
            if index < 0 or index >= length:
                raise ValueError(f"allowed-residue position {index} outside range")
            if not allowed or not set(allowed).issubset(STANDARD_AA):
                raise ValueError(f"invalid allowed-residue set at {index}")
        return self


class ConstraintPosition(BaseModel):
    model_config = ConfigDict(extra="forbid")

    index_0: int = Field(ge=0)
    biological_index_1: int = Field(ge=1)
    pdb_residue_number: int | None = None
    insertion_code: str = ""
    wt_amino_acid: str = Field(min_length=1, max_length=1)
    subtype: str | None = None
    domain: str | None = None
    structural_state: str | None = None
    solvent_accessibility: float | None = None
    burial: Literal["buried", "surface", "unknown"] = "unknown"
    rna_contact_distance: float | None = None
    target_rna_contact: bool = False
    crrna_contact: bool = False
    conservation: float | None = None
    entropy: float | None = None
    sequence_weight: float | None = None
    coevolution_score: float | None = None
    hepn_motif_annotation: str | None = None
    interface_annotation: str | None = None
    manual_annotation: str | None = None
    final_class: Literal["hard_fixed", "soft_constrained", "free"]
    decision_reasons: list[str]
