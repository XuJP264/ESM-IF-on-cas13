"""Causal constrained autoregressive sampling.

Fixed residues are inserted at their decoding step and become prefix context for
later residues. Earlier residues cannot observe future fixed tokens. This is not
bidirectional global conditional sampling.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cas13_if.schemas import STANDARD_AA

LogitsFunction = Callable[[str, int], NDArray[np.float64]]


@dataclass(frozen=True)
class DecodedPosition:
    index: int
    logits: tuple[float, ...]
    probabilities: tuple[float, ...]
    selected_token: str
    fixed: bool
    temperature: float
    seed: int


@dataclass(frozen=True)
class DecodedSequence:
    sequence: str
    trace: tuple[DecodedPosition, ...]
    fixed_position_violations: int
    semantics: str = "left_to_right_causal_hard_fixed"


def constrained_autoregressive_sample(
    *,
    length: int,
    logits_function: LogitsFunction,
    alphabet: Sequence[str] = tuple(sorted(STANDARD_AA)),
    fixed_positions: dict[int, str] | None = None,
    allowed_residues: dict[int, set[str]] | None = None,
    temperature: float = 1.0,
    seed: int = 20260731,
) -> DecodedSequence:
    if length < 1:
        raise ValueError("length must be positive")
    if temperature <= 0:
        raise ValueError("temperature must be positive")
    tokens = tuple(token.upper() for token in alphabet)
    if len(tokens) != len(set(tokens)) or not set(tokens).issubset(STANDARD_AA):
        raise ValueError("alphabet must contain unique standard amino acids")
    token_to_index = {token: index for index, token in enumerate(tokens)}
    fixed = {index: token.upper() for index, token in (fixed_positions or {}).items()}
    allowed = {
        index: {token.upper() for token in values}
        for index, values in (allowed_residues or {}).items()
    }
    _validate_constraints(length, token_to_index, fixed, allowed)
    rng = np.random.default_rng(seed)
    prefix: list[str] = []
    trace: list[DecodedPosition] = []
    for index in range(length):
        logits = np.asarray(logits_function("".join(prefix), index), dtype=np.float64)
        if logits.shape != (len(tokens),):
            raise ValueError(
                f"logits at position {index} have shape {logits.shape}; "
                f"expected {(len(tokens),)}"
            )
        if not np.all(np.isfinite(logits)):
            raise ValueError(f"non-finite logits at position {index}")
        filtered = logits.copy()
        if index in allowed:
            disallowed = [
                token_index
                for token_index, token in enumerate(tokens)
                if token not in allowed[index]
            ]
            filtered[disallowed] = -np.inf
        probabilities = _softmax(filtered / temperature)
        is_fixed = index in fixed
        if is_fixed:
            selected = fixed[index]
        else:
            selected = str(rng.choice(tokens, p=probabilities))
        prefix.append(selected)
        trace.append(
            DecodedPosition(
                index=index,
                logits=tuple(float(value) for value in logits),
                probabilities=tuple(float(value) for value in probabilities),
                selected_token=selected,
                fixed=is_fixed,
                temperature=temperature,
                seed=seed,
            )
        )
    sequence = "".join(prefix)
    violations = sum(sequence[index] != token for index, token in fixed.items())
    return DecodedSequence(
        sequence=sequence,
        trace=tuple(trace),
        fixed_position_violations=violations,
    )


def _validate_constraints(
    length: int,
    token_to_index: dict[str, int],
    fixed: dict[int, str],
    allowed: dict[int, set[str]],
) -> None:
    for index, token in fixed.items():
        if not 0 <= index < length:
            raise ValueError(f"fixed position {index} outside [0, {length})")
        if token not in token_to_index:
            raise ValueError(f"fixed token {token!r} is not in the alphabet")
        if index in allowed and token not in allowed[index]:
            raise ValueError(f"fixed token {token!r} is disallowed at position {index}")
    for index, tokens in allowed.items():
        if not 0 <= index < length:
            raise ValueError(f"allowed-residue position {index} outside range")
        if not tokens:
            raise ValueError(f"empty allowed-residue set at position {index}")
        unknown = tokens.difference(token_to_index)
        if unknown:
            raise ValueError(f"unknown allowed tokens at {index}: {sorted(unknown)}")


def _softmax(values: NDArray[np.float64]) -> NDArray[np.float64]:
    maximum = float(np.max(values))
    if not np.isfinite(maximum):
        raise ValueError("all amino acids were filtered out")
    exponential = np.exp(values - maximum)
    total = float(exponential.sum())
    if not np.isfinite(total) or total <= 0:
        raise ValueError("invalid probability normalization")
    return exponential / total
