"""Backend interface shared by genuine and test-only inverse-folding engines."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

from cas13_if.schemas import (
    BackendCapabilities,
    Candidate,
    SampleRequest,
    ScoreRequest,
    ScoreResult,
)


class InverseFoldingBackend(ABC):
    @abstractmethod
    def capabilities(self) -> BackendCapabilities:
        """Describe supported operations and conditioning context."""

    @abstractmethod
    def load(self) -> None:
        """Load local assets without silently using the network."""

    @abstractmethod
    def score(self, request: ScoreRequest) -> ScoreResult:
        """Score a sequence conditioned on a supplied structure."""

    @abstractmethod
    def sample(self, request: SampleRequest) -> list[Candidate]:
        """Sample validated candidates."""

    @abstractmethod
    def metadata(self) -> dict[str, Any]:
        """Return backend, checkpoint, license, and runtime metadata."""
