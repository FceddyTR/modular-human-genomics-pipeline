"""Stable contracts for comparative candidate ranking."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from genomics_platform.ranking.ranking_modes import (
    RankingMode,
)


SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class ComponentContribution:
    name: str
    contribution: float
    max_contribution: float
    status: str
    reasons: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.max_contribution < 0:
            raise ValueError(
                "max_contribution cannot be negative."
            )

        if self.contribution < 0:
            raise ValueError(
                "contribution cannot be negative."
            )

        if self.contribution > self.max_contribution:
            raise ValueError(
                "contribution cannot exceed "
                "max_contribution."
            )

        allowed = {
            "SUPPORTING",
            "NEUTRAL",
            "UNAVAILABLE",
            "BLOCKED",
            "ERROR",
        }

        if self.status not in allowed:
            raise ValueError(
                f"Invalid component status: "
                f"{self.status!r}"
            )

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "contribution": self.contribution,
            "max_contribution": self.max_contribution,
            "status": self.status,
            "reasons": list(self.reasons),
            "limitations": list(self.limitations),
        }


@dataclass(frozen=True)
class RankedCandidate:
    candidate_id: str
    variant_key: str
    gene_symbol: str | None
    rank: int
    prioritization_score: float
    components: tuple[
        ComponentContribution,
        ...
    ]
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.rank < 1:
            raise ValueError(
                "rank must be >= 1."
            )

        if self.prioritization_score < 0:
            raise ValueError(
                "prioritization_score cannot "
                "be negative."
            )

    def to_dict(self) -> dict[str, Any]:
        return {
            "candidate_id": self.candidate_id,
            "variant_key": self.variant_key,
            "gene_symbol": self.gene_symbol,
            "rank": self.rank,
            "prioritization_score": (
                self.prioritization_score
            ),
            "components": [
                component.to_dict()
                for component in self.components
            ],
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class RankingResult:
    mode: RankingMode
    candidates: tuple[RankedCandidate, ...]
    engine_version: str
    schema_version: str = SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "engine": {
                "name": "candidate_prioritization",
                "version": self.engine_version,
            },
            "mode": self.mode.value,
            "candidate_count": len(
                self.candidates
            ),
            "candidates": [
                candidate.to_dict()
                for candidate in self.candidates
            ],
            "semantic_boundaries": {
                "relative_prioritization_only": True,
                "pathogenicity_probability": False,
                "acmg_amp_classification": False,
                "variant_classification": False,
                "diagnosis": False,
                "candidate_tier": False,
            },
        }
