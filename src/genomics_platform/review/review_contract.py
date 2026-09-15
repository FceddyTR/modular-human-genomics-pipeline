"""Contracts for post-ranking technical variant review.

Technical review is deliberately separated from:
- evidence aggregation,
- candidate prioritization,
- pathogenicity assessment,
- ACMG/AMP classification,
- diagnosis.

Review observations never modify ranking scores.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from enum import Enum
from typing import Any


SCHEMA_VERSION = "1.0"


class ReviewStatus(str, Enum):
    NOT_REVIEWED = "NOT_REVIEWED"
    SUPPORTS_CALL = "SUPPORTS_CALL"
    AMBIGUOUS = "AMBIGUOUS"
    SUSPECTED_ARTIFACT = "SUSPECTED_ARTIFACT"


@dataclass(frozen=True)
class TechnicalReviewRecord:
    """Human technical review of one ranked variant."""

    candidate_id: str
    variant_key: str
    status: ReviewStatus = ReviewStatus.NOT_REVIEWED
    reviewer: str | None = None
    reviewed_at: str | None = None
    notes: str | None = None
    schema_version: str = SCHEMA_VERSION

    def __post_init__(self) -> None:
        candidate_id = self.candidate_id.strip()
        variant_key = self.variant_key.strip()

        if not candidate_id:
            raise ValueError(
                "candidate_id is required."
            )

        if not variant_key:
            raise ValueError(
                "variant_key is required."
            )

        object.__setattr__(
            self,
            "candidate_id",
            candidate_id,
        )

        object.__setattr__(
            self,
            "variant_key",
            variant_key,
        )

        reviewer = (
            self.reviewer.strip()
            if self.reviewer
            else None
        )

        notes = (
            self.notes.strip()
            if self.notes
            else None
        )

        object.__setattr__(
            self,
            "reviewer",
            reviewer,
        )

        object.__setattr__(
            self,
            "notes",
            notes,
        )

        if (
            self.status
            is ReviewStatus.NOT_REVIEWED
            and self.reviewed_at is not None
        ):
            raise ValueError(
                "NOT_REVIEWED cannot have reviewed_at."
            )

    @classmethod
    def reviewed(
        cls,
        *,
        candidate_id: str,
        variant_key: str,
        status: ReviewStatus,
        reviewer: str | None = None,
        notes: str | None = None,
    ) -> "TechnicalReviewRecord":
        if status is ReviewStatus.NOT_REVIEWED:
            raise ValueError(
                "reviewed() requires a completed "
                "review status."
            )

        timestamp = (
            datetime.now(timezone.utc)
            .isoformat()
        )

        return cls(
            candidate_id=candidate_id,
            variant_key=variant_key,
            status=status,
            reviewer=reviewer,
            reviewed_at=timestamp,
            notes=notes,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "candidate_id": self.candidate_id,
            "variant_key": self.variant_key,
            "status": self.status.value,
            "reviewer": self.reviewer,
            "reviewed_at": self.reviewed_at,
            "notes": self.notes,
            "semantic_boundaries": {
                "post_ranking_review": True,
                "technical_review_only": True,
                "affects_ranking": False,
                "pathogenicity_probability": False,
                "acmg_amp_classification": False,
                "variant_classification": False,
                "diagnosis": False,
            },
        }
