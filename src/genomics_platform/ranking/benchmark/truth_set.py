"""Truth labels used only after candidate ranking."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class RankingTruth:
    case_id: str
    causal_candidate_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        case_id = self.case_id.strip()

        if not case_id:
            raise ValueError(
                "case_id is required."
            )

        causal = tuple(
            dict.fromkeys(
                candidate.strip()
                for candidate
                in self.causal_candidate_ids
                if candidate.strip()
            )
        )

        if not causal:
            raise ValueError(
                "At least one causal candidate "
                "is required."
            )

        object.__setattr__(
            self,
            "case_id",
            case_id,
        )

        object.__setattr__(
            self,
            "causal_candidate_ids",
            causal,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "case_id": self.case_id,
            "causal_candidate_ids": list(
                self.causal_candidate_ids
            ),
        }
