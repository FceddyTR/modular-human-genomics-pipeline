"""Evaluation metrics for candidate ranking."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

from genomics_platform.ranking.ranking_contract import (
    RankingResult,
)
from genomics_platform.ranking.benchmark.truth_set import (
    RankingTruth,
)


DEFAULT_K_VALUES = (
    1,
    3,
    5,
    10,
    20,
)


@dataclass(frozen=True)
class RankingMetrics:
    case_id: str
    causal_ranks: tuple[int, ...]
    best_causal_rank: int | None
    reciprocal_rank: float
    hit_at_k: dict[int, bool]
    precision_at_k: dict[int, float]
    recall_at_k: dict[int, float]

    def to_dict(self) -> dict[str, Any]:
        return {
            "case_id": self.case_id,
            "causal_ranks": list(
                self.causal_ranks
            ),
            "best_causal_rank": (
                self.best_causal_rank
            ),
            "reciprocal_rank": (
                self.reciprocal_rank
            ),
            "hit_at_k": {
                str(k): value
                for k, value
                in self.hit_at_k.items()
            },
            "precision_at_k": {
                str(k): value
                for k, value
                in self.precision_at_k.items()
            },
            "recall_at_k": {
                str(k): value
                for k, value
                in self.recall_at_k.items()
            },
        }


def evaluate_ranking(
    result: RankingResult,
    truth: RankingTruth,
    k_values: Iterable[int] = DEFAULT_K_VALUES,
) -> RankingMetrics:
    k_values = tuple(
        sorted(
            set(
                int(k)
                for k in k_values
            )
        )
    )

    if any(k < 1 for k in k_values):
        raise ValueError(
            "All k values must be >= 1."
        )

    truth_ids = set(
        truth.causal_candidate_ids
    )

    ranked_ids = [
        candidate.candidate_id
        for candidate in result.candidates
    ]

    causal_ranks = tuple(
        index
        for index, candidate_id
        in enumerate(
            ranked_ids,
            start=1,
        )
        if candidate_id in truth_ids
    )

    best_rank = (
        min(causal_ranks)
        if causal_ranks
        else None
    )

    reciprocal_rank = (
        1.0 / best_rank
        if best_rank is not None
        else 0.0
    )

    hit_at_k = {}
    precision_at_k = {}
    recall_at_k = {}

    total_relevant = len(truth_ids)

    for k in k_values:
        top_k = ranked_ids[:k]

        relevant_retrieved = sum(
            candidate_id in truth_ids
            for candidate_id in top_k
        )

        hit_at_k[k] = (
            relevant_retrieved > 0
        )

        # Precision denominator is the actual
        # number returned when fewer than k
        # candidates exist.
        denominator = len(top_k)

        precision_at_k[k] = (
            relevant_retrieved / denominator
            if denominator
            else 0.0
        )

        recall_at_k[k] = (
            relevant_retrieved
            / total_relevant
        )

    return RankingMetrics(
        case_id=truth.case_id,
        causal_ranks=causal_ranks,
        best_causal_rank=best_rank,
        reciprocal_rank=reciprocal_rank,
        hit_at_k=hit_at_k,
        precision_at_k=precision_at_k,
        recall_at_k=recall_at_k,
    )
