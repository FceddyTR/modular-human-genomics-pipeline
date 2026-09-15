"""Tie-aware benchmark metrics for Borda candidate ranking."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

from genomics_platform.ranking.borda import (
    BordaResult,
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
class BordaRankingMetrics:
    case_id: str
    causal_ranks: tuple[int, ...]
    best_causal_rank: int | None
    reciprocal_rank: float
    hit_at_k: dict[int, bool]
    unique_hit_at_k: dict[int, bool]
    top_tie_size: int
    causal_top_tie: bool

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
            "unique_hit_at_k": {
                str(k): value
                for k, value
                in self.unique_hit_at_k.items()
            },
            "top_tie_size": self.top_tie_size,
            "causal_top_tie": (
                self.causal_top_tie
            ),
            "metric_semantics": {
                "competition_rank": True,
                "ties_preserved": True,
                "reciprocal_rank_uses_competition_rank": True,
                "unique_hit_requires_no_noncausal_tie": True,
            },
        }


def evaluate_borda_ranking(
    result: BordaResult,
    truth: RankingTruth,
    k_values: Iterable[int] = DEFAULT_K_VALUES,
) -> BordaRankingMetrics:
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

    causal_ranks = tuple(
        candidate.rank
        for candidate in result.candidates
        if candidate.candidate_id
        in truth_ids
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
    unique_hit_at_k = {}

    for k in k_values:
        eligible = tuple(
            candidate
            for candidate in result.candidates
            if candidate.rank <= k
        )

        relevant = tuple(
            candidate
            for candidate in eligible
            if candidate.candidate_id
            in truth_ids
        )

        hit_at_k[k] = bool(
            relevant
        )

        # A unique hit means at least one causal
        # candidate is within the requested rank
        # boundary and there is no non-causal
        # candidate sharing the best causal rank.
        if not relevant:
            unique_hit_at_k[k] = False
            continue

        best_relevant_rank = min(
            candidate.rank
            for candidate in relevant
        )

        tied_at_best = tuple(
            candidate
            for candidate in result.candidates
            if candidate.rank
            == best_relevant_rank
        )

        unique_hit_at_k[k] = all(
            candidate.candidate_id
            in truth_ids
            for candidate in tied_at_best
        )

    top_candidates = tuple(
        candidate
        for candidate in result.candidates
        if candidate.rank == 1
    )

    top_tie_size = len(
        top_candidates
    )

    causal_top_tie = (
        top_tie_size > 1
        and any(
            candidate.candidate_id
            in truth_ids
            for candidate in top_candidates
        )
    )

    return BordaRankingMetrics(
        case_id=truth.case_id,
        causal_ranks=causal_ranks,
        best_causal_rank=best_rank,
        reciprocal_rank=reciprocal_rank,
        hit_at_k=hit_at_k,
        unique_hit_at_k=unique_hit_at_k,
        top_tie_size=top_tie_size,
        causal_top_tie=causal_top_tie,
    )
