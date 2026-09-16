"""Tie-aware benchmark metrics for experimental rank aggregation."""

from __future__ import annotations

from dataclasses import dataclass

from genomics_platform.ranking.rank_aggregation import (
    RankAggregationResult,
)
from genomics_platform.ranking.benchmark.truth_set import (
    RankingTruth,
)


@dataclass(frozen=True)
class RankAggregationMetrics:
    case_id: str
    causal_ranks: tuple[int, ...]
    best_causal_rank: int | None
    reciprocal_rank: float
    hit_at_k: dict[int, bool]
    unique_hit_at_k: dict[int, bool]
    top_tie_size: int
    causal_top_tie: bool


def evaluate_rank_aggregation(
    result: RankAggregationResult,
    truth: RankingTruth,
    *,
    k_values=(1, 3, 5),
) -> RankAggregationMetrics:
    k_values = tuple(
        sorted(set(int(k) for k in k_values))
    )

    if any(k < 1 for k in k_values):
        raise ValueError("All k values must be >= 1.")

    causal_ids = set(
        truth.causal_candidate_ids
    )

    causal_ranks = tuple(
        sorted(
            candidate.rank
            for candidate in result.candidates
            if candidate.candidate_id in causal_ids
        )
    )

    best_causal_rank = (
        causal_ranks[0]
        if causal_ranks
        else None
    )

    reciprocal_rank = (
        1.0 / best_causal_rank
        if best_causal_rank is not None
        else 0.0
    )

    top_rank = (
        min(
            candidate.rank
            for candidate in result.candidates
        )
        if result.candidates
        else None
    )

    top_candidates = (
        [
            candidate
            for candidate in result.candidates
            if candidate.rank == top_rank
        ]
        if top_rank is not None
        else []
    )

    top_tie_size = len(top_candidates)

    causal_top_tie = (
        top_tie_size > 1
        and any(
            candidate.candidate_id in causal_ids
            for candidate in top_candidates
        )
    )

    hit_at_k = {}
    unique_hit_at_k = {}

    for k in k_values:
        hit = (
            best_causal_rank is not None
            and best_causal_rank <= k
        )

        hit_at_k[k] = hit

        if not hit:
            unique_hit_at_k[k] = False
            continue

        causal_at_best_rank = {
            candidate.candidate_id
            for candidate in result.candidates
            if (
                candidate.rank == best_causal_rank
                and candidate.candidate_id
                in causal_ids
            )
        }

        noncausal_at_best_rank = {
            candidate.candidate_id
            for candidate in result.candidates
            if (
                candidate.rank == best_causal_rank
                and candidate.candidate_id
                not in causal_ids
            )
        }

        unique_hit_at_k[k] = bool(
            causal_at_best_rank
            and not noncausal_at_best_rank
        )

    return RankAggregationMetrics(
        case_id=truth.case_id,
        causal_ranks=causal_ranks,
        best_causal_rank=best_causal_rank,
        reciprocal_rank=reciprocal_rank,
        hit_at_k=hit_at_k,
        unique_hit_at_k=unique_hit_at_k,
        top_tie_size=top_tie_size,
        causal_top_tie=causal_top_tie,
    )
