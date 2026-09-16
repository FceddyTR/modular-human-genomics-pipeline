from genomics_platform.ranking.benchmark.rank_aggregation_metrics import (
    evaluate_rank_aggregation,
)
from genomics_platform.ranking.benchmark.truth_set import (
    RankingTruth,
)
from genomics_platform.ranking.rank_aggregation import (
    RankAggregatedCandidate,
    RankAggregationMethod,
    RankAggregationResult,
)


COMPONENT_NAMES = (
    "clinical",
    "population",
    "functional",
    "gene_disease",
    "inheritance",
    "phenotype",
)


def _candidate(
    candidate_id,
    rank,
    score,
):
    return RankAggregatedCandidate(
        candidate_id=candidate_id,
        variant_key=f"1:{rank}:A:G",
        gene_symbol=f"GENE{rank}",
        rank=rank,
        aggregate_rank_score=score,
        components_participated=3,
        components_available=6,
        evidence_coverage=0.5,
        component_entries=(),
        warnings=(),
    )


def _result(*candidates):
    return RankAggregationResult(
        source_mode="phenotype_semantic",
        method=RankAggregationMethod.GEOMETRIC_MEAN,
        candidates=tuple(candidates),
        component_names=COMPONENT_NAMES,
    )


def _truth(case_id, *causal_ids):
    return RankingTruth(
        case_id=case_id,
        causal_candidate_ids=tuple(causal_ids),
    )


def test_rank_aggregation_metrics_unique_top_hit():
    result = _result(
        _candidate("causal", 1, 1.0),
        _candidate("decoy", 2, 2.0),
    )

    metrics = evaluate_rank_aggregation(
        result,
        _truth("case-1", "causal"),
        k_values=(1, 3, 5),
    )

    assert metrics.best_causal_rank == 1
    assert metrics.reciprocal_rank == 1.0
    assert metrics.hit_at_k[1] is True
    assert metrics.unique_hit_at_k[1] is True
    assert metrics.top_tie_size == 1
    assert metrics.causal_top_tie is False


def test_rank_aggregation_metrics_preserve_top_tie():
    result = _result(
        _candidate("causal", 1, 1.0),
        _candidate("decoy", 1, 1.0),
        _candidate("other", 3, 3.0),
    )

    metrics = evaluate_rank_aggregation(
        result,
        _truth("case-2", "causal"),
        k_values=(1, 3),
    )

    assert metrics.best_causal_rank == 1
    assert metrics.reciprocal_rank == 1.0
    assert metrics.hit_at_k[1] is True
    assert metrics.unique_hit_at_k[1] is False
    assert metrics.top_tie_size == 2
    assert metrics.causal_top_tie is True


def test_rank_aggregation_metrics_rank_two_causal():
    result = _result(
        _candidate("decoy", 1, 1.0),
        _candidate("causal", 2, 2.0),
    )

    metrics = evaluate_rank_aggregation(
        result,
        _truth("case-3", "causal"),
        k_values=(1, 3),
    )

    assert metrics.best_causal_rank == 2
    assert metrics.reciprocal_rank == 0.5
    assert metrics.hit_at_k[1] is False
    assert metrics.unique_hit_at_k[1] is False
    assert metrics.hit_at_k[3] is True
    assert metrics.unique_hit_at_k[3] is True


def test_rank_aggregation_metrics_multiple_causal_candidates():
    result = _result(
        _candidate("decoy", 1, 1.0),
        _candidate("causal-a", 2, 2.0),
        _candidate("causal-b", 3, 3.0),
    )

    metrics = evaluate_rank_aggregation(
        result,
        _truth(
            "case-4",
            "causal-a",
            "causal-b",
        ),
        k_values=(1, 3),
    )

    assert metrics.causal_ranks == (2, 3)
    assert metrics.best_causal_rank == 2
    assert metrics.reciprocal_rank == 0.5
    assert metrics.hit_at_k[1] is False
    assert metrics.hit_at_k[3] is True
