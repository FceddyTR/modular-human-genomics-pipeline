import pytest

from genomics_platform.ranking.ranking_modes import (
    RankingMode,
    normalize_ranking_mode,
)
from genomics_platform.ranking.ranking_contract import (
    RankedCandidate,
    RankingResult,
)
from genomics_platform.ranking.benchmark.truth_set import (
    RankingTruth,
)
from genomics_platform.ranking.benchmark.metrics import (
    evaluate_ranking,
)
from genomics_platform.ranking.benchmark.comparison import (
    RankingComparison,
)


def candidate(candidate_id, rank):
    return RankedCandidate(
        candidate_id=candidate_id,
        variant_key=f"1:{rank}:A:T",
        gene_symbol=f"GENE{rank}",
        rank=rank,
        prioritization_score=float(
            100 - rank
        ),
        components=(),
    )


def result(mode, ids):
    return RankingResult(
        mode=mode,
        candidates=tuple(
            candidate(
                candidate_id,
                rank,
            )
            for rank, candidate_id
            in enumerate(
                ids,
                start=1,
            )
        ),
        engine_version="test",
    )


def test_modes():
    assert (
        normalize_ranking_mode(
            "variant_first"
        )
        == RankingMode.VARIANT_FIRST
    )

    assert (
        normalize_ranking_mode(
            "phenotype"
        )
        == RankingMode.PHENOTYPE
    )


def test_invalid_mode():
    with pytest.raises(ValueError):
        normalize_ranking_mode(
            "diagnosis"
        )


def test_truth_requires_candidate():
    with pytest.raises(ValueError):
        RankingTruth(
            case_id="case",
            causal_candidate_ids=(),
        )


def test_hit_at_k():
    ranking = result(
        RankingMode.VARIANT_FIRST,
        [
            "decoy1",
            "decoy2",
            "causal",
            "decoy3",
        ],
    )

    truth = RankingTruth(
        case_id="case-1",
        causal_candidate_ids=(
            "causal",
        ),
    )

    metrics = evaluate_ranking(
        ranking,
        truth,
        k_values=(1, 3, 5),
    )

    assert metrics.hit_at_k[1] is False
    assert metrics.hit_at_k[3] is True
    assert metrics.hit_at_k[5] is True


def test_reciprocal_rank():
    ranking = result(
        RankingMode.VARIANT_FIRST,
        [
            "a",
            "b",
            "causal",
        ],
    )

    truth = RankingTruth(
        case_id="case",
        causal_candidate_ids=(
            "causal",
        ),
    )

    metrics = evaluate_ranking(
        ranking,
        truth,
        k_values=(3,),
    )

    assert metrics.best_causal_rank == 3
    assert metrics.reciprocal_rank == pytest.approx(
        1 / 3
    )


def test_precision_recall_at_k():
    ranking = result(
        RankingMode.VARIANT_FIRST,
        [
            "causal1",
            "decoy",
            "causal2",
            "decoy2",
        ],
    )

    truth = RankingTruth(
        case_id="case",
        causal_candidate_ids=(
            "causal1",
            "causal2",
        ),
    )

    metrics = evaluate_ranking(
        ranking,
        truth,
        k_values=(1, 3),
    )

    assert (
        metrics.precision_at_k[1]
        == 1.0
    )

    assert (
        metrics.recall_at_k[1]
        == 0.5
    )

    assert (
        metrics.precision_at_k[3]
        == pytest.approx(2 / 3)
    )

    assert (
        metrics.recall_at_k[3]
        == 1.0
    )


def test_phenotype_rank_improvement():
    truth = RankingTruth(
        case_id="case",
        causal_candidate_ids=(
            "causal",
        ),
    )

    baseline = evaluate_ranking(
        result(
            RankingMode.VARIANT_FIRST,
            [
                "a",
                "b",
                "c",
                "causal",
            ],
        ),
        truth,
        k_values=(1, 5),
    )

    phenotype = evaluate_ranking(
        result(
            RankingMode.PHENOTYPE,
            [
                "causal",
                "a",
                "b",
                "c",
            ],
        ),
        truth,
        k_values=(1, 5),
    )

    comparison = RankingComparison(
        baseline=baseline,
        phenotype=phenotype,
    )

    assert (
        comparison.causal_rank_improvement
        == 3
    )

    assert (
        baseline.hit_at_k[1]
        is False
    )

    assert (
        phenotype.hit_at_k[1]
        is True
    )


def test_truth_not_present():
    truth = RankingTruth(
        case_id="case",
        causal_candidate_ids=(
            "missing",
        ),
    )

    metrics = evaluate_ranking(
        result(
            RankingMode.VARIANT_FIRST,
            [
                "a",
                "b",
            ],
        ),
        truth,
        k_values=(1, 5),
    )

    assert metrics.best_causal_rank is None
    assert metrics.reciprocal_rank == 0.0
    assert metrics.hit_at_k[5] is False


def test_result_semantic_boundaries():
    data = result(
        RankingMode.PHENOTYPE,
        ["candidate"],
    ).to_dict()

    boundaries = data[
        "semantic_boundaries"
    ]

    assert (
        boundaries[
            "relative_prioritization_only"
        ]
        is True
    )

    assert (
        boundaries[
            "pathogenicity_probability"
        ]
        is False
    )

    assert boundaries["diagnosis"] is False
