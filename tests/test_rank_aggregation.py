import math

from genomics_platform.ranking.rank_aggregation import (
    RankAggregationMethod,
    _average_competition_ranks,
    aggregate_component_ranks,
)
from genomics_platform.ranking.ranking_engine import (
    rank_candidates,
)
from genomics_platform.ranking.ranking_modes import (
    RankingMode,
)
from genomics_platform.ranking.benchmark.synthetic_factory import (
    make_gene_context,
    make_synthetic_candidate,
)


def test_average_component_ranks_preserve_ties():
    ranks = _average_competition_ranks(
        [
            ("a", 1.0),
            ("b", 1.0),
            ("c", 0.5),
        ]
    )

    assert ranks["a"] == 1.5
    assert ranks["b"] == 1.5
    assert ranks["c"] == 3.0


def test_geometric_rank_aggregation_is_deterministic():
    candidates = (
        make_synthetic_candidate(
            candidate_id="a",
            gene="A",
            position=90001,
            af=1e-6,
        ),
        make_synthetic_candidate(
            candidate_id="b",
            gene="B",
            position=90002,
            af=0.01,
        ),
    )

    weighted = rank_candidates(
        candidates,
        mode=RankingMode.VARIANT_FIRST,
    )

    first = aggregate_component_ranks(
        weighted
    )
    second = aggregate_component_ranks(
        weighted
    )

    assert first == second


def test_missing_population_abstains():
    candidates = (
        make_synthetic_candidate(
            candidate_id="measured",
            gene="MEASURED",
            position=90011,
            af=1e-6,
        ),
        make_synthetic_candidate(
            candidate_id="missing",
            gene="MISSING",
            position=90012,
            af=None,
        ),
    )

    weighted = rank_candidates(
        candidates,
        mode=RankingMode.VARIANT_FIRST,
    )

    result = aggregate_component_ranks(
        weighted
    )

    missing = next(
        c
        for c in result.candidates
        if c.candidate_id == "missing"
    )

    population = next(
        e
        for e in missing.component_entries
        if e.component_name == "population"
    )

    assert population.abstained is True
    assert population.component_rank is None
    assert population.normalized_evidence is None


def test_population_error_does_not_become_support():
    candidates = (
        make_synthetic_candidate(
            candidate_id="error",
            gene="ERROR",
            position=90021,
            af=None,
            population_status="ERROR",
        ),
        make_synthetic_candidate(
            candidate_id="measured",
            gene="MEASURED",
            position=90022,
            af=1e-6,
        ),
    )

    weighted = rank_candidates(
        candidates,
        mode=RankingMode.VARIANT_FIRST,
    )

    result = aggregate_component_ranks(
        weighted
    )

    error = next(
        c
        for c in result.candidates
        if c.candidate_id == "error"
    )

    population = next(
        e
        for e in error.component_entries
        if e.component_name == "population"
    )

    assert population.status == "ERROR"
    assert population.abstained is True
    assert population.normalized_evidence is None


def test_gene_context_can_contribute_without_singleton_zero_points():
    context = make_gene_context(
        gene="CAUSAL",
        gencc_classification="Definitive",
        panelapp_confidence="Green",
    )

    candidates = (
        make_synthetic_candidate(
            candidate_id="causal",
            gene="CAUSAL",
            position=90031,
            gene_context=context,
        ),
        make_synthetic_candidate(
            candidate_id="decoy",
            gene="DECOY",
            position=90032,
        ),
    )

    weighted = rank_candidates(
        candidates,
        mode=RankingMode.PHENOTYPE,
    )

    result = aggregate_component_ranks(
        weighted
    )

    causal = next(
        c
        for c in result.candidates
        if c.candidate_id == "causal"
    )

    gene_entry = next(
        e
        for e in causal.component_entries
        if e.component_name == "gene_disease"
    )

    assert gene_entry.abstained is False
    assert gene_entry.component_rank == 1.0


def test_arithmetic_method_is_supported():
    candidates = (
        make_synthetic_candidate(
            candidate_id="a",
            gene="A",
            position=90041,
        ),
        make_synthetic_candidate(
            candidate_id="b",
            gene="B",
            position=90042,
        ),
    )

    weighted = rank_candidates(
        candidates,
        mode=RankingMode.VARIANT_FIRST,
    )

    result = aggregate_component_ranks(
        weighted,
        method=(
            RankAggregationMethod.ARITHMETIC_MEAN
        ),
    )

    assert (
        result.method
        == RankAggregationMethod.ARITHMETIC_MEAN
    )

    assert all(
        math.isfinite(
            candidate.aggregate_rank_score
        )
        for candidate in result.candidates
    )


def test_result_declares_semantic_boundaries():
    candidate = make_synthetic_candidate(
        candidate_id="candidate",
        gene="GENE",
        position=90051,
    )

    weighted = rank_candidates(
        (candidate,),
        mode=RankingMode.VARIANT_FIRST,
    )

    result = aggregate_component_ranks(
        weighted
    ).to_dict()

    boundaries = result[
        "semantic_boundaries"
    ]

    assert (
        boundaries[
            "benchmark_truth_used"
        ]
        is False
    )
    assert (
        boundaries[
            "pathogenicity_probability"
        ]
        is False
    )
    assert (
        boundaries[
            "acmg_classification"
        ]
        is False
    )
    assert boundaries["diagnosis"] is False
