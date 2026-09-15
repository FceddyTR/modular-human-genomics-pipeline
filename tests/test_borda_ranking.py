"""Tests for evidence-component Borda rank aggregation."""

from __future__ import annotations

import pytest

from genomics_platform.ranking.borda import (
    aggregate_borda,
)
from genomics_platform.ranking.ranking_contract import (
    ComponentContribution,
    RankedCandidate,
    RankingResult,
)
from genomics_platform.ranking.ranking_modes import (
    RankingMode,
)


def component(
    name,
    contribution,
    maximum=10.0,
    status="SUPPORTING",
):
    return ComponentContribution(
        name=name,
        contribution=contribution,
        max_contribution=maximum,
        status=status,
    )


def candidate(
    candidate_id,
    components,
):
    return RankedCandidate(
        candidate_id=candidate_id,
        variant_key=f"1:100:A:{candidate_id}",
        gene_symbol=candidate_id.upper(),
        rank=1,
        prioritization_score=0.0,
        components=tuple(components),
    )


def result(*candidates):
    return RankingResult(
        mode=RankingMode.PHENOTYPE_SEMANTIC,
        candidates=tuple(candidates),
        engine_version="synthetic",
    )


def entry(
    borda_candidate,
    component_name,
):
    return next(
        item
        for item in borda_candidate.ballot_entries
        if item.component_name == component_name
    )


def test_single_ballot_orders_candidates():
    source = result(
        candidate(
            "a",
            [component("clinical", 10)],
        ),
        candidate(
            "b",
            [component("clinical", 5)],
        ),
        candidate(
            "c",
            [component("clinical", 0, status="NEUTRAL")],
        ),
    )

    borda = aggregate_borda(source)

    assert [
        item.candidate_id
        for item in borda.candidates
    ] == ["a", "b", "c"]

    assert [
        item.borda_score
        for item in borda.candidates
    ] == [2.0, 1.0, 0.0]


def test_component_values_are_normalized_before_ballot():
    source = result(
        candidate(
            "a",
            [
                component(
                    "clinical",
                    contribution=5,
                    maximum=10,
                )
            ],
        ),
        candidate(
            "b",
            [
                component(
                    "clinical",
                    contribution=6,
                    maximum=20,
                )
            ],
        ),
    )

    borda = aggregate_borda(source)

    # a = 0.50 normalized evidence
    # b = 0.30 normalized evidence
    assert (
        borda.candidates[0].candidate_id
        == "a"
    )


def test_ties_receive_average_rank_and_points():
    source = result(
        candidate(
            "a",
            [component("clinical", 10)],
        ),
        candidate(
            "b",
            [component("clinical", 10)],
        ),
        candidate(
            "c",
            [component("clinical", 0, status="NEUTRAL")],
        ),
    )

    borda = aggregate_borda(source)

    by_id = {
        item.candidate_id: item
        for item in borda.candidates
    }

    a = entry(
        by_id["a"],
        "clinical",
    )

    b = entry(
        by_id["b"],
        "clinical",
    )

    c = entry(
        by_id["c"],
        "clinical",
    )

    assert a.rank == pytest.approx(1.5)
    assert b.rank == pytest.approx(1.5)

    assert (
        a.borda_points
        == pytest.approx(1.5)
    )

    assert (
        b.borda_points
        == pytest.approx(1.5)
    )

    assert c.rank == pytest.approx(3.0)
    assert c.borda_points == 0.0


@pytest.mark.parametrize(
    "status",
    [
        "UNAVAILABLE",
        "BLOCKED",
        "ERROR",
    ],
)
def test_missing_or_blocked_evidence_abstains(
    status,
):
    source = result(
        candidate(
            "a",
            [
                component(
                    "clinical",
                    10,
                )
            ],
        ),
        candidate(
            "b",
            [
                component(
                    "clinical",
                    0,
                    status=status,
                )
            ],
        ),
    )

    borda = aggregate_borda(source)

    by_id = {
        item.candidate_id: item
        for item in borda.candidates
    }

    b = entry(
        by_id["b"],
        "clinical",
    )

    assert b.abstained is True
    assert b.rank is None
    assert b.normalized_evidence is None
    assert b.borda_points == 0.0


def test_zero_weight_component_is_not_a_ballot():
    source = result(
        candidate(
            "a",
            [
                component(
                    "phenotype",
                    contribution=0,
                    maximum=0,
                    status="UNAVAILABLE",
                ),
                component(
                    "clinical",
                    contribution=10,
                    maximum=10,
                ),
            ],
        )
    )

    borda = aggregate_borda(source)

    assert borda.component_names == (
        "clinical",
    )


def test_truth_has_no_place_in_borda_contract():
    source = result(
        candidate(
            "a",
            [component("clinical", 10)],
        )
    )

    borda = aggregate_borda(source)

    data = borda.to_dict()

    assert "truth" not in data

    assert (
        data["semantic_boundaries"][
            "benchmark_truth_used"
        ]
        is False
    )

    assert (
        data["semantic_boundaries"][
            "pathogenicity_probability"
        ]
        is False
    )

    assert (
        data["semantic_boundaries"][
            "diagnosis"
        ]
        is False
    )


def test_output_is_deterministic_for_complete_tie():
    source = result(
        candidate(
            "z",
            [component("clinical", 5)],
        ),
        candidate(
            "a",
            [component("clinical", 5)],
        ),
    )

    borda = aggregate_borda(source)

    assert [
        item.candidate_id
        for item in borda.candidates
    ] == ["a", "z"]


def test_mean_score_avoids_direct_missing_ballot_penalty():
    source = result(
        candidate(
            "complete",
            [
                component("clinical", 10),
                component("functional", 0, status="NEUTRAL"),
            ],
        ),
        candidate(
            "partial",
            [
                component("clinical", 10),
                component(
                    "functional",
                    0,
                    status="UNAVAILABLE",
                ),
            ],
        ),
    )

    borda = aggregate_borda(source)

    by_id = {
        item.candidate_id: item
        for item in borda.candidates
    }

    assert (
        by_id["partial"].ballots_participated
        == 1
    )

    assert (
        by_id["complete"].ballots_participated
        == 2
    )

    assert any(
        "BORDA_PARTIAL_EVIDENCE"
        in warning
        for warning
        in by_id["partial"].warnings
    )


def test_borda_score_is_not_weighted_prioritization_score():
    source = result(
        candidate(
            "a",
            [
                component("clinical", 10),
                component("functional", 5),
            ],
        ),
        candidate(
            "b",
            [
                component("clinical", 5),
                component("functional", 10),
            ],
        ),
    )

    borda = aggregate_borda(source)

    data = borda.to_dict()

    assert (
        data["strategy"]
        == "borda"
    )

    assert (
        data["engine"]["name"]
        == "borda_rank_aggregation"
    )

    assert all(
        "borda_score"
        in item
        for item in data["candidates"]
    )


def test_partial_one_hit_cannot_win_by_missingness_inflation():
    """Sparse evidence must not inflate aggregate Borda score.

    UNAVAILABLE evidence is not treated as negative biological
    evidence. Instead, evidence coverage adjusts confidence in
    the aggregate score.
    """

    component_names = (
        "clinical",
        "population",
        "functional",
        "gene_disease",
        "inheritance",
        "phenotype",
    )

    def evidence_component(
        name,
        value=None,
    ):
        if value is None:
            return ComponentContribution(
                name=name,
                contribution=0.0,
                max_contribution=10.0,
                status="UNAVAILABLE",
            )

        return ComponentContribution(
            name=name,
            contribution=value * 10.0,
            max_contribution=10.0,
            status=(
                "SUPPORTING"
                if value > 0
                else "NEUTRAL"
            ),
        )

    def evidence_candidate(
        candidate_id,
        values,
    ):
        return RankedCandidate(
            candidate_id=candidate_id,
            variant_key=(
                f"1:200:A:{candidate_id}"
            ),
            gene_symbol=(
                candidate_id.upper()
            ),
            rank=1,
            prioritization_score=0.0,
            components=tuple(
                evidence_component(
                    name,
                    values.get(name),
                )
                for name
                in component_names
            ),
        )

    complete = evidence_candidate(
        "complete-balanced",
        {
            name: 0.80
            for name in component_names
        },
    )

    partial = evidence_candidate(
        "partial-one-hit",
        {
            "phenotype": 1.00,
        },
    )

    weak = evidence_candidate(
        "weak-complete",
        {
            name: 0.20
            for name in component_names
        },
    )

    source = RankingResult(
        mode=(
            RankingMode.PHENOTYPE_SEMANTIC
        ),
        candidates=(
            complete,
            partial,
            weak,
        ),
        engine_version="synthetic",
    )

    borda = aggregate_borda(
        source
    )

    by_id = {
        item.candidate_id: item
        for item in borda.candidates
    }

    assert (
        borda.candidates[0].candidate_id
        == "complete-balanced"
    )

    assert (
        by_id["complete-balanced"]
        .ballots_participated
        == 6
    )

    assert (
        by_id["partial-one-hit"]
        .ballots_participated
        == 1
    )

    assert (
        by_id["complete-balanced"]
        .borda_score
        > by_id["partial-one-hit"]
        .borda_score
    )

    assert (
        by_id["partial-one-hit"]
        .borda_score
        == pytest.approx(
            2.0 / 6.0,
            abs=1e-6,
        )
    )


def test_equal_coverage_adjustment_preserves_relative_order():
    """Equal evidence coverage must preserve Borda ordering.

    Coverage adjustment may rescale scores, but when candidates
    participate in the same number of ballots it must not reverse
    their relative Borda order.
    """

    names = (
        "clinical",
        "population",
        "functional",
        "phenotype",
        "gene_disease",
        "inheritance",
    )

    def make_component(
        name,
        value,
    ):
        if value is None:
            return ComponentContribution(
                name=name,
                contribution=0.0,
                max_contribution=10.0,
                status="UNAVAILABLE",
            )

        return ComponentContribution(
            name=name,
            contribution=value * 10.0,
            max_contribution=10.0,
            status=(
                "SUPPORTING"
                if value > 0
                else "NEUTRAL"
            ),
        )

    def make_candidate(
        candidate_id,
        values,
    ):
        return RankedCandidate(
            candidate_id=candidate_id,
            variant_key=(
                f"1:300:A:{candidate_id}"
            ),
            gene_symbol=candidate_id.upper(),
            rank=1,
            prioritization_score=0.0,
            components=tuple(
                make_component(
                    name,
                    values.get(name),
                )
                for name in names
            ),
        )

    first = make_candidate(
        "first",
        {
            "clinical": 0.9,
            "population": 0.9,
            "functional": 0.9,
            "phenotype": 0.4,
        },
    )

    second = make_candidate(
        "second",
        {
            "clinical": 0.4,
            "population": 0.4,
            "functional": 0.4,
            "phenotype": 0.9,
        },
    )

    source = RankingResult(
        mode=(
            RankingMode.PHENOTYPE_SEMANTIC
        ),
        candidates=(
            first,
            second,
        ),
        engine_version="synthetic",
    )

    borda = aggregate_borda(
        source
    )

    by_id = {
        item.candidate_id: item
        for item in borda.candidates
    }

    assert (
        by_id["first"].ballots_participated
        == 4
    )

    assert (
        by_id["second"].ballots_participated
        == 4
    )

    assert (
        by_id["first"].borda_score
        > by_id["second"].borda_score
    )

    assert (
        borda.candidates[0].candidate_id
        == "first"
    )


def test_equal_final_borda_scores_share_competition_rank():
    source = result(
        candidate(
            "a",
            [
                component("clinical", 10),
                component("functional", 0, status="NEUTRAL"),
            ],
        ),
        candidate(
            "b",
            [
                component("clinical", 0, status="NEUTRAL"),
                component("functional", 10),
            ],
        ),
        candidate(
            "c",
            [
                component("clinical", 0, status="NEUTRAL"),
                component("functional", 0, status="NEUTRAL"),
            ],
        ),
    )

    borda = aggregate_borda(source)

    by_id = {
        item.candidate_id: item
        for item in borda.candidates
    }

    assert (
        by_id["a"].borda_score
        == by_id["b"].borda_score
    )

    assert by_id["a"].rank == 1
    assert by_id["b"].rank == 1

    assert by_id["c"].rank == 3


def test_candidate_id_does_not_break_final_rank_tie():
    source = result(
        candidate(
            "z-candidate",
            [component("clinical", 5)],
        ),
        candidate(
            "a-candidate",
            [component("clinical", 5)],
        ),
    )

    borda = aggregate_borda(source)

    # Presentation remains deterministic.
    assert [
        item.candidate_id
        for item in borda.candidates
    ] == [
        "a-candidate",
        "z-candidate",
    ]

    # But candidate ID must not manufacture
    # a ranking distinction.
    assert [
        item.rank
        for item in borda.candidates
    ] == [1, 1]
