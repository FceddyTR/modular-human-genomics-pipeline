from __future__ import annotations

import pytest

from genomics_platform.phenotype.hpo_semantic_similarity import (
    HPOSemanticSimilarity,
)
from genomics_platform.ranking.benchmark.controlled_cases import (
    build_initial_registry,
)
from genomics_platform.ranking.benchmark.controlled_suite import (
    validate_suite,
)
from genomics_platform.ranking.benchmark.synthetic_factory import (
    make_synthetic_candidate,
)
from genomics_platform.ranking.ranking_engine import (
    rank_candidates,
)
from genomics_platform.ranking.ranking_modes import (
    RankingMode,
)


@pytest.fixture(scope="module")
def semantic():
    return HPOSemanticSimilarity(
        ontology_database=(
            "data/evidence/hpo_ontology/"
            "hpo_ontology.sqlite"
        ),
        information_content_database=(
            "data/evidence/hpo_ontology/"
            "hpo_ic.sqlite"
        ),
    )


def test_initial_registry_contains_12_cases():
    registry = build_initial_registry()

    assert len(
        registry.case_ids()
    ) == 12

    assert registry.family_counts() == {
        "functional": 4,
        "phenotype": 4,
        "population": 4,
    }


def test_initial_suite_is_structurally_valid():
    registry = build_initial_registry()

    cases = registry.build_all()

    validate_suite(cases)


def test_candidates_do_not_contain_truth_labels():
    registry = build_initial_registry()

    for case in registry.build_all():
        for candidate in case.ranking_candidates():
            assert not hasattr(
                candidate,
                "is_causal",
            )


def test_all_cases_run_weighted_exact():
    registry = build_initial_registry()

    for case in registry.build_all():
        result = rank_candidates(
            case.ranking_candidates(),
            mode=RankingMode.PHENOTYPE,
        )

        assert len(
            result.candidates
        ) == len(
            case.candidates
        )


def test_all_cases_run_variant_first():
    registry = build_initial_registry()

    for case in registry.build_all():
        result = rank_candidates(
            case.ranking_candidates(),
            mode=RankingMode.VARIANT_FIRST,
        )

        assert len(
            result.candidates
        ) == len(
            case.candidates
        )


def test_all_cases_run_semantic(
    semantic,
):
    registry = build_initial_registry()

    for case in registry.build_all():
        result = rank_candidates(
            case.ranking_candidates(),
            mode=(
                RankingMode.PHENOTYPE_SEMANTIC
            ),
            semantic_engine=semantic,
        )

        assert len(
            result.candidates
        ) == len(
            case.candidates
        )


def test_population_not_found_is_not_af_zero():
    registry = build_initial_registry()

    case = registry.build(
        "population-not-found-001"
    )

    result = rank_candidates(
        case.ranking_candidates(),
        mode=RankingMode.VARIANT_FIRST,
    )

    by_id = {
        candidate.candidate_id:
        candidate
        for candidate in result.candidates
    }

    decoy = by_id[
        "not-found-decoy"
    ]

    population = next(
        component
        for component in decoy.components
        if component.name == "population"
    )

    assert population.contribution == 0.0
    assert population.status == "UNAVAILABLE"


def test_synthetic_factory_preserves_gene_consistency():
    candidate = make_synthetic_candidate(
        candidate_id="factory-check",
        gene="TESTGENE",
        position=79999,
    )

    assert (
        candidate.gene_consistency
        == "MATCH"
    )


def test_truth_ids_are_resolved_only_from_case():
    registry = build_initial_registry()

    case = registry.build(
        "population-rare-001"
    )

    result = rank_candidates(
        case.ranking_candidates(),
        mode=RankingMode.VARIANT_FIRST,
    )

    ranked_ids = {
        candidate.candidate_id
        for candidate in result.candidates
    }

    assert set(
        case.causal_candidate_ids()
    ) <= ranked_ids
