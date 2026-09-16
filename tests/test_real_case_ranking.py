import inspect

from genomics_platform.phenotype.hpo_semantic_similarity import (
    HPOSemanticSimilarity,
)
from genomics_platform.ranking.benchmark.truth_set import (
    RankingTruth,
)
from genomics_platform.ranking.benchmark.real_cases.ranking import (
    RealCaseRankingBundle,
    rank_real_case,
)
from genomics_platform.ranking.benchmark.real_cases.validation import (
    validate_case_truth_pair,
)
from genomics_platform.ranking.benchmark.real_cases.contract import (
    CaseProvenance,
    RealCaseInput,
)
from test_case_model import make_candidate


def _semantic():
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


def _case():
    candidate = make_candidate()

    return RealCaseInput(
        case_id="public-case-001",
        assembly="GRCh38",
        present_hpo_terms=(
            candidate.case_context
            .present_hpo_terms
        ),
        absent_hpo_terms=(
            candidate.case_context
            .absent_hpo_terms
        ),
        candidates=(candidate,),
        provenance=CaseProvenance(
            source_name="Public source",
            source_reference="reference-001",
            source_type="published_case",
            accessed_date="2026-09-16",
        ),
    )


def test_rank_real_case_has_no_truth_parameter():
    signature = inspect.signature(
        rank_real_case
    )

    assert "truth" not in signature.parameters
    assert "causal_candidate_ids" not in (
        signature.parameters
    )

    assert tuple(
        signature.parameters
    ) == (
        "case",
        "semantic_engine",
    )


def test_real_case_runs_five_expected_strategies():
    bundle = rank_real_case(
        _case(),
        semantic_engine=_semantic(),
    )

    assert isinstance(
        bundle,
        RealCaseRankingBundle,
    )

    assert tuple(
        run.strategy
        for run in bundle.runs
    ) == (
        "weighted_variant_first",
        "weighted_exact",
        "weighted_semantic",
        "borda_semantic",
        "rank_geomean_semantic",
    )


def test_real_case_ranking_is_deterministic():
    case = _case()
    semantic = _semantic()

    first = rank_real_case(
        case,
        semantic_engine=semantic,
    )

    second = rank_real_case(
        case,
        semantic_engine=semantic,
    )

    assert first == second
    assert first.to_dict() == second.to_dict()


def test_ranking_bundle_declares_truth_absent():
    bundle = rank_real_case(
        _case(),
        semantic_engine=_semantic(),
    )

    serialized = bundle.to_dict()

    boundaries = serialized[
        "semantic_boundaries"
    ]

    assert (
        boundaries["truth_used_for_ranking"]
        is False
    )

    assert (
        boundaries["benchmark_truth_present"]
        is False
    )


def test_truth_is_introduced_only_after_ranking():
    case = _case()

    bundle = rank_real_case(
        case,
        semantic_engine=_semantic(),
    )

    frozen = bundle.to_dict()

    truth = RankingTruth(
        case_id=case.case_id,
        causal_candidate_ids=(
            case.candidates[0].candidate_id,
        ),
    )

    validate_case_truth_pair(
        case,
        truth,
    )

    assert bundle.to_dict() == frozen


def test_truth_object_cannot_change_frozen_ranking():
    case = _case()

    bundle = rank_real_case(
        case,
        semantic_engine=_semantic(),
    )

    frozen = bundle.to_dict()

    truth_a = RankingTruth(
        case_id=case.case_id,
        causal_candidate_ids=(
            case.candidates[0].candidate_id,
        ),
    )

    validate_case_truth_pair(
        case,
        truth_a,
    )

    assert bundle.to_dict() == frozen

    # A different truth object is benchmark metadata only.
    # It is deliberately not accepted by rank_real_case.
    truth_b = RankingTruth(
        case_id=case.case_id,
        causal_candidate_ids=(
            case.candidates[0].candidate_id,
        ),
    )

    validate_case_truth_pair(
        case,
        truth_b,
    )

    assert bundle.to_dict() == frozen
