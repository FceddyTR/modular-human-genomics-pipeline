from dataclasses import dataclass, fields

import pytest

from genomics_platform.ranking.benchmark.real_cases import (
    CaseProvenance,
    RealCaseInput,
    validate_case_truth_pair,
    validate_real_case_input,
)
from genomics_platform.ranking.benchmark.real_cases.validation import (
    FORBIDDEN_RANKING_FIELDS,
)
from genomics_platform.ranking.benchmark.truth_set import (
    RankingTruth,
)


@dataclass(frozen=True)
class CandidateStub:
    candidate_id: str


def _provenance():
    return CaseProvenance(
        source_name="Public benchmark source",
        source_reference="public-case-reference",
        source_type="published_case",
        accessed_date="2026-09-16",
    )


def _case():
    return RealCaseInput(
        case_id="public-case-001",
        assembly="GRCh38",
        present_hpo_terms=(
            "HP:0001250",
            "HP:0001250",
        ),
        absent_hpo_terms=(
            "HP:0001263",
        ),
        candidates=(
            CandidateStub("candidate-a"),
            CandidateStub("candidate-b"),
        ),
        provenance=_provenance(),
    )


def test_real_case_input_is_truth_blind():
    case = _case()

    field_names = {
        field.name
        for field in fields(case)
    }

    assert (
        field_names
        & FORBIDDEN_RANKING_FIELDS
    ) == set()

    assert not hasattr(
        case,
        "causal_candidate_ids",
    )

    assert not hasattr(
        case,
        "truth",
    )

    validate_real_case_input(case)


def test_real_case_normalizes_hpo_terms():
    case = _case()

    assert case.present_hpo_terms == (
        "HP:0001250",
    )

    assert case.absent_hpo_terms == (
        "HP:0001263",
    )


def test_real_case_rejects_present_absent_overlap():
    with pytest.raises(
        ValueError,
        match="both present and absent",
    ):
        RealCaseInput(
            case_id="public-case-overlap",
            assembly="GRCh38",
            present_hpo_terms=(
                "HP:0001250",
            ),
            absent_hpo_terms=(
                "HP:0001250",
            ),
            candidates=(
                CandidateStub("candidate-a"),
            ),
            provenance=_provenance(),
        )


def test_real_case_rejects_duplicate_candidate_ids():
    with pytest.raises(
        ValueError,
        match="unique",
    ):
        RealCaseInput(
            case_id="public-case-duplicate",
            assembly="GRCh38",
            present_hpo_terms=(),
            absent_hpo_terms=(),
            candidates=(
                CandidateStub("duplicate"),
                CandidateStub("duplicate"),
            ),
            provenance=_provenance(),
        )


def test_truth_pair_is_validated_only_at_benchmark_boundary():
    case = _case()

    truth = RankingTruth(
        case_id=case.case_id,
        causal_candidate_ids=(
            "candidate-b",
        ),
    )

    validate_case_truth_pair(
        case,
        truth,
    )


def test_truth_pair_rejects_case_id_mismatch():
    case = _case()

    truth = RankingTruth(
        case_id="different-case",
        causal_candidate_ids=(
            "candidate-b",
        ),
    )

    with pytest.raises(
        ValueError,
        match="does not match",
    ):
        validate_case_truth_pair(
            case,
            truth,
        )


def test_truth_pair_rejects_unknown_causal_candidate():
    case = _case()

    truth = RankingTruth(
        case_id=case.case_id,
        causal_candidate_ids=(
            "not-in-cohort",
        ),
    )

    with pytest.raises(
        ValueError,
        match="absent from the ranked cohort",
    ):
        validate_case_truth_pair(
            case,
            truth,
        )


def test_serialized_real_case_contains_no_truth_labels():
    payload = _case().to_dict()

    serialized = repr(
        payload
    ).lower()

    assert "causal_candidate_ids" not in serialized
    assert "is_causal" not in serialized
    assert "expected_rank" not in serialized

    assert (
        payload["semantic_boundaries"][
            "truth_present_in_ranking_input"
        ]
        is False
    )

    assert (
        payload["semantic_boundaries"][
            "truth_used_for_ranking"
        ]
        is False
    )
