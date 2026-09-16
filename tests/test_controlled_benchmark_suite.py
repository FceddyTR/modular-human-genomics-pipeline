from __future__ import annotations

from dataclasses import dataclass

import pytest

from genomics_platform.ranking.benchmark.controlled_suite import (
    ControlledCase,
    ControlledCaseSpec,
    ControlledSuiteRegistry,
    validate_suite,
)


@dataclass(frozen=True)
class DummyCandidate:
    candidate_id: str


def make_case(
    case_id: str = "case-001",
    family: str = "phenotype",
) -> ControlledCase:
    spec = ControlledCaseSpec(
        case_id=case_id,
        family=family,
        description="Synthetic controlled case.",
        causal_candidate_ids=("causal",),
        expected_behavior=(
            "Causal candidate should remain "
            "independently benchmarkable."
        ),
        tags=("synthetic",),
    )

    return ControlledCase(
        spec=spec,
        candidates=(
            DummyCandidate("decoy"),
            DummyCandidate("causal"),
        ),
    )


def test_case_keeps_truth_separate_from_candidates():
    case = make_case()

    candidates = case.ranking_candidates()

    assert len(candidates) == 2
    assert not hasattr(
        candidates[0],
        "is_causal",
    )
    assert not hasattr(
        candidates[1],
        "is_causal",
    )

    assert (
        case.causal_candidate_ids()
        == ("causal",)
    )


def test_case_rejects_missing_causal_candidate():
    spec = ControlledCaseSpec(
        case_id="bad-case",
        family="phenotype",
        description="Bad case.",
        causal_candidate_ids=("missing",),
        expected_behavior="Validation failure.",
    )

    with pytest.raises(
        ValueError,
        match="absent from cohort",
    ):
        ControlledCase(
            spec=spec,
            candidates=(
                DummyCandidate("decoy"),
            ),
        )


def test_case_rejects_duplicate_candidate_ids():
    spec = ControlledCaseSpec(
        case_id="duplicate",
        family="population",
        description="Duplicate candidate case.",
        causal_candidate_ids=("causal",),
        expected_behavior="Validation failure.",
    )

    with pytest.raises(
        ValueError,
        match="unique",
    ):
        ControlledCase(
            spec=spec,
            candidates=(
                DummyCandidate("causal"),
                DummyCandidate("causal"),
            ),
        )


def test_registry_builds_cases():
    registry = ControlledSuiteRegistry()

    registry.register(
        "case-001",
        lambda: make_case("case-001"),
    )

    registry.register(
        "case-002",
        lambda: make_case(
            "case-002",
            family="population",
        ),
    )

    assert registry.case_ids() == (
        "case-001",
        "case-002",
    )

    cases = registry.build_all()

    assert len(cases) == 2

    assert registry.family_counts() == {
        "phenotype": 1,
        "population": 1,
    }


def test_registry_rejects_duplicate_case():
    registry = ControlledSuiteRegistry()

    registry.register(
        "case-001",
        lambda: make_case(),
    )

    with pytest.raises(
        ValueError,
        match="Duplicate",
    ):
        registry.register(
            "case-001",
            lambda: make_case(),
        )


def test_registry_detects_builder_id_mismatch():
    registry = ControlledSuiteRegistry()

    registry.register(
        "expected-id",
        lambda: make_case(
            "different-id"
        ),
    )

    with pytest.raises(
        ValueError,
        match="does not match",
    ):
        registry.build(
            "expected-id"
        )


def test_suite_validation():
    cases = (
        make_case(
            "phenotype-001",
            "phenotype",
        ),
        make_case(
            "population-001",
            "population",
        ),
    )

    validate_suite(cases)


def test_spec_declares_nonclinical_boundaries():
    data = make_case().spec.to_dict()

    boundaries = data[
        "semantic_boundaries"
    ]

    assert (
        boundaries[
            "synthetic_controlled_case"
        ]
        is True
    )

    assert (
        boundaries[
            "clinical_validation"
        ]
        is False
    )

    assert (
        boundaries[
            "truth_used_for_ranking"
        ]
        is False
    )

    assert (
        boundaries[
            "truth_used_for_benchmark"
        ]
        is True
    )
