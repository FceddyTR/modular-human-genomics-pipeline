from dataclasses import FrozenInstanceError

import pytest

from genomics_platform.interpretation.case_model import (
    CaseContext,
    CandidateCase,
    build_candidate_case,
)
from genomics_platform.interpretation.interpretation_engine import (
    interpret_variant,
)

from test_interpretation_engine import make_contract


def make_candidate():
    contract = make_contract()
    profile = interpret_variant(contract)

    return build_candidate_case(
        candidate_id="candidate-001",
        evidence_contract=contract,
        interpretation_profile=profile,
        gene_symbol="NCSTN",
        present_hpo_terms=[
            "HP:0040154",
            "HP:0000987",
        ],
        genotype_state="heterozygous",
        proband_sex="female",
    )


def test_build_candidate_case():
    candidate = make_candidate()

    assert (
        candidate.candidate_id
        == "candidate-001"
    )
    assert candidate.gene_symbol == "NCSTN"


def test_variant_key_matches_contract():
    candidate = make_candidate()

    assert (
        candidate.variant_key
        == candidate.evidence_contract
        .identity.variant_key
    )


def test_case_context_normalizes_hpo():
    context = CaseContext(
        present_hpo_terms=(
            "hp:0040154",
            "HP:0040154",
            " HP:0000987 ",
        )
    )

    assert context.present_hpo_terms == (
        "HP:0040154",
        "HP:0000987",
    )


def test_present_absent_conflict_rejected():
    with pytest.raises(ValueError):
        CaseContext(
            present_hpo_terms=(
                "HP:0040154",
            ),
            absent_hpo_terms=(
                "HP:0040154",
            ),
        )


def test_invalid_hpo_rejected():
    with pytest.raises(ValueError):
        CaseContext(
            present_hpo_terms=(
                "BAD",
            )
        )


def test_invalid_genotype_rejected():
    with pytest.raises(ValueError):
        CaseContext(
            genotype_state="bad-state"
        )


def test_invalid_sex_rejected():
    with pytest.raises(ValueError):
        CaseContext(
            proband_sex="invalid"
        )


def test_candidate_is_frozen():
    candidate = make_candidate()

    with pytest.raises(FrozenInstanceError):
        candidate.gene_symbol = "OTHER"


def test_case_context_is_frozen():
    context = CaseContext()

    with pytest.raises(FrozenInstanceError):
        context.proband_sex = "male"


def test_variant_mismatch_rejected():
    contract = make_contract()
    profile = interpret_variant(contract)

    object.__setattr__(
        profile,
        "variant_key",
        "GRCh38:1:1:A:T",
    )

    with pytest.raises(ValueError):
        CandidateCase(
            candidate_id="bad",
            evidence_contract=contract,
            interpretation_profile=profile,
        )


def test_serialization_contains_contract():
    candidate = make_candidate()

    data = candidate.to_dict()

    assert "evidence_contract" in data
    assert "interpretation_profile" in data
    assert "case_context" in data


def test_serialization_preserves_context():
    data = make_candidate().to_dict()

    assert data["case_context"][
        "genotype_state"
    ] == "heterozygous"

    assert data["case_context"][
        "present_hpo_terms"
    ] == [
        "HP:0040154",
        "HP:0000987",
    ]


def test_available_dimensions_forwarded():
    candidate = make_candidate()

    assert set(
        candidate.available_dimensions
    ) == {
        "population",
        "clinical",
        "functional",
    }


def test_no_error_dimensions():
    candidate = make_candidate()

    assert candidate.error_dimensions == ()


def test_no_ranking_or_classification():
    data = make_candidate().to_dict()

    boundaries = data[
        "semantic_boundaries"
    ]

    assert (
        boundaries["candidate_ranking"]
        is False
    )
    assert (
        boundaries["candidate_tier"]
        is False
    )
    assert (
        boundaries["variant_classification"]
        is False
    )
    assert (
        boundaries[
            "acmg_amp_classification"
        ]
        is False
    )
    assert (
        boundaries[
            "pathogenicity_probability"
        ]
        is False
    )
    assert boundaries["diagnosis"] is False


def test_candidate_id_required():
    contract = make_contract()
    profile = interpret_variant(contract)

    with pytest.raises(ValueError):
        build_candidate_case(
            candidate_id=" ",
            evidence_contract=contract,
            interpretation_profile=profile,
        )


def test_wrong_contract_type_rejected():
    contract = make_contract()
    profile = interpret_variant(contract)

    with pytest.raises(TypeError):
        CandidateCase(
            candidate_id="candidate",
            evidence_contract={},
            interpretation_profile=profile,
        )


def test_wrong_profile_type_rejected():
    contract = make_contract()

    with pytest.raises(TypeError):
        CandidateCase(
            candidate_id="candidate",
            evidence_contract=contract,
            interpretation_profile={},
        )
