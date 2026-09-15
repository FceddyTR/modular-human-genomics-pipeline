import pytest

from genomics_platform.interpretation.phenotype_interpreter import (
    interpret_phenotype,
)


def make_bundle():
    return {
        "schema_version": "1.0",
        "query": {
            "gene_symbol": "NCSTN",
        },
        "records": [
            {
                "evidence": {
                    "source": "HPO",
                    "evidence_type": "gene_phenotype",
                    "gene_symbol": "NCSTN",
                    "disease_id": "OMIM:142690",
                    "phenotype_ids": [
                        "HP:0000987"
                    ],
                    "source_record_id": (
                        "NCSTN|HP:0000987"
                    ),
                    "evidence_description": (
                        "Atypical scarring of skin"
                    ),
                },
                "disease_normalization": {
                    "status": "FOUND",
                    "mondo_id": "MONDO:0007728",
                    "mondo_label": (
                        "acne inversa, familial, 1"
                    ),
                },
            },
            {
                "evidence": {
                    "source": "HPO",
                    "evidence_type": "gene_phenotype",
                    "gene_symbol": "NCSTN",
                    "disease_id": "OMIM:142690",
                    "phenotype_ids": [
                        "HP:0040154"
                    ],
                    "source_record_id": (
                        "NCSTN|HP:0040154"
                    ),
                    "evidence_description": (
                        "Acne inversa"
                    ),
                },
                "disease_normalization": {
                    "status": "FOUND",
                    "mondo_id": "MONDO:0007728",
                    "mondo_label": (
                        "acne inversa, familial, 1"
                    ),
                },
            },
            {
                "evidence": {
                    "source": "GenCC",
                    "gene_symbol": "NCSTN",
                },
                "disease_normalization": {
                    "status": "FOUND",
                },
            },
        ],
    }


def get_observation(dimension, code):
    return next(
        observation
        for observation in dimension.observations
        if observation.code == code
    )


def test_none_bundle_is_unavailable():
    dimension = interpret_phenotype(
        None,
        present_hpo_terms=["HP:0040154"],
    )

    assert dimension.status == "UNAVAILABLE"


def test_no_patient_terms_is_unavailable():
    dimension = interpret_phenotype(
        make_bundle()
    )

    assert dimension.status == "UNAVAILABLE"


def test_exact_present_match_is_available():
    dimension = interpret_phenotype(
        make_bundle(),
        present_hpo_terms=[
            "HP:0040154",
            "HP:0000987",
        ],
    )

    assert dimension.status == "AVAILABLE"

    summary = get_observation(
        dimension,
        "PHENOTYPE_EXACT_MATCH_SUMMARY",
    )

    assert summary.value["matched_present"] == [
        "HP:0040154",
        "HP:0000987",
    ]

    assert summary.value["unmatched_present"] == []


def test_unmatched_present_term_is_partial():
    dimension = interpret_phenotype(
        make_bundle(),
        present_hpo_terms=[
            "HP:0040154",
            "HP:1234567",
        ],
    )

    assert dimension.status == "PARTIAL"

    summary = get_observation(
        dimension,
        "PHENOTYPE_EXACT_MATCH_SUMMARY",
    )

    assert summary.value["matched_present"] == [
        "HP:0040154"
    ]

    assert summary.value["unmatched_present"] == [
        "HP:1234567"
    ]


def test_matched_negative_term_is_preserved():
    dimension = interpret_phenotype(
        make_bundle(),
        present_hpo_terms=[
            "HP:0040154"
        ],
        absent_hpo_terms=[
            "HP:0000987"
        ],
    )

    assert dimension.status == "PARTIAL"

    summary = get_observation(
        dimension,
        "PHENOTYPE_EXACT_MATCH_SUMMARY",
    )

    assert summary.value["matched_absent"] == [
        "HP:0000987"
    ]


def test_unmatched_negative_does_not_force_partial():
    dimension = interpret_phenotype(
        make_bundle(),
        present_hpo_terms=[
            "HP:0040154"
        ],
        absent_hpo_terms=[
            "HP:1234567"
        ],
    )

    assert dimension.status == "AVAILABLE"


def test_duplicate_terms_are_deduplicated():
    dimension = interpret_phenotype(
        make_bundle(),
        present_hpo_terms=[
            "HP:0040154",
            "hp:0040154",
            " HP:0040154 ",
        ],
    )

    context = get_observation(
        dimension,
        "PHENOTYPE_PATIENT_CONTEXT",
    )

    assert context.value[
        "present_hpo_terms"
    ] == ["HP:0040154"]


def test_same_term_present_and_absent_rejected():
    with pytest.raises(ValueError):
        interpret_phenotype(
            make_bundle(),
            present_hpo_terms=[
                "HP:0040154"
            ],
            absent_hpo_terms=[
                "HP:0040154"
            ],
        )


def test_invalid_hpo_identifier_rejected():
    with pytest.raises(ValueError):
        interpret_phenotype(
            make_bundle(),
            present_hpo_terms=[
                "NOT_AN_HPO_TERM"
            ],
        )


def test_wrong_bundle_type_rejected():
    with pytest.raises(TypeError):
        interpret_phenotype(
            ["wrong"],
            present_hpo_terms=[
                "HP:0040154"
            ],
        )


def test_only_hpo_source_records_are_used():
    dimension = interpret_phenotype(
        make_bundle(),
        present_hpo_terms=[
            "HP:0040154"
        ],
    )

    source_records = [
        observation
        for observation in dimension.observations
        if observation.code
        == "PHENOTYPE_SOURCE_ASSOCIATION"
    ]

    assert len(source_records) == 2

    assert all(
        observation.source == "HPO"
        for observation in source_records
    )


def test_no_ranking_or_classification_output():
    dimension = interpret_phenotype(
        make_bundle(),
        present_hpo_terms=[
            "HP:0040154"
        ],
    )

    data = dimension.to_dict()

    forbidden = {
        "score",
        "ranking_score",
        "tier",
        "acmg_classification",
        "pathogenicity_probability",
        "diagnosis",
    }

    assert forbidden.isdisjoint(data.keys())


def test_source_record_preserves_normalization():
    dimension = interpret_phenotype(
        make_bundle(),
        present_hpo_terms=[
            "HP:0040154"
        ],
    )

    source_records = [
        observation
        for observation in dimension.observations
        if observation.code
        == "PHENOTYPE_SOURCE_ASSOCIATION"
    ]

    assert any(
        observation.details[
            "disease_normalization"
        ]["mondo_id"]
        == "MONDO:0007728"
        for observation in source_records
    )
