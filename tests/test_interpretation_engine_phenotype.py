import pytest

from genomics_platform.interpretation.interpretation_engine import (
    interpret_variant,
)

from test_interpretation_engine import (
    make_contract,
)


def make_gene_context():
    return {
        "schema_version": "1.0",
        "engine": {
            "name": "gene_evidence_bundle",
            "version": "0.2.0",
        },
        "query": {
            "gene_symbol": "NCSTN",
        },
        "scope": {
            "clinical_classification": False,
            "variant_pathogenicity": False,
            "candidate_ranking": False,
        },
        "summary": {
            "record_count": 2,
            "source_counts": {
                "HPO": 2,
            },
            "disease_normalization": {
                "FOUND": 2,
            },
            "normalized_disease_count": 1,
        },
        "normalized_diseases": [
            {
                "mondo_id": "MONDO:0007728",
                "mondo_label": (
                    "acne inversa, familial, 1"
                ),
            }
        ],
        "records": [
            {
                "evidence": {
                    "schema_version": "1.1",
                    "source": "HPO",
                    "evidence_type": (
                        "gene_phenotype"
                    ),
                    "gene_symbol": "NCSTN",
                    "gene_id": (
                        "NCBIGene:23385"
                    ),
                    "disease_id": (
                        "OMIM:142690"
                    ),
                    "phenotype_ids": [
                        "HP:0040154"
                    ],
                    "source_record_id": (
                        "NCBIGene:23385|"
                        "OMIM:142690|"
                        "HP:0040154"
                    ),
                    "evidence_description": (
                        "HP:0040154 "
                        "Acne inversa"
                    ),
                },
                "disease_normalization": {
                    "status": "FOUND",
                    "input_id": (
                        "OMIM:142690"
                    ),
                    "matched_source_id": (
                        "OMIM:142690"
                    ),
                    "mondo_id": (
                        "MONDO:0007728"
                    ),
                    "mondo_label": (
                        "acne inversa, "
                        "familial, 1"
                    ),
                },
            },
            {
                "evidence": {
                    "schema_version": "1.1",
                    "source": "HPO",
                    "evidence_type": (
                        "gene_phenotype"
                    ),
                    "gene_symbol": "NCSTN",
                    "gene_id": (
                        "NCBIGene:23385"
                    ),
                    "disease_id": (
                        "OMIM:142690"
                    ),
                    "phenotype_ids": [
                        "HP:0000987"
                    ],
                    "source_record_id": (
                        "NCBIGene:23385|"
                        "OMIM:142690|"
                        "HP:0000987"
                    ),
                    "evidence_description": (
                        "HP:0000987 "
                        "Atypical scarring "
                        "of skin"
                    ),
                },
                "disease_normalization": {
                    "status": "FOUND",
                    "input_id": (
                        "OMIM:142690"
                    ),
                    "matched_source_id": (
                        "OMIM:142690"
                    ),
                    "mondo_id": (
                        "MONDO:0007728"
                    ),
                    "mondo_label": (
                        "acne inversa, "
                        "familial, 1"
                    ),
                },
            },
        ],
    }


def test_old_behavior_has_no_phenotype_dimension():
    profile = interpret_variant(
        make_contract()
    )

    assert profile.phenotype is None


def test_gene_context_alone_does_not_request_phenotype():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=(
            make_gene_context()
        ),
    )

    assert profile.gene_disease is not None
    assert profile.phenotype is None


def test_present_hpo_runs_phenotype_interpreter():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=(
            make_gene_context()
        ),
        present_hpo_terms=[
            "HP:0040154"
        ],
    )

    assert profile.phenotype is not None
    assert (
        profile.phenotype.status
        == "AVAILABLE"
    )


def test_exact_match_is_preserved():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=(
            make_gene_context()
        ),
        present_hpo_terms=[
            "HP:0040154"
        ],
    )

    summary = next(
        observation
        for observation
        in profile.phenotype.observations
        if observation.code
        == "PHENOTYPE_EXACT_MATCH_SUMMARY"
    )

    assert summary.value[
        "matched_present"
    ] == ["HP:0040154"]


def test_unmatched_present_is_partial():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=(
            make_gene_context()
        ),
        present_hpo_terms=[
            "HP:1234567"
        ],
    )

    assert (
        profile.phenotype.status
        == "PARTIAL"
    )


def test_negative_match_is_partial():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=(
            make_gene_context()
        ),
        absent_hpo_terms=[
            "HP:0000987"
        ],
    )

    assert (
        profile.phenotype.status
        == "PARTIAL"
    )


def test_phenotype_requested_without_gene_context():
    profile = interpret_variant(
        make_contract(),
        present_hpo_terms=[
            "HP:0040154"
        ],
    )

    assert profile.phenotype is not None
    assert (
        profile.phenotype.status
        == "UNAVAILABLE"
    )


def test_phenotype_serializes_in_profile():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=(
            make_gene_context()
        ),
        present_hpo_terms=[
            "HP:0040154"
        ],
    )

    data = profile.to_dict()

    assert (
        data["dimensions"]["phenotype"]
        is not None
    )

    assert (
        data["dimensions"]["phenotype"][
            "status"
        ]
        == "AVAILABLE"
    )


def test_phenotype_appears_available_dimensions():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=(
            make_gene_context()
        ),
        present_hpo_terms=[
            "HP:0040154"
        ],
    )

    assert (
        "phenotype"
        in profile.available_dimensions
    )


def test_no_ranking_or_classification_added():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=(
            make_gene_context()
        ),
        present_hpo_terms=[
            "HP:0040154"
        ],
    )

    data = profile.to_dict()

    forbidden = {
        "score",
        "ranking_score",
        "tier",
        "acmg_classification",
        "pathogenicity_probability",
        "diagnosis",
    }

    assert forbidden.isdisjoint(
        data.keys()
    )

    assert forbidden.isdisjoint(
        data["dimensions"][
            "phenotype"
        ].keys()
    )


def test_present_and_absent_conflict_propagates():
    with pytest.raises(ValueError):
        interpret_variant(
            make_contract(),
            gene_evidence_context=(
                make_gene_context()
            ),
            present_hpo_terms=[
                "HP:0040154"
            ],
            absent_hpo_terms=[
                "HP:0040154"
            ],
        )
