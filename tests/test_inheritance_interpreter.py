import pytest

from genomics_platform.interpretation.inheritance_interpreter import (
    interpret_inheritance,
)


def make_bundle():
    return {
        "query": {
            "gene_symbol": "NCSTN",
        },
        "records": [
            {
                "evidence": {
                    "source": "GenCC",
                    "gene_symbol": "NCSTN",
                    "disease_id": "MONDO:0007728",
                    "inheritance": [
                        "Autosomal dominant"
                    ],
                    "phenotype_ids": [],
                    "source_record_id": "GENCC-1",
                    "evidence_description": (
                        "Gene-disease validity"
                    ),
                },
                "disease_normalization": {
                    "status": "FOUND",
                    "mondo_id": "MONDO:0007728",
                },
            },
            {
                "evidence": {
                    "source": "PanelApp",
                    "gene_symbol": "NCSTN",
                    "disease_id": None,
                    "inheritance": [
                        "MONOALLELIC, autosomal or "
                        "pseudoautosomal, imprinted "
                        "status unknown"
                    ],
                    "phenotype_ids": [],
                    "source_record_id": "PANEL-1",
                },
                "disease_normalization": {
                    "status": "NO_DISEASE_ID",
                },
            },
            {
                "evidence": {
                    "source": "HPO",
                    "gene_symbol": "NCSTN",
                    "disease_id": "OMIM:142690",
                    "inheritance": [],
                    "phenotype_ids": [
                        "HP:0000006"
                    ],
                    "source_record_id": "HPO-1",
                    "evidence_description": (
                        "HP:0000006 Autosomal "
                        "dominant inheritance"
                    ),
                },
                "disease_normalization": {
                    "status": "FOUND",
                    "mondo_id": "MONDO:0007728",
                },
            },
        ],
    }


def source_observations(dimension):
    return [
        observation
        for observation in dimension.observations
        if observation.code
        == "INHERITANCE_SOURCE_EVIDENCE"
    ]


def summary_observation(dimension):
    return next(
        observation
        for observation in dimension.observations
        if observation.code
        == "INHERITANCE_SUMMARY"
    )


def test_none_context_unavailable():
    result = interpret_inheritance(
        None,
        genotype_state="heterozygous",
    )

    assert result.status == "UNAVAILABLE"


def test_wrong_context_type_rejected():
    with pytest.raises(TypeError):
        interpret_inheritance(
            [],
            genotype_state="heterozygous",
        )


def test_invalid_genotype_rejected():
    with pytest.raises(ValueError):
        interpret_inheritance(
            make_bundle(),
            genotype_state="bad-state",
        )


def test_gencc_ad_normalized():
    result = interpret_inheritance(
        make_bundle(),
        genotype_state="heterozygous",
    )

    records = source_observations(result)

    gencc = next(
        record
        for record in records
        if record.source == "GenCC"
    )

    assert (
        gencc.value[
            "normalized_inheritance"
        ]
        == "AUTOSOMAL_DOMINANT"
    )


def test_panelapp_monoallelic_preserved_conservatively():
    result = interpret_inheritance(
        make_bundle(),
        genotype_state="heterozygous",
    )

    panel = next(
        record
        for record in source_observations(result)
        if record.source == "PanelApp"
    )

    assert (
        panel.value[
            "normalized_inheritance"
        ]
        == "AUTOSOMAL_DOMINANT_OR_MONOALLELIC"
    )


def test_hpo_ad_term_recognized():
    result = interpret_inheritance(
        make_bundle(),
        genotype_state="heterozygous",
    )

    hpo = next(
        record
        for record in source_observations(result)
        if record.source == "HPO"
    )

    assert (
        hpo.value[
            "normalized_inheritance"
        ]
        == "AUTOSOMAL_DOMINANT"
    )


def test_heterozygous_ad_is_compatible():
    result = interpret_inheritance(
        make_bundle(),
        genotype_state="heterozygous",
    )

    assert result.status == "AVAILABLE"

    states = {
        record.value[
            "genotype_compatibility"
        ]
        for record in source_observations(
            result
        )
    }

    assert "COMPATIBLE" in states


def test_missing_genotype_is_partial():
    result = interpret_inheritance(
        make_bundle()
    )

    assert result.status == "PARTIAL"


def test_case_context_preserves_de_novo():
    result = interpret_inheritance(
        make_bundle(),
        genotype_state="heterozygous",
        proband_sex="female",
        de_novo_status=True,
    )

    context = next(
        observation
        for observation in result.observations
        if observation.code
        == "INHERITANCE_CASE_CONTEXT"
    )

    assert context.value[
        "de_novo_status"
    ] is True


def test_invalid_sex_rejected():
    with pytest.raises(ValueError):
        interpret_inheritance(
            make_bundle(),
            genotype_state="heterozygous",
            proband_sex="invalid",
        )


def test_no_inheritance_records_unavailable():
    bundle = {
        "records": [
            {
                "evidence": {
                    "source": "GenCC",
                    "inheritance": [],
                    "phenotype_ids": [],
                }
            }
        ]
    }

    result = interpret_inheritance(
        bundle,
        genotype_state="heterozygous",
    )

    assert result.status == "UNAVAILABLE"


def test_summary_contains_models():
    result = interpret_inheritance(
        make_bundle(),
        genotype_state="heterozygous",
    )

    summary = summary_observation(result)

    assert "AUTOSOMAL_DOMINANT" in (
        summary.value[
            "normalized_models"
        ]
    )

    assert (
        "AUTOSOMAL_DOMINANT_OR_MONOALLELIC"
        in summary.value[
            "normalized_models"
        ]
    )


def test_no_ranking_or_classification():
    result = interpret_inheritance(
        make_bundle(),
        genotype_state="heterozygous",
    )

    data = result.to_dict()

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
