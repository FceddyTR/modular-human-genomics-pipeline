from genomics_platform.interpretation.interpretation_engine import (
    interpret_variant,
)

from test_interpretation_engine import make_contract


def make_gene_context():
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


def test_old_behavior_has_no_inheritance():
    profile = interpret_variant(
        make_contract()
    )

    assert profile.inheritance is None


def test_gene_context_alone_does_not_request_inheritance():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=make_gene_context(),
    )

    assert profile.inheritance is None


def test_genotype_requests_inheritance():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=make_gene_context(),
        genotype_state="heterozygous",
    )

    assert profile.inheritance is not None
    assert profile.inheritance.status == "AVAILABLE"


def test_inheritance_without_gene_context_unavailable():
    profile = interpret_variant(
        make_contract(),
        genotype_state="heterozygous",
    )

    assert profile.inheritance is not None
    assert profile.inheritance.status == "UNAVAILABLE"


def test_proband_sex_requests_inheritance():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=make_gene_context(),
        proband_sex="female",
    )

    assert profile.inheritance is not None
    assert profile.inheritance.status == "PARTIAL"


def test_de_novo_status_requests_inheritance():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=make_gene_context(),
        de_novo_status=True,
    )

    assert profile.inheritance is not None


def test_inheritance_serializes():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=make_gene_context(),
        genotype_state="heterozygous",
    )

    data = profile.to_dict()

    assert data["dimensions"]["inheritance"] is not None

    assert (
        data["dimensions"]["inheritance"]["status"]
        == "AVAILABLE"
    )


def test_inheritance_in_available_dimensions():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=make_gene_context(),
        genotype_state="heterozygous",
    )

    assert (
        "inheritance"
        in profile.available_dimensions
    )


def test_ncstn_ad_is_compatible():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=make_gene_context(),
        genotype_state="heterozygous",
    )

    records = [
        observation
        for observation
        in profile.inheritance.observations
        if observation.code
        == "INHERITANCE_SOURCE_EVIDENCE"
    ]

    assert any(
        record.value[
            "normalized_inheritance"
        ] == "AUTOSOMAL_DOMINANT"
        and record.value[
            "genotype_compatibility"
        ] == "COMPATIBLE"
        for record in records
    )


def test_phenotype_and_inheritance_can_coexist():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=make_gene_context(),
        present_hpo_terms=[
            "HP:0000006"
        ],
        genotype_state="heterozygous",
    )

    assert profile.phenotype is not None
    assert profile.inheritance is not None


def test_no_ranking_or_classification_added():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=make_gene_context(),
        genotype_state="heterozygous",
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

    assert forbidden.isdisjoint(data.keys())

    assert forbidden.isdisjoint(
        data["dimensions"]["inheritance"].keys()
    )
