from genomics_platform.interpretation.interpretation_engine import (
    interpret_variant,
)

from test_interpretation_engine import make_contract


def make_gene_context():
    return {
        "schema_version": "1.0",
        "engine": {
            "name": "gene_evidence_bundle",
            "version": "0.2.0",
        },
        "query": {
            "gene_symbol": "OR4F5",
        },
        "scope": {
            "clinical_classification": False,
            "variant_pathogenicity": False,
            "candidate_ranking": False,
        },
        "summary": {
            "record_count": 2,
            "source_counts": {
                "GenCC": 1,
                "PanelApp": 1,
            },
            "disease_normalization": {
                "FOUND": 2,
            },
            "normalized_disease_count": 1,
        },
        "normalized_diseases": [
            {
                "mondo_id": "MONDO:0000001",
                "mondo_label": "test disease",
            }
        ],
        "records": [
            {
                "evidence": {
                    "source": "GenCC",
                    "gene_symbol": "OR4F5",
                    "disease_id": "MONDO:0000001",
                    "disease_label": "test disease",
                    "classification": "Definitive",
                },
                "disease_normalization": {
                    "status": "FOUND",
                    "mondo_id": "MONDO:0000001",
                    "mondo_label": "test disease",
                },
            },
            {
                "evidence": {
                    "source": "PanelApp",
                    "gene_symbol": "OR4F5",
                    "disease_id": "MONDO:0000001",
                    "disease_label": "test disease",
                    "confidence_level": "Green",
                },
                "disease_normalization": {
                    "status": "FOUND",
                    "mondo_id": "MONDO:0000001",
                    "mondo_label": "test disease",
                },
            },
        ],
    }


def test_engine_without_gene_context_preserves_old_behavior():
    profile = interpret_variant(
        make_contract()
    )

    assert profile.gene_disease is None


def test_engine_accepts_gene_evidence_context():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=make_gene_context(),
    )

    assert profile.gene_disease is not None
    assert profile.gene_disease.name == "gene_disease"
    assert profile.gene_disease.status == "AVAILABLE"


def test_engine_serializes_gene_disease_dimension():
    data = interpret_variant(
        make_contract(),
        gene_evidence_context=make_gene_context(),
    ).to_dict()

    gene_disease = data["dimensions"]["gene_disease"]

    assert gene_disease is not None
    assert gene_disease["name"] == "gene_disease"
    assert gene_disease["status"] == "AVAILABLE"


def test_engine_gene_disease_is_available_dimension():
    profile = interpret_variant(
        make_contract(),
        gene_evidence_context=make_gene_context(),
    )

    assert "gene_disease" in profile.available_dimensions


def test_gene_context_does_not_create_variant_classification():
    data = interpret_variant(
        make_contract(),
        gene_evidence_context=make_gene_context(),
    ).to_dict()

    forbidden = {
        "score",
        "ranking_score",
        "tier",
        "acmg_classification",
        "pathogenicity_probability",
        "diagnosis",
    }

    assert forbidden.isdisjoint(data.keys())

    gene_disease = data["dimensions"]["gene_disease"]

    assert forbidden.isdisjoint(
        gene_disease.keys()
    )


def test_gencc_definitive_remains_source_evidence():
    data = interpret_variant(
        make_contract(),
        gene_evidence_context=make_gene_context(),
    ).to_dict()

    records = [
        observation
        for observation in data[
            "dimensions"
        ]["gene_disease"]["observations"]
        if observation["code"]
        == "GENE_DISEASE_SOURCE_RECORD"
    ]

    gencc = next(
        record
        for record in records
        if record["value"]["source"] == "GenCC"
    )

    assert (
        gencc["value"]["evidence"]["classification"]
        == "Definitive"
    )

    assert (
        gencc["details"]["variant_pathogenicity"]
        is False
    )


def test_panelapp_green_remains_gene_level_context():
    data = interpret_variant(
        make_contract(),
        gene_evidence_context=make_gene_context(),
    ).to_dict()

    records = [
        observation
        for observation in data[
            "dimensions"
        ]["gene_disease"]["observations"]
        if observation["code"]
        == "GENE_DISEASE_SOURCE_RECORD"
    ]

    panelapp = next(
        record
        for record in records
        if record["value"]["source"] == "PanelApp"
    )

    assert (
        panelapp["value"]["evidence"]["confidence_level"]
        == "Green"
    )

    assert (
        panelapp["details"]["variant_pathogenicity"]
        is False
    )
