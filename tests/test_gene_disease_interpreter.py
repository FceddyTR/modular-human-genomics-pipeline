import pytest

from genomics_platform.interpretation.gene_disease_interpreter import (
    interpret_gene_disease,
)


def make_bundle():
    return {
        "schema_version": "1.0",
        "engine": {
            "name": "gene_evidence_bundle",
            "version": "0.2.0",
        },
        "query": {
            "gene_symbol": "TEST1",
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
                    "gene_symbol": "TEST1",
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
                    "gene_symbol": "TEST1",
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


def test_available_when_records_are_present_and_normalized():
    dimension = interpret_gene_disease(make_bundle())

    assert dimension.name == "gene_disease"
    assert dimension.status == "AVAILABLE"


def test_preserves_gene_symbol():
    data = interpret_gene_disease(
        make_bundle()
    ).to_dict()

    observation = next(
        item
        for item in data["observations"]
        if item["code"] == "GENE_DISEASE_GENE"
    )

    assert observation["value"] == "TEST1"


def test_preserves_normalized_disease():
    data = interpret_gene_disease(
        make_bundle()
    ).to_dict()

    observation = next(
        item
        for item in data["observations"]
        if item["code"]
        == "GENE_DISEASE_NORMALIZED_DISEASES"
    )

    assert observation["value"][0]["mondo_id"] == (
        "MONDO:0000001"
    )


def test_preserves_source_specific_records():
    data = interpret_gene_disease(
        make_bundle()
    ).to_dict()

    records = [
        item
        for item in data["observations"]
        if item["code"] == "GENE_DISEASE_SOURCE_RECORD"
    ]

    assert len(records) == 2

    sources = {
        item["value"]["source"]
        for item in records
    }

    assert sources == {
        "GenCC",
        "PanelApp",
    }


def test_gencc_definitive_is_not_variant_classification():
    data = interpret_gene_disease(
        make_bundle()
    ).to_dict()

    gencc = next(
        item
        for item in data["observations"]
        if (
            item["code"] == "GENE_DISEASE_SOURCE_RECORD"
            and item["value"]["source"] == "GenCC"
        )
    )

    assert (
        gencc["value"]["evidence"]["classification"]
        == "Definitive"
    )

    assert gencc["details"]["variant_pathogenicity"] is False
    assert gencc["details"]["candidate_ranking"] is False


def test_panelapp_green_is_not_variant_pathogenicity():
    data = interpret_gene_disease(
        make_bundle()
    ).to_dict()

    panelapp = next(
        item
        for item in data["observations"]
        if (
            item["code"] == "GENE_DISEASE_SOURCE_RECORD"
            and item["value"]["source"] == "PanelApp"
        )
    )

    assert (
        panelapp["value"]["evidence"]["confidence_level"]
        == "Green"
    )

    assert panelapp["details"]["variant_pathogenicity"] is False


def test_no_records_is_unavailable_not_benign():
    bundle = make_bundle()

    bundle["summary"]["record_count"] = 0
    bundle["summary"]["source_counts"] = {}
    bundle["summary"]["normalized_disease_count"] = 0
    bundle["summary"]["disease_normalization"] = {}
    bundle["records"] = []
    bundle["normalized_diseases"] = []

    dimension = interpret_gene_disease(bundle)

    assert dimension.status == "UNAVAILABLE"

    assert any(
        "does not mean" in limitation
        for limitation in dimension.limitations
    )


def test_partial_when_disease_normalization_is_incomplete():
    bundle = make_bundle()

    bundle["summary"]["disease_normalization"] = {
        "FOUND": 1,
        "NOT_FOUND": 1,
    }

    bundle["records"][1]["disease_normalization"] = {
        "status": "NOT_FOUND",
        "mondo_id": None,
        "mondo_label": None,
    }

    dimension = interpret_gene_disease(bundle)

    assert dimension.status == "PARTIAL"

    data = dimension.to_dict()

    panelapp = next(
        item
        for item in data["observations"]
        if (
            item["code"] == "GENE_DISEASE_SOURCE_RECORD"
            and item["value"]["source"] == "PanelApp"
        )
    )

    assert (
        panelapp["value"]["normalization_status"]
        == "NOT_FOUND"
    )

    assert (
        panelapp["value"]["evidence"]["disease_id"]
        == "MONDO:0000001"
    )


def test_none_bundle_is_unavailable():
    dimension = interpret_gene_disease(None)

    assert dimension.status == "UNAVAILABLE"


def test_wrong_input_type_rejected():
    with pytest.raises(TypeError):
        interpret_gene_disease(
            ["not", "a", "bundle"]
        )


def test_no_ranking_or_classification_output():
    data = interpret_gene_disease(
        make_bundle()
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
