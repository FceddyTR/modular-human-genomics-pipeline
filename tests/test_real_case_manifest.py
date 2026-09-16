import json

import pytest

from genomics_platform.ranking.benchmark.real_cases import (
    load_case_manifest,
    load_truth_manifest,
)


def _write_json(path, payload):
    path.write_text(
        json.dumps(payload),
        encoding="utf-8",
    )


def _valid_case_payload():
    return {
        "schema_version": "0.1",
        "case_id": "public-case-001",
        "assembly": "GRCh38",
        "phenotype": {
            "present_hpo_terms": [
                "HP:0001250",
            ],
            "absent_hpo_terms": [],
        },
        "candidates": [
            {
                "candidate_id": "candidate-a",
                "variant_key": "1:100:A:G",
                "gene_symbol": "GENE1",
            },
            {
                "candidate_id": "candidate-b",
                "variant_key": "1:200:C:T",
                "gene_symbol": "GENE2",
            },
        ],
        "provenance": {
            "source_name": "Public source",
            "source_reference": "reference-001",
            "source_type": "published_case",
        },
    }


def test_case_manifest_loads_without_truth(tmp_path):
    path = tmp_path / "case.json"

    _write_json(
        path,
        _valid_case_payload(),
    )

    payload = load_case_manifest(path)

    assert payload["case_id"] == "public-case-001"
    assert len(payload["candidates"]) == 2


@pytest.mark.parametrize(
    "forbidden_key",
    (
        "truth",
        "causal_candidate_id",
        "causal_candidate_ids",
        "is_causal",
        "expected_rank",
        "expected_ranking",
        "benchmark_label",
    ),
)
def test_case_manifest_rejects_root_truth_keys(
    tmp_path,
    forbidden_key,
):
    payload = _valid_case_payload()
    payload[forbidden_key] = "leak"

    path = tmp_path / "case.json"
    _write_json(path, payload)

    with pytest.raises(
        ValueError,
        match="forbidden",
    ):
        load_case_manifest(path)


def test_case_manifest_rejects_nested_truth_leakage(
    tmp_path,
):
    payload = _valid_case_payload()

    payload["candidates"][0]["metadata"] = {
        "is_causal": True,
    }

    path = tmp_path / "case.json"
    _write_json(path, payload)

    with pytest.raises(
        ValueError,
        match=r"\$\.candidates\[0\]\.metadata\.is_causal",
    ):
        load_case_manifest(path)


def test_case_manifest_rejects_wrong_schema(
    tmp_path,
):
    payload = _valid_case_payload()
    payload["schema_version"] = "999"

    path = tmp_path / "case.json"
    _write_json(path, payload)

    with pytest.raises(
        ValueError,
        match="schema_version",
    ):
        load_case_manifest(path)


def test_truth_manifest_loads_separately(
    tmp_path,
):
    path = tmp_path / "truth.json"

    _write_json(
        path,
        {
            "schema_version": "0.1",
            "case_id": "public-case-001",
            "causal_candidate_ids": [
                "candidate-b",
            ],
        },
    )

    truth = load_truth_manifest(path)

    assert truth.case_id == "public-case-001"
    assert truth.causal_candidate_ids == (
        "candidate-b",
    )


def test_truth_manifest_rejects_extra_fields(
    tmp_path,
):
    path = tmp_path / "truth.json"

    _write_json(
        path,
        {
            "schema_version": "0.1",
            "case_id": "public-case-001",
            "causal_candidate_ids": [
                "candidate-b",
            ],
            "expected_rank": 1,
        },
    )

    with pytest.raises(
        ValueError,
        match="Unexpected",
    ):
        load_truth_manifest(path)
