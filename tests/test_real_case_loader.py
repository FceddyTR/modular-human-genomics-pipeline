import json

import pytest

from genomics_platform.interpretation.case_model import (
    CandidateCase,
)
from genomics_platform.ranking.benchmark.real_cases import (
    load_real_case_input,
)
from genomics_platform.ranking.benchmark.real_cases.snapshot import (
    EvidenceSnapshotMetadata,
    snapshot_candidate,
)
from test_case_model import make_candidate


def _case_payload():
    candidate = make_candidate()

    metadata = EvidenceSnapshotMetadata(
        snapshot_date="2026-09-16",
        pipeline_version="6e1a5e0",
        source_versions={
            "ClinVar": "test",
            "gnomAD": "test",
        },
    )

    return {
        "schema_version": "0.1",
        "case_id": "public-case-001",
        "assembly": "GRCh38",
        "phenotype": {
            "present_hpo_terms": list(
                candidate.case_context
                .present_hpo_terms
            ),
            "absent_hpo_terms": list(
                candidate.case_context
                .absent_hpo_terms
            ),
        },
        "candidates": [
            snapshot_candidate(
                candidate,
                metadata,
            ),
        ],
        "provenance": {
            "source_name": "Public source",
            "source_reference": "reference-001",
            "source_type": "published_case",
            "accessed_date": "2026-09-16",
        },
    }


def _write_case(tmp_path, payload):
    path = tmp_path / "case.json"

    path.write_text(
        json.dumps(payload),
        encoding="utf-8",
    )

    return path


def test_load_real_case_input_end_to_end(
    tmp_path,
):
    path = _write_case(
        tmp_path,
        _case_payload(),
    )

    case = load_real_case_input(path)

    assert case.case_id == "public-case-001"
    assert case.assembly == "GRCh38"

    assert len(case.candidates) == 1

    assert isinstance(
        case.candidates[0],
        CandidateCase,
    )

    assert (
        case.candidates[0].to_dict()
        == make_candidate().to_dict()
    )


def test_real_case_loader_has_no_truth_parameter():
    import inspect

    signature = inspect.signature(
        load_real_case_input
    )

    assert tuple(
        signature.parameters
    ) == ("path",)


def test_loader_rejects_truth_leakage(
    tmp_path,
):
    payload = _case_payload()

    payload["candidates"][0][
        "candidate"
    ][
        "is_causal"
    ] = True

    path = _write_case(
        tmp_path,
        payload,
    )

    with pytest.raises(
        ValueError,
        match="forbidden",
    ):
        load_real_case_input(path)


def test_loader_rejects_assembly_mismatch(
    tmp_path,
):
    payload = _case_payload()

    payload["assembly"] = "GRCh37"

    path = _write_case(
        tmp_path,
        payload,
    )

    with pytest.raises(
        ValueError,
        match="assembly",
    ):
        load_real_case_input(path)


def test_loader_rejects_present_phenotype_mismatch(
    tmp_path,
):
    payload = _case_payload()

    payload["phenotype"][
        "present_hpo_terms"
    ] = ["HP:0001250"]

    path = _write_case(
        tmp_path,
        payload,
    )

    with pytest.raises(
        ValueError,
        match="present HPO",
    ):
        load_real_case_input(path)


def test_loader_rejects_absent_phenotype_mismatch(
    tmp_path,
):
    payload = _case_payload()

    payload["phenotype"][
        "absent_hpo_terms"
    ] = ["HP:0001263"]

    path = _write_case(
        tmp_path,
        payload,
    )

    with pytest.raises(
        ValueError,
        match="absent HPO",
    ):
        load_real_case_input(path)


def test_loader_rejects_direct_unsnapshotted_candidate(
    tmp_path,
):
    payload = _case_payload()

    payload["candidates"] = [
        make_candidate().to_dict()
    ]

    path = _write_case(
        tmp_path,
        payload,
    )

    with pytest.raises(
        ValueError,
        match="snapshot_metadata",
    ):
        load_real_case_input(path)


def test_loader_rejects_snapshot_truth_boundary(
    tmp_path,
):
    payload = _case_payload()

    payload["candidates"][0][
        "semantic_boundaries"
    ][
        "benchmark_truth_present"
    ] = True

    path = _write_case(
        tmp_path,
        payload,
    )

    with pytest.raises(
        ValueError,
        match="benchmark_truth_present",
    ):
        load_real_case_input(path)


def test_loader_rejects_snapshot_without_pipeline_version(
    tmp_path,
):
    payload = _case_payload()

    payload["candidates"][0][
        "snapshot_metadata"
    ][
        "pipeline_version"
    ] = ""

    path = _write_case(
        tmp_path,
        payload,
    )

    with pytest.raises(
        ValueError,
        match="pipeline_version",
    ):
        load_real_case_input(path)
