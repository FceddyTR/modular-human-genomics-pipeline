from dataclasses import dataclass

import pytest

from genomics_platform.ranking.benchmark.real_cases.snapshot import (
    EvidenceSnapshotMetadata,
    snapshot_candidate,
)


@dataclass(frozen=True)
class NotCandidate:
    candidate_id: str = "fake"


def test_snapshot_metadata_normalizes_versions():
    metadata = EvidenceSnapshotMetadata(
        snapshot_date="2026-09-16",
        pipeline_version="1d1b8f7",
        source_versions={
            " ClinVar ": " 2026-09 ",
            "gnomAD": "4.1",
        },
    )

    assert metadata.source_versions == {
        "ClinVar": "2026-09",
        "gnomAD": "4.1",
    }


def test_snapshot_metadata_rejects_invalid_date():
    with pytest.raises(
        ValueError,
        match="ISO",
    ):
        EvidenceSnapshotMetadata(
            snapshot_date="16/09/2026",
            pipeline_version="1d1b8f7",
            source_versions={},
        )


def test_snapshot_metadata_requires_pipeline_version():
    with pytest.raises(
        ValueError,
        match="pipeline_version",
    ):
        EvidenceSnapshotMetadata(
            snapshot_date="2026-09-16",
            pipeline_version=" ",
            source_versions={},
        )


def test_snapshot_metadata_rejects_empty_source_version():
    with pytest.raises(
        ValueError,
        match="source_versions",
    ):
        EvidenceSnapshotMetadata(
            snapshot_date="2026-09-16",
            pipeline_version="1d1b8f7",
            source_versions={
                "ClinVar": "",
            },
        )


def test_snapshot_candidate_rejects_non_candidate():
    metadata = EvidenceSnapshotMetadata(
        snapshot_date="2026-09-16",
        pipeline_version="1d1b8f7",
        source_versions={},
    )

    with pytest.raises(
        TypeError,
        match="CandidateCase",
    ):
        snapshot_candidate(
            NotCandidate(),
            metadata,
        )
