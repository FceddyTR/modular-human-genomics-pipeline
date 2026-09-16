"""Frozen production-state snapshots for real-case benchmarks."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date
from typing import Any

from genomics_platform.interpretation.case_model import (
    CandidateCase,
)


SNAPSHOT_SCHEMA_VERSION = "0.1"


@dataclass(frozen=True)
class EvidenceSnapshotMetadata:
    """Metadata required to audit a frozen evidence snapshot."""

    snapshot_date: str
    pipeline_version: str
    source_versions: dict[str, str]

    def __post_init__(self) -> None:
        snapshot_date = self.snapshot_date.strip()
        pipeline_version = self.pipeline_version.strip()

        if not snapshot_date:
            raise ValueError(
                "snapshot_date is required."
            )

        try:
            date.fromisoformat(snapshot_date)
        except ValueError as exc:
            raise ValueError(
                "snapshot_date must use ISO YYYY-MM-DD format."
            ) from exc

        if not pipeline_version:
            raise ValueError(
                "pipeline_version is required."
            )

        normalized_versions: dict[str, str] = {}

        for source, version in self.source_versions.items():
            source = str(source).strip()
            version = str(version).strip()

            if not source or not version:
                raise ValueError(
                    "source_versions cannot contain "
                    "empty source names or versions."
                )

            normalized_versions[source] = version

        object.__setattr__(
            self,
            "snapshot_date",
            snapshot_date,
        )
        object.__setattr__(
            self,
            "pipeline_version",
            pipeline_version,
        )
        object.__setattr__(
            self,
            "source_versions",
            normalized_versions,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "snapshot_schema_version": (
                SNAPSHOT_SCHEMA_VERSION
            ),
            "snapshot_date": self.snapshot_date,
            "pipeline_version": self.pipeline_version,
            "source_versions": dict(
                sorted(
                    self.source_versions.items()
                )
            ),
        }


def snapshot_candidate(
    candidate: CandidateCase,
    metadata: EvidenceSnapshotMetadata,
) -> dict[str, Any]:
    """Serialize one production CandidateCase without benchmark truth."""

    if not isinstance(
        candidate,
        CandidateCase,
    ):
        raise TypeError(
            "candidate must be a CandidateCase."
        )

    return {
        "snapshot_metadata": metadata.to_dict(),
        "candidate": candidate.to_dict(),
        "semantic_boundaries": {
            "frozen_production_state": True,
            "benchmark_truth_present": False,
            "truth_used_for_ranking": False,
            "clinical_validation": False,
        },
    }
