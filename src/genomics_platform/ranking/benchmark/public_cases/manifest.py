"""Reproducible public benchmark cohort manifest."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .cohort import PublicCaseRecord
from .phenopacket import (
    PHENOPACKET_STORE_RELEASE,
    PUBLIC_CASE_PROTOCOL_VERSION,
)
from .split import PublicBenchmarkSplit


PUBLIC_COHORT_MANIFEST_SCHEMA_VERSION = "0.1"


@dataclass(frozen=True)
class PublicCohortManifest:
    development: tuple[PublicCaseRecord, ...]
    locked_evaluation: tuple[PublicCaseRecord, ...]
    source_release: str = PHENOPACKET_STORE_RELEASE
    protocol_version: str = PUBLIC_CASE_PROTOCOL_VERSION
    schema_version: str = (
        PUBLIC_COHORT_MANIFEST_SCHEMA_VERSION
    )

    def __post_init__(self) -> None:
        # Reuse the leakage contract.
        PublicBenchmarkSplit(
            development=self.development,
            locked_evaluation=self.locked_evaluation,
        )

        if len(self.development) != 300:
            raise ValueError(
                "Public benchmark v0.1 requires "
                "300 development cases."
            )

        if len(self.locked_evaluation) != 200:
            raise ValueError(
                "Public benchmark v0.1 requires "
                "200 locked evaluation cases."
            )

        if (
            self.source_release
            != PHENOPACKET_STORE_RELEASE
        ):
            raise ValueError(
                "Unexpected Phenopacket Store "
                "release."
            )

        if (
            self.protocol_version
            != PUBLIC_CASE_PROTOCOL_VERSION
        ):
            raise ValueError(
                "Unexpected public benchmark "
                "protocol version."
            )

    @property
    def cases(
        self,
    ) -> tuple[PublicCaseRecord, ...]:
        return (
            self.development
            + self.locked_evaluation
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "protocol_version": (
                self.protocol_version
            ),
            "source": {
                "name": (
                    "Monarch Initiative "
                    "Phenopacket Store"
                ),
                "release": self.source_release,
            },
            "cohort": {
                "total_cases": len(self.cases),
                "development_cases": len(
                    self.development
                ),
                "locked_evaluation_cases": len(
                    self.locked_evaluation
                ),
            },
            "development": [
                _record_to_dict(case)
                for case in self.development
            ],
            "locked_evaluation": [
                _record_to_dict(case)
                for case
                in self.locked_evaluation
            ],
            "semantic_boundaries": {
                "contains_benchmark_truth": True,
                "ranking_input": False,
                "truth_used_for_ranking": False,
                "clinical_validation": False,
            },
        }


def _record_to_dict(
    case: PublicCaseRecord,
) -> dict[str, Any]:
    return {
        "case_id": case.case_id,
        "source_reference": (
            case.source_reference
        ),
        "publication_id": (
            case.publication_id
        ),
        "causal_genes": list(
            case.causal_genes
        ),
        "causal_variant_keys": list(
            case.causal_variant_keys
        ),
    }


def build_public_cohort_manifest(
    split: PublicBenchmarkSplit,
) -> PublicCohortManifest:
    return PublicCohortManifest(
        development=split.development,
        locked_evaluation=(
            split.locked_evaluation
        ),
    )
