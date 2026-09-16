"""Load truth-blind frozen real-case benchmark inputs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .contract import (
    CaseProvenance,
    RealCaseInput,
)
from .manifest import load_case_manifest
from .replay import replay_candidate_case
from .snapshot import SNAPSHOT_SCHEMA_VERSION
from .validation import validate_real_case_input


def _mapping(
    value: Any,
    *,
    field_name: str,
) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(
            f"{field_name} must be an object."
        )
    return value


def _sequence(
    value: Any,
    *,
    field_name: str,
) -> list[Any]:
    if not isinstance(value, list):
        raise ValueError(
            f"{field_name} must be an array."
        )
    return value


def _string_list(
    value: Any,
    *,
    field_name: str,
) -> tuple[str, ...]:
    values = _sequence(
        value,
        field_name=field_name,
    )

    if not all(
        isinstance(item, str)
        for item in values
    ):
        raise ValueError(
            f"{field_name} must contain strings."
        )

    return tuple(values)


def _load_provenance(
    payload: Any,
) -> CaseProvenance:
    payload = _mapping(
        payload,
        field_name="provenance",
    )

    return CaseProvenance(
        source_name=str(
            payload.get(
                "source_name",
                "",
            )
        ),
        source_reference=str(
            payload.get(
                "source_reference",
                "",
            )
        ),
        source_type=str(
            payload.get(
                "source_type",
                "",
            )
        ),
        accessed_date=payload.get(
            "accessed_date"
        ),
        notes=payload.get("notes"),
    )


def _load_snapshot_candidate(
    payload: Any,
):
    wrapper = _mapping(
        payload,
        field_name="candidate_snapshot",
    )

    metadata = _mapping(
        wrapper.get("snapshot_metadata"),
        field_name=(
            "candidate_snapshot."
            "snapshot_metadata"
        ),
    )

    if (
        metadata.get(
            "snapshot_schema_version"
        )
        != SNAPSHOT_SCHEMA_VERSION
    ):
        raise ValueError(
            "Unsupported candidate snapshot "
            "schema version."
        )

    snapshot_date = metadata.get(
        "snapshot_date"
    )
    pipeline_version = metadata.get(
        "pipeline_version"
    )
    source_versions = metadata.get(
        "source_versions"
    )

    if not isinstance(
        snapshot_date,
        str,
    ) or not snapshot_date.strip():
        raise ValueError(
            "Candidate snapshot requires "
            "snapshot_date."
        )

    if not isinstance(
        pipeline_version,
        str,
    ) or not pipeline_version.strip():
        raise ValueError(
            "Candidate snapshot requires "
            "pipeline_version."
        )

    if not isinstance(
        source_versions,
        dict,
    ):
        raise ValueError(
            "Candidate snapshot requires "
            "source_versions object."
        )

    candidate_payload = _mapping(
        wrapper.get("candidate"),
        field_name=(
            "candidate_snapshot.candidate"
        ),
    )

    boundaries = _mapping(
        wrapper.get(
            "semantic_boundaries",
            {},
        ),
        field_name=(
            "candidate_snapshot."
            "semantic_boundaries"
        ),
    )

    if (
        boundaries.get(
            "benchmark_truth_present"
        )
        is not False
    ):
        raise ValueError(
            "Candidate snapshot must declare "
            "benchmark_truth_present=false."
        )

    if (
        boundaries.get(
            "truth_used_for_ranking"
        )
        is not False
    ):
        raise ValueError(
            "Candidate snapshot must declare "
            "truth_used_for_ranking=false."
        )

    return replay_candidate_case(
        candidate_payload
    )


def load_real_case_input(
    path: str | Path,
) -> RealCaseInput:
    """Load one frozen case without accessing benchmark truth."""

    payload = load_case_manifest(path)

    phenotype = _mapping(
        payload.get(
            "phenotype",
            {},
        ),
        field_name="phenotype",
    )

    candidate_payloads = _sequence(
        payload.get(
            "candidates",
            [],
        ),
        field_name="candidates",
    )

    candidates = tuple(
        _load_snapshot_candidate(item)
        for item in candidate_payloads
    )

    case = RealCaseInput(
        case_id=str(
            payload.get(
                "case_id",
                "",
            )
        ),
        assembly=str(
            payload.get(
                "assembly",
                "",
            )
        ),
        present_hpo_terms=_string_list(
            phenotype.get(
                "present_hpo_terms",
                [],
            ),
            field_name=(
                "phenotype.present_hpo_terms"
            ),
        ),
        absent_hpo_terms=_string_list(
            phenotype.get(
                "absent_hpo_terms",
                [],
            ),
            field_name=(
                "phenotype.absent_hpo_terms"
            ),
        ),
        candidates=candidates,
        provenance=_load_provenance(
            payload.get("provenance")
        ),
    )

    validate_real_case_input(case)

    for candidate in case.candidates:
        if (
            candidate.evidence_contract
            .identity
            .assembly
            != case.assembly
        ):
            raise ValueError(
                "Candidate assembly does not match "
                "case assembly: "
                f"{candidate.candidate_id!r}"
            )

        if (
            candidate.case_context
            .present_hpo_terms
            != case.present_hpo_terms
        ):
            raise ValueError(
                "Candidate present HPO context does "
                "not match case phenotype: "
                f"{candidate.candidate_id!r}"
            )

        if (
            candidate.case_context
            .absent_hpo_terms
            != case.absent_hpo_terms
        ):
            raise ValueError(
                "Candidate absent HPO context does "
                "not match case phenotype: "
                f"{candidate.candidate_id!r}"
            )

    return case
