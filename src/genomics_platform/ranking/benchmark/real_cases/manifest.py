"""Disk manifest contracts for truth-blind real-case benchmarks."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from genomics_platform.ranking.benchmark.truth_set import (
    RankingTruth,
)

from .contract import REAL_CASE_SCHEMA_VERSION


FORBIDDEN_CASE_KEYS = frozenset(
    {
        "truth",
        "causal_candidate_id",
        "causal_candidate_ids",
        "is_causal",
        "expected_rank",
        "expected_ranking",
        "benchmark_label",
    }
)


def _load_json(path: str | Path) -> dict[str, Any]:
    path = Path(path)

    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        payload = json.load(handle)

    if not isinstance(payload, dict):
        raise ValueError(
            f"Manifest must contain a JSON object: {path}"
        )

    return payload


def _find_forbidden_keys(
    value: Any,
    *,
    path: str = "$",
) -> tuple[str, ...]:
    """Find forbidden truth-bearing keys recursively."""

    found: list[str] = []

    if isinstance(value, dict):
        for key, child in value.items():
            child_path = f"{path}.{key}"

            if key in FORBIDDEN_CASE_KEYS:
                found.append(child_path)

            found.extend(
                _find_forbidden_keys(
                    child,
                    path=child_path,
                )
            )

    elif isinstance(value, list):
        for index, child in enumerate(value):
            found.extend(
                _find_forbidden_keys(
                    child,
                    path=f"{path}[{index}]",
                )
            )

    return tuple(found)


def load_case_manifest(
    path: str | Path,
) -> dict[str, Any]:
    """Load ranking-time manifest while rejecting truth leakage."""

    payload = _load_json(path)

    leaked = _find_forbidden_keys(payload)

    if leaked:
        raise ValueError(
            "Truth-bearing keys are forbidden from "
            "case manifests: "
            f"{list(leaked)!r}"
        )

    schema_version = str(
        payload.get("schema_version", "")
    ).strip()

    if schema_version != REAL_CASE_SCHEMA_VERSION:
        raise ValueError(
            "Unsupported real-case schema_version: "
            f"{schema_version!r}"
        )

    case_id = str(
        payload.get("case_id", "")
    ).strip()

    if not case_id:
        raise ValueError(
            "case_id is required."
        )

    return payload


def load_truth_manifest(
    path: str | Path,
) -> RankingTruth:
    """Load benchmark truth separately from ranking-time input."""

    payload = _load_json(path)

    schema_version = str(
        payload.get("schema_version", "")
    ).strip()

    if schema_version != REAL_CASE_SCHEMA_VERSION:
        raise ValueError(
            "Unsupported truth schema_version: "
            f"{schema_version!r}"
        )

    allowed_keys = {
        "schema_version",
        "case_id",
        "causal_candidate_ids",
    }

    unexpected = (
        set(payload)
        - allowed_keys
    )

    if unexpected:
        raise ValueError(
            "Unexpected truth manifest keys: "
            f"{sorted(unexpected)!r}"
        )

    causal = payload.get(
        "causal_candidate_ids",
        (),
    )

    if not isinstance(causal, list):
        raise ValueError(
            "causal_candidate_ids must be a JSON array."
        )

    return RankingTruth(
        case_id=str(
            payload.get("case_id", "")
        ),
        causal_candidate_ids=tuple(
            str(candidate_id)
            for candidate_id in causal
        ),
    )
