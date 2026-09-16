"""Phenopacket Store source ingestion for benchmark construction."""

from __future__ import annotations

import json
import re
from pathlib import Path

from .cohort import PublicCaseRecord
from .phenopacket import evaluate_phenopacket


_PMID_RE = re.compile(
    r"PMID[_:]?(\d+)",
    re.IGNORECASE,
)


def publication_id_from_path(
    path: Path,
) -> str:
    match = _PMID_RE.search(path.name)

    if match:
        return f"PMID:{match.group(1)}"

    return f"SOURCE:{path.parent.name}"


def load_eligible_source_records(
    root: str | Path,
) -> tuple[PublicCaseRecord, ...]:
    """Load strict eligible records directly from a frozen release."""

    root = Path(root)

    if not root.is_dir():
        raise ValueError(
            f"Source root does not exist: {root}"
        )

    records: list[
        PublicCaseRecord
    ] = []

    for path in sorted(
        root.rglob("*.json")
    ):
        payload = json.loads(
            path.read_text(
                encoding="utf-8"
            )
        )

        eligibility = (
            evaluate_phenopacket(
                payload
            )
        )

        if not eligibility.eligible:
            continue

        records.append(
            PublicCaseRecord(
                eligibility=eligibility,
                publication_id=(
                    publication_id_from_path(
                        path
                    )
                ),
                source_reference=str(
                    path.relative_to(root)
                ),
            )
        )

    return tuple(records)
