#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from typing import Dict, Iterable, List, Optional

from genomics_platform.evidence.disease_normalizer import (
    normalize_disease,
)


ADAPTER_VERSION = "0.1.0"


def canonicalize_clinvar_identifier(
    identifier: Optional[str],
) -> Optional[str]:
    """
    Apply only safe, source-specific identifier cleanup.

    Raw ClinVar identifiers are never modified in the ClinVarRecord.
    This function produces a derived identifier used for MONDO lookup.

    Example:
        MONDO:MONDO:0019200 -> MONDO:0019200

    No fuzzy disease-name matching is performed.
    """

    if identifier is None:
        return None

    identifier = identifier.strip()

    if not identifier:
        return None

    if identifier.startswith("MONDO:MONDO:"):
        return "MONDO:" + identifier.split(
            "MONDO:MONDO:",
            1,
        )[1]

    return identifier


def normalize_identifiers(
    database_path: str,
    identifiers: Iterable[str],
) -> List[Dict]:
    results: List[Dict] = []

    for raw_identifier in identifiers:
        lookup_identifier = (
            canonicalize_clinvar_identifier(
                raw_identifier
            )
        )

        if not lookup_identifier:
            continue

        normalized = normalize_disease(
            database_path,
            lookup_identifier,
        )

        results.append(
            {
                "raw_identifier": raw_identifier,
                "lookup_identifier": lookup_identifier,
                "status": normalized["status"],
                "matched_source_id": normalized[
                    "matched_source_id"
                ],
                "mondo_id": normalized["mondo_id"],
                "mondo_label": normalized[
                    "mondo_label"
                ],
            }
        )

    return results


def summarize_condition_normalization(
    normalized_identifiers: List[Dict],
) -> Dict:
    """
    Summarize identifier-level MONDO matches without implying
    variant-disease causality.

    Multiple identifiers may resolve to the same MONDO concept.
    If identifiers resolve to different MONDO concepts, the result
    is marked AMBIGUOUS rather than selecting one silently.
    """

    found = [
        item
        for item in normalized_identifiers
        if item["status"] == "FOUND"
        and item["mondo_id"]
    ]

    mondo_ids = list(
        dict.fromkeys(
            item["mondo_id"]
            for item in found
        )
    )

    if not mondo_ids:
        return {
            "status": "NOT_FOUND",
            "mondo_id": None,
            "mondo_label": None,
            "matched_via": [],
        }

    if len(mondo_ids) > 1:
        return {
            "status": "AMBIGUOUS",
            "mondo_id": None,
            "mondo_label": None,
            "candidate_mondo_ids": mondo_ids,
            "matched_via": [
                item["raw_identifier"]
                for item in found
            ],
        }

    mondo_id = mondo_ids[0]

    label = next(
        (
            item["mondo_label"]
            for item in found
            if item["mondo_id"] == mondo_id
        ),
        None,
    )

    return {
        "status": "FOUND",
        "mondo_id": mondo_id,
        "mondo_label": label,
        "matched_via": [
            item["raw_identifier"]
            for item in found
            if item["mondo_id"] == mondo_id
        ],
    }


def normalize_clinvar_condition(
    database_path: str,
    condition: Dict,
) -> Dict:
    raw_name = condition.get("name")
    raw_identifiers = condition.get(
        "identifiers",
        [],
    )

    identifier_results = normalize_identifiers(
        database_path,
        raw_identifiers,
    )

    summary = summarize_condition_normalization(
        identifier_results
    )

    return {
        "raw_name": raw_name,
        "raw_identifiers": list(
            raw_identifiers
        ),
        "normalization": summary,
        "identifier_results": identifier_results,
        "adapter_version": ADAPTER_VERSION,
        "semantics": (
            "Disease identifier normalization only; "
            "does not establish variant-disease causality."
        ),
    }


def normalize_clinvar_record_conditions(
    database_path: str,
    clinvar_record: Dict,
) -> List[Dict]:
    return [
        normalize_clinvar_condition(
            database_path,
            condition,
        )
        for condition in clinvar_record.get(
            "conditions",
            [],
        )
    ]


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Normalize ClinVar condition identifiers "
            "to MONDO concepts."
        )
    )

    parser.add_argument(
        "--database",
        required=True,
    )

    parser.add_argument(
        "--condition-name",
        required=True,
    )

    parser.add_argument(
        "--identifier",
        action="append",
        default=[],
    )

    args = parser.parse_args()

    result = normalize_clinvar_condition(
        args.database,
        {
            "name": args.condition_name,
            "identifiers": args.identifier,
        },
    )

    print(
        json.dumps(
            result,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
