#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sqlite3
from typing import Dict, List, Optional


def clean_id(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None

    value = value.strip()

    return value if value else None


def lookup_candidates(disease_id: str) -> List[str]:
    """
    Generate safe identifier aliases.

    This performs identifier normalization only.
    It does NOT perform fuzzy disease-name matching.
    """

    disease_id = clean_id(disease_id)

    if not disease_id:
        return []

    candidates = [disease_id]

    if disease_id.startswith("ORPHA:"):
        local_id = disease_id.split(":", 1)[1]
        candidates.append(f"Orphanet:{local_id}")

    elif disease_id.startswith("Orphanet:"):
        local_id = disease_id.split(":", 1)[1]
        candidates.append(f"ORPHA:{local_id}")

    # Preserve order while removing duplicates.
    return list(dict.fromkeys(candidates))


def normalize_disease(
    database_path: str,
    disease_id: Optional[str],
) -> Dict:

    disease_id = clean_id(disease_id)

    if not disease_id:
        return {
            "status": "NO_DISEASE_ID",
            "input_id": None,
            "matched_source_id": None,
            "mondo_id": None,
            "mondo_label": None,
        }

    con = sqlite3.connect(database_path)

    try:
        if disease_id.startswith("MONDO:"):
            row = con.execute(
                """
                SELECT mondo_id, label
                FROM disease
                WHERE mondo_id = ?
                """,
                (disease_id,),
            ).fetchone()

            if not row:
                return {
                    "status": "NOT_FOUND",
                    "input_id": disease_id,
                    "matched_source_id": None,
                    "mondo_id": None,
                    "mondo_label": None,
                }

            return {
                "status": "FOUND",
                "input_id": disease_id,
                "matched_source_id": disease_id,
                "mondo_id": row[0],
                "mondo_label": row[1],
            }

        for candidate in lookup_candidates(disease_id):
            row = con.execute(
                """
                SELECT
                    x.source_id,
                    d.mondo_id,
                    d.label
                FROM xref x
                JOIN disease d
                  ON d.mondo_id = x.mondo_id
                WHERE x.source_id = ?
                ORDER BY d.mondo_id
                LIMIT 1
                """,
                (candidate,),
            ).fetchone()

            if row:
                return {
                    "status": "FOUND",
                    "input_id": disease_id,
                    "matched_source_id": row[0],
                    "mondo_id": row[1],
                    "mondo_label": row[2],
                }

        return {
            "status": "NOT_FOUND",
            "input_id": disease_id,
            "matched_source_id": None,
            "mondo_id": None,
            "mondo_label": None,
        }

    finally:
        con.close()


def main():
    parser = argparse.ArgumentParser(
        description="Normalize external disease identifiers to MONDO."
    )

    parser.add_argument(
        "--database",
        required=True,
    )

    parser.add_argument(
        "--disease-id",
        required=True,
    )

    args = parser.parse_args()

    result = normalize_disease(
        args.database,
        args.disease_id,
    )

    print(
        json.dumps(
            result,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
