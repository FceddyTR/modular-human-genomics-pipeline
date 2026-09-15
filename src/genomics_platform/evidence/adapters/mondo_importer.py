#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple


IMPORTER_VERSION = "0.3.0"


def normalize_curie(value: Optional[str]) -> Optional[str]:
    if not value:
        return None

    value = value.strip()

    if not value:
        return None

    if value.startswith("http://purl.obolibrary.org/obo/"):
        tail = value.rsplit("/", 1)[-1]

        if "_" in tail:
            prefix, local = tail.split("_", 1)
            return f"{prefix}:{local}"

    return value


def extract_label(node: Dict[str, Any]) -> Optional[str]:
    lbl = node.get("lbl")

    if isinstance(lbl, str) and lbl.strip():
        return lbl.strip()

    meta = node.get("meta") or {}
    label = meta.get("label")

    if isinstance(label, str) and label.strip():
        return label.strip()

    return None


def iter_xrefs(node: Dict[str, Any]) -> Iterable[str]:
    meta = node.get("meta") or {}
    xrefs = meta.get("xrefs") or []

    for xref in xrefs:
        if isinstance(xref, str):
            value = xref
        elif isinstance(xref, dict):
            value = xref.get("val") or xref.get("id")
        else:
            value = None

        value = normalize_curie(value)

        if value:
            yield value


def iter_mondo_nodes(payload: Dict[str, Any]) -> Iterable[Dict[str, Any]]:
    for graph in payload.get("graphs") or []:
        for node in graph.get("nodes") or []:
            if isinstance(node, dict):
                yield node


def build_database(
    input_json: str,
    database_path: str,
) -> Dict[str, int]:

    database = Path(database_path)
    database.parent.mkdir(parents=True, exist_ok=True)

    if database.exists():
        database.unlink()

    with open(input_json, encoding="utf-8") as handle:
        payload = json.load(handle)

    con = sqlite3.connect(database)

    try:
        con.executescript(
            """
            CREATE TABLE metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );

            CREATE TABLE disease (
                mondo_id TEXT PRIMARY KEY,
                label TEXT
            );

            CREATE TABLE xref (
                source_id TEXT NOT NULL,
                mondo_id TEXT NOT NULL,
                PRIMARY KEY (source_id, mondo_id)
            );

            CREATE INDEX idx_xref_source
            ON xref(source_id);

            CREATE INDEX idx_xref_mondo
            ON xref(mondo_id);
            """
        )

        con.execute(
            "INSERT INTO metadata(key, value) VALUES (?, ?)",
            ("importer_version", IMPORTER_VERSION),
        )

        con.execute(
            "INSERT INTO metadata(key, value) VALUES (?, ?)",
            ("source", "MONDO"),
        )

        for node in iter_mondo_nodes(payload):
            mondo_id = normalize_curie(node.get("id"))

            if not mondo_id or not mondo_id.startswith("MONDO:"):
                continue

            label = extract_label(node)

            con.execute(
                """
                INSERT OR REPLACE INTO disease(mondo_id, label)
                VALUES (?, ?)
                """,
                (mondo_id, label),
            )

            for source_id in iter_xrefs(node):
                if source_id.startswith("MONDO:"):
                    continue

                con.execute(
                    """
                    INSERT OR IGNORE INTO xref(source_id, mondo_id)
                    VALUES (?, ?)
                    """,
                    (source_id, mondo_id),
                )

        con.commit()

        disease_count = con.execute(
            "SELECT COUNT(*) FROM disease"
        ).fetchone()[0]

        xref_count = con.execute(
            "SELECT COUNT(*) FROM xref"
        ).fetchone()[0]

        source_count = con.execute(
            "SELECT COUNT(DISTINCT source_id) FROM xref"
        ).fetchone()[0]

        return {
            "unique_mondo_diseases": disease_count,
            "unique_xrefs": xref_count,
            "unique_source_ids": source_count,
        }

    finally:
        con.close()


def lookup(
    database_path: str,
    disease_id: str,
) -> Tuple[str, list]:

    query_id = normalize_curie(disease_id)

    con = sqlite3.connect(database_path)

    try:
        if query_id and query_id.startswith("MONDO:"):
            row = con.execute(
                """
                SELECT mondo_id, label
                FROM disease
                WHERE mondo_id = ?
                """,
                (query_id,),
            ).fetchone()

            if not row:
                return "NOT_FOUND", []

            return "FOUND", [
                {
                    "input_id": query_id,
                    "mondo_id": row[0],
                    "label": row[1],
                }
            ]

        rows = con.execute(
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
            """,
            (query_id,),
        ).fetchall()

        if not rows:
            return "NOT_FOUND", []

        return "FOUND", [
            {
                "input_id": row[0],
                "mondo_id": row[1],
                "label": row[2],
            }
            for row in rows
        ]

    finally:
        con.close()


def main():
    parser = argparse.ArgumentParser(
        description="Build/query MONDO disease normalization store."
    )

    sub = parser.add_subparsers(
        dest="command",
        required=True,
    )

    build = sub.add_parser("build")
    build.add_argument("--input-json", required=True)
    build.add_argument("--database", required=True)

    query = sub.add_parser("query")
    query.add_argument("--database", required=True)
    query.add_argument("--disease-id", required=True)

    args = parser.parse_args()

    if args.command == "build":
        summary = build_database(
            args.input_json,
            args.database,
        )

        print("MONDO normalizer v0.3")
        print(
            "  unique MONDO diseases:",
            summary["unique_mondo_diseases"],
        )
        print(
            "  unique xrefs:",
            summary["unique_xrefs"],
        )
        print(
            "  unique source IDs:",
            summary["unique_source_ids"],
        )
        return

    status, records = lookup(
        args.database,
        args.disease_id,
    )

    print("Disease:", args.disease_id)
    print("Status:", status)
    print("Matches:", len(records))

    for record in records:
        print(
            "-",
            record["mondo_id"],
            "|",
            record["label"],
        )


if __name__ == "__main__":
    main()
