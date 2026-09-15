#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path
from typing import Dict, List

from src.genomics_platform.evidence.adapters.orphanet_importer import (
    IMPORTER_VERSION,
    iter_orphanet_records,
)


STORE_VERSION = "0.1.0"


def build_database(
    input_xml: str,
    database_path: str,
) -> Dict[str, int]:
    database = Path(database_path)
    database.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if database.exists():
        database.unlink()

    con = sqlite3.connect(database)

    try:
        con.executescript(
            """
            CREATE TABLE metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );

            CREATE TABLE evidence (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                gene_symbol TEXT NOT NULL,
                disease_id TEXT,
                disease_name TEXT,
                record_json TEXT NOT NULL
            );

            CREATE INDEX idx_orphanet_gene
            ON evidence(gene_symbol);

            CREATE INDEX idx_orphanet_disease
            ON evidence(disease_id);
            """
        )

        metadata = {
            "source": "Orphanet",
            "product": "Product 6",
            "license": "CC-BY-4.0",
            "importer_version": IMPORTER_VERSION,
            "store_version": STORE_VERSION,
        }

        con.executemany(
            """
            INSERT INTO metadata(key, value)
            VALUES (?, ?)
            """,
            metadata.items(),
        )

        record_count = 0
        genes = set()
        diseases = set()

        for record in iter_orphanet_records(
            input_xml
        ):
            payload = record.to_dict()

            gene_symbol = (
                payload.get("gene_symbol")
                or ""
            ).upper()

            disease_id = payload.get(
                "disease_id"
            )

            con.execute(
                """
                INSERT INTO evidence(
                    gene_symbol,
                    disease_id,
                    disease_name,
                    record_json
                )
                VALUES (?, ?, ?, ?)
                """,
                (
                    gene_symbol,
                    disease_id,
                    payload.get("disease_name"),
                    json.dumps(
                        payload,
                        separators=(",", ":"),
                    ),
                ),
            )

            record_count += 1

            if gene_symbol:
                genes.add(gene_symbol)

            if disease_id:
                diseases.add(disease_id)

        con.commit()

        return {
            "record_count": record_count,
            "unique_genes": len(genes),
            "unique_diseases": len(diseases),
        }

    finally:
        con.close()


def query_gene(
    database_path: str,
    gene_symbol: str,
) -> List[Dict]:
    con = sqlite3.connect(database_path)

    try:
        rows = con.execute(
            """
            SELECT record_json
            FROM evidence
            WHERE gene_symbol = ?
            ORDER BY disease_id
            """,
            (gene_symbol.upper(),),
        ).fetchall()

        return [
            json.loads(row[0])
            for row in rows
        ]

    finally:
        con.close()


def query_disease(
    database_path: str,
    disease_id: str,
) -> List[Dict]:
    disease_id = disease_id.strip()

    if not disease_id.upper().startswith(
        "ORPHA:"
    ):
        disease_id = f"ORPHA:{disease_id}"

    con = sqlite3.connect(database_path)

    try:
        rows = con.execute(
            """
            SELECT record_json
            FROM evidence
            WHERE disease_id = ?
            ORDER BY gene_symbol
            """,
            (disease_id,),
        ).fetchall()

        return [
            json.loads(row[0])
            for row in rows
        ]

    finally:
        con.close()


def main():
    parser = argparse.ArgumentParser(
        description="Build/query Orphanet evidence store."
    )

    sub = parser.add_subparsers(
        dest="command",
        required=True,
    )

    build = sub.add_parser("build")
    build.add_argument(
        "--input-xml",
        required=True,
    )
    build.add_argument(
        "--database",
        required=True,
    )

    gene = sub.add_parser("gene")
    gene.add_argument(
        "--database",
        required=True,
    )
    gene.add_argument(
        "--gene",
        required=True,
    )

    disease = sub.add_parser("disease")
    disease.add_argument(
        "--database",
        required=True,
    )
    disease.add_argument(
        "--disease-id",
        required=True,
    )

    args = parser.parse_args()

    if args.command == "build":
        summary = build_database(
            args.input_xml,
            args.database,
        )

        print("Orphanet store")
        print(
            "  records:",
            summary["record_count"],
        )
        print(
            "  unique genes:",
            summary["unique_genes"],
        )
        print(
            "  unique diseases:",
            summary["unique_diseases"],
        )
        return

    if args.command == "gene":
        records = query_gene(
            args.database,
            args.gene,
        )
    else:
        records = query_disease(
            args.database,
            args.disease_id,
        )

    print(
        json.dumps(
            records,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
