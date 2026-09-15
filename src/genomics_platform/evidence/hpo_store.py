#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
from pathlib import Path
from typing import Dict, List

from src.genomics_platform.evidence.evidence_record import EvidenceRecord


STORE_VERSION = "0.1.0"


SCHEMA = """
CREATE TABLE IF NOT EXISTS metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS gene_phenotype (
    id INTEGER PRIMARY KEY,
    ncbi_gene_id TEXT,
    gene_symbol TEXT NOT NULL,
    hpo_id TEXT NOT NULL,
    hpo_name TEXT,
    frequency TEXT,
    disease_id TEXT
);

CREATE INDEX IF NOT EXISTS idx_hpo_gene_symbol
ON gene_phenotype(gene_symbol);

CREATE INDEX IF NOT EXISTS idx_hpo_hpo_id
ON gene_phenotype(hpo_id);

CREATE INDEX IF NOT EXISTS idx_hpo_disease_id
ON gene_phenotype(disease_id);

CREATE INDEX IF NOT EXISTS idx_hpo_gene_disease
ON gene_phenotype(gene_symbol, disease_id);
"""


def build_database(input_tsv: str, database_path: str) -> Dict:
    database = Path(database_path)
    database.parent.mkdir(parents=True, exist_ok=True)

    if database.exists():
        database.unlink()

    con = sqlite3.connect(database)

    try:
        con.executescript(SCHEMA)

        con.execute(
            "INSERT INTO metadata(key, value) VALUES (?, ?)",
            ("store_version", STORE_VERSION),
        )

        con.execute(
            "INSERT INTO metadata(key, value) VALUES (?, ?)",
            ("source", "Human Phenotype Ontology"),
        )

        con.execute(
            "INSERT INTO metadata(key, value) VALUES (?, ?)",
            ("source_file", "genes_to_phenotype.txt"),
        )

        batch = []

        with open(
            input_tsv,
            "r",
            encoding="utf-8-sig",
            newline="",
        ) as handle:

            reader = csv.DictReader(
                handle,
                delimiter="\t",
            )

            required = {
                "ncbi_gene_id",
                "gene_symbol",
                "hpo_id",
                "hpo_name",
                "frequency",
                "disease_id",
            }

            headers = set(reader.fieldnames or [])

            missing = required - headers

            if missing:
                raise ValueError(
                    "Missing required HPO columns: "
                    f"{sorted(missing)}"
                )

            count = 0

            for row in reader:
                gene_symbol = (
                    row.get("gene_symbol") or ""
                ).strip()

                hpo_id = (
                    row.get("hpo_id") or ""
                ).strip()

                if not gene_symbol or not hpo_id:
                    continue

                batch.append(
                    (
                        (row.get("ncbi_gene_id") or "").strip() or None,
                        gene_symbol.upper(),
                        hpo_id,
                        (row.get("hpo_name") or "").strip() or None,
                        (row.get("frequency") or "").strip() or None,
                        (row.get("disease_id") or "").strip() or None,
                    )
                )

                if len(batch) >= 5000:
                    con.executemany(
                        """
                        INSERT INTO gene_phenotype (
                            ncbi_gene_id,
                            gene_symbol,
                            hpo_id,
                            hpo_name,
                            frequency,
                            disease_id
                        )
                        VALUES (?, ?, ?, ?, ?, ?)
                        """,
                        batch,
                    )

                    count += len(batch)
                    batch.clear()

            if batch:
                con.executemany(
                    """
                    INSERT INTO gene_phenotype (
                        ncbi_gene_id,
                        gene_symbol,
                        hpo_id,
                        hpo_name,
                        frequency,
                        disease_id
                    )
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    batch,
                )

                count += len(batch)

        con.commit()

        summary = con.execute(
            """
            SELECT
                COUNT(*) AS records,
                COUNT(DISTINCT gene_symbol) AS genes,
                COUNT(DISTINCT hpo_id) AS phenotypes,
                COUNT(DISTINCT disease_id) AS diseases
            FROM gene_phenotype
            """
        ).fetchone()

        return {
            "records": summary[0],
            "unique_genes": summary[1],
            "unique_phenotypes": summary[2],
            "unique_diseases": summary[3],
        }

    finally:
        con.close()


def row_to_record(row) -> EvidenceRecord:
    (
        ncbi_gene_id,
        gene_symbol,
        hpo_id,
        hpo_name,
        frequency,
        disease_id,
    ) = row

    gene_id = (
        f"NCBIGene:{ncbi_gene_id}"
        if ncbi_gene_id
        else None
    )

    description = hpo_id

    if hpo_name:
        description += f" {hpo_name}"

    return EvidenceRecord(
        source="HPO",
        evidence_type="gene_phenotype",
        gene_symbol=gene_symbol,
        gene_id=gene_id,
        disease_id=disease_id,
        disease_name=None,
        phenotype_ids=[hpo_id],
        inheritance=[],
        classification=None,
        confidence=None,
        frequency=frequency,
        source_record_id="|".join(
            x for x in (
                gene_id,
                disease_id,
                hpo_id,
            )
            if x
        ),
        citations=[],
        evidence_description=description,
        license=None,
    )


def query_gene(
    database_path: str,
    gene_symbol: str,
) -> List[EvidenceRecord]:

    con = sqlite3.connect(database_path)

    try:
        rows = con.execute(
            """
            SELECT
                ncbi_gene_id,
                gene_symbol,
                hpo_id,
                hpo_name,
                frequency,
                disease_id
            FROM gene_phenotype
            WHERE gene_symbol = ?
            ORDER BY disease_id, hpo_id
            """,
            (gene_symbol.upper(),),
        ).fetchall()

        return [
            row_to_record(row)
            for row in rows
        ]

    finally:
        con.close()


def main():
    parser = argparse.ArgumentParser(
        description="Build/query compact local HPO SQLite evidence store."
    )

    sub = parser.add_subparsers(
        dest="command",
        required=True,
    )

    build = sub.add_parser("build")

    build.add_argument(
        "--input",
        required=True,
    )

    build.add_argument(
        "--database",
        required=True,
    )

    query = sub.add_parser("query")

    query.add_argument(
        "--database",
        required=True,
    )

    query.add_argument(
        "--gene",
        required=True,
    )

    query.add_argument(
        "--json",
        action="store_true",
    )

    args = parser.parse_args()

    if args.command == "build":
        summary = build_database(
            args.input,
            args.database,
        )

        print("HPO SQLite store")
        print(f"  records:           {summary['records']}")
        print(f"  unique genes:      {summary['unique_genes']}")
        print(f"  unique phenotypes: {summary['unique_phenotypes']}")
        print(f"  unique diseases:   {summary['unique_diseases']}")

        return

    records = query_gene(
        args.database,
        args.gene,
    )

    if args.json:
        print(
            json.dumps(
                {
                    "gene": args.gene.upper(),
                    "status": (
                        "FOUND"
                        if records
                        else "NOT_FOUND"
                    ),
                    "record_count": len(records),
                    "records": [
                        record.to_dict()
                        for record in records
                    ],
                },
                indent=2,
            )
        )

        return

    print(f"Gene: {args.gene.upper()}")
    print(f"Records: {len(records)}")
    print()

    for record in records:
        phenotype = (
            record.evidence_description
            or "-"
        )

        print(
            f"- {record.disease_id} | "
            f"{phenotype} | "
            f"frequency={record.frequency or '-'}"
        )


if __name__ == "__main__":
    main()
