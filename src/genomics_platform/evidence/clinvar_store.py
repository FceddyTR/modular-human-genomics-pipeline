#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sqlite3
import time
from pathlib import Path
from typing import Dict, List, Optional

from genomics_platform.evidence.adapters.clinvar_vcf_parser import (
    PARSER_VERSION,
    iter_clinvar_records,
)


STORE_VERSION = "0.1.0"


def canonical_chrom(chrom: str) -> str:
    chrom = chrom.strip()

    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]

    if chrom == "M":
        return "MT"

    return chrom


def build_database(
    input_vcf: str,
    database_path: str,
) -> Dict:
    database = Path(database_path)

    database.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if database.exists():
        database.unlink()

    con = sqlite3.connect(database)

    try:
        # Faster bulk import. Database is reproducible from source VCF.
        con.execute("PRAGMA journal_mode=OFF")
        con.execute("PRAGMA synchronous=OFF")
        con.execute("PRAGMA temp_store=MEMORY")
        con.execute("PRAGMA cache_size=-131072")

        con.executescript(
            """
            CREATE TABLE metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );

            CREATE TABLE variant (
                id INTEGER PRIMARY KEY,
                assembly TEXT NOT NULL,
                chrom TEXT NOT NULL,
                pos INTEGER NOT NULL,
                ref TEXT NOT NULL,
                alt TEXT NOT NULL,
                variation_id TEXT,
                allele_id TEXT,
                clinical_significance TEXT,
                review_stars INTEGER,
                record_json TEXT NOT NULL
            );
            """
        )

        metadata = {
            "source": "ClinVar",
            "assembly": "GRCh38",
            "parser_version": PARSER_VERSION,
            "store_version": STORE_VERSION,
        }

        con.executemany(
            """
            INSERT INTO metadata(key, value)
            VALUES (?, ?)
            """,
            metadata.items(),
        )

        started = time.time()
        count = 0
        source_release: Optional[str] = None

        batch = []

        for record in iter_clinvar_records(input_vcf):
            payload = record.to_dict()

            if source_release is None:
                source_release = record.source_release

            chrom = canonical_chrom(record.chrom)

            batch.append(
                (
                    record.assembly,
                    chrom,
                    record.pos,
                    record.ref,
                    record.alt,
                    record.variation_id,
                    record.allele_id,
                    "|".join(
                        record.clinical_significance
                    ),
                    record.review_stars,
                    json.dumps(
                        payload,
                        separators=(",", ":"),
                    ),
                )
            )

            if len(batch) >= 10000:
                con.executemany(
                    """
                    INSERT INTO variant(
                        assembly,
                        chrom,
                        pos,
                        ref,
                        alt,
                        variation_id,
                        allele_id,
                        clinical_significance,
                        review_stars,
                        record_json
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    batch,
                )

                count += len(batch)
                batch.clear()

                if count % 500000 == 0:
                    elapsed = time.time() - started
                    print(
                        f"  imported {count:,} records "
                        f"({elapsed:.1f}s)"
                    )

        if batch:
            con.executemany(
                """
                INSERT INTO variant(
                    assembly,
                    chrom,
                    pos,
                    ref,
                    alt,
                    variation_id,
                    allele_id,
                    clinical_significance,
                    review_stars,
                    record_json
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                batch,
            )

            count += len(batch)

        print("  creating coordinate index...")

        con.execute(
            """
            CREATE INDEX idx_clinvar_variant
            ON variant(
                assembly,
                chrom,
                pos,
                ref,
                alt
            )
            """
        )

        con.execute(
            """
            CREATE INDEX idx_clinvar_variation_id
            ON variant(variation_id)
            """
        )

        if source_release:
            con.execute(
                """
                INSERT OR REPLACE INTO metadata(key, value)
                VALUES ('source_release', ?)
                """,
                (source_release,),
            )

        con.execute(
            """
            INSERT OR REPLACE INTO metadata(key, value)
            VALUES ('record_count', ?)
            """,
            (str(count),),
        )

        con.commit()

        elapsed = time.time() - started

        return {
            "record_count": count,
            "source_release": source_release,
            "elapsed_seconds": round(elapsed, 2),
        }

    finally:
        con.close()


def query_variant(
    database_path: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    assembly: str = "GRCh38",
) -> Dict:
    chrom = canonical_chrom(chrom)

    con = sqlite3.connect(database_path)

    try:
        rows = con.execute(
            """
            SELECT record_json
            FROM variant
            WHERE assembly = ?
              AND chrom = ?
              AND pos = ?
              AND ref = ?
              AND alt = ?
            ORDER BY variation_id
            """,
            (
                assembly,
                chrom,
                int(pos),
                ref,
                alt,
            ),
        ).fetchall()

    finally:
        con.close()

    records = [
        json.loads(row[0])
        for row in rows
    ]

    if not records:
        return {
            "status": "NOT_FOUND",
            "interpretation": "NO_CLINVAR_RECORD",
            "query": {
                "assembly": assembly,
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
            },
            "records": [],
        }

    return {
        "status": "FOUND",
        "interpretation": "CLINVAR_RECORD_PRESENT",
        "query": {
            "assembly": assembly,
            "chrom": chrom,
            "pos": int(pos),
            "ref": ref,
            "alt": alt,
        },
        "records": records,
    }


def query_variation_id(
    database_path: str,
    variation_id: str,
) -> List[Dict]:
    con = sqlite3.connect(database_path)

    try:
        rows = con.execute(
            """
            SELECT record_json
            FROM variant
            WHERE variation_id = ?
            """,
            (str(variation_id),),
        ).fetchall()

    finally:
        con.close()

    return [
        json.loads(row[0])
        for row in rows
    ]


def main():
    parser = argparse.ArgumentParser(
        description="Build/query local ClinVar GRCh38 store."
    )

    sub = parser.add_subparsers(
        dest="command",
        required=True,
    )

    build = sub.add_parser("build")

    build.add_argument(
        "--input-vcf",
        required=True,
    )

    build.add_argument(
        "--database",
        required=True,
    )

    query = sub.add_parser("variant")

    query.add_argument(
        "--database",
        required=True,
    )
    query.add_argument(
        "--chrom",
        required=True,
    )
    query.add_argument(
        "--pos",
        required=True,
        type=int,
    )
    query.add_argument(
        "--ref",
        required=True,
    )
    query.add_argument(
        "--alt",
        required=True,
    )

    vid = sub.add_parser("variation-id")

    vid.add_argument(
        "--database",
        required=True,
    )
    vid.add_argument(
        "--id",
        required=True,
    )

    args = parser.parse_args()

    if args.command == "build":
        result = build_database(
            args.input_vcf,
            args.database,
        )

        print("ClinVar store")
        print(
            "  records:",
            f'{result["record_count"]:,}',
        )
        print(
            "  release:",
            result["source_release"],
        )
        print(
            "  elapsed:",
            result["elapsed_seconds"],
            "seconds",
        )

        return

    if args.command == "variant":
        result = query_variant(
            args.database,
            args.chrom,
            args.pos,
            args.ref,
            args.alt,
        )
    else:
        result = query_variation_id(
            args.database,
            args.id,
        )

    print(
        json.dumps(
            result,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
