#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path


INDEX_SCHEMA_VERSION = "1.0"


def build_index(gencc_json_path: str):
    with open(gencc_json_path, encoding="utf-8") as handle:
        data = json.load(handle)

    by_gene = defaultdict(list)

    for record in data.get("records", []):
        gene_symbol = record.get("gene_symbol")

        if not gene_symbol:
            continue

        by_gene[gene_symbol.upper()].append(record)

    index = {
        "schema_version": INDEX_SCHEMA_VERSION,
        "source": "GenCC",
        "source_schema_version": data.get("schema_version"),
        "gene_count": len(by_gene),
        "by_gene": dict(sorted(by_gene.items())),
    }

    return index


def query_gene(index, gene_symbol: str):
    return index.get("by_gene", {}).get(gene_symbol.upper(), [])


def main():
    parser = argparse.ArgumentParser(
        description="Build or query a GenCC gene evidence index."
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    build_parser = subparsers.add_parser("build")

    build_parser.add_argument(
        "--input-json",
        required=True,
    )

    build_parser.add_argument(
        "--output-json",
        required=True,
    )

    query_parser = subparsers.add_parser("query")

    query_parser.add_argument(
        "--index-json",
        required=True,
    )

    query_parser.add_argument(
        "--gene",
        required=True,
    )

    args = parser.parse_args()

    if args.command == "build":
        index = build_index(args.input_json)

        output_path = Path(args.output_json)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with output_path.open("w", encoding="utf-8") as handle:
            json.dump(index, handle, indent=2)
            handle.write("\n")

        print("GenCC gene index")
        print(f"  indexed genes: {index['gene_count']}")

    elif args.command == "query":
        with open(args.index_json, encoding="utf-8") as handle:
            index = json.load(handle)

        records = query_gene(index, args.gene)

        print(f"Gene: {args.gene}")
        print(f"Records: {len(records)}")
        print()

        for record in records:
            print(
                f"- {record.get('disease_name')} | "
                f"{record.get('classification')} | "
                f"{', '.join(record.get('inheritance', []))}"
            )


if __name__ == "__main__":
    main()
