#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from genomics_platform.evidence.evidence_record import EvidenceRecord


IMPORTER_VERSION = "0.1.0"


def clean(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    value = value.strip()
    return value if value else None


def normalize_gene_id(value: Optional[str]) -> Optional[str]:
    value = clean(value)

    if not value:
        return None

    if value.startswith("NCBIGene:"):
        return value

    return f"NCBIGene:{value}"


def build_description(
    hpo_id: Optional[str],
    hpo_name: Optional[str],
    frequency: Optional[str],
) -> Optional[str]:

    parts = []

    if hpo_id or hpo_name:
        if hpo_id and hpo_name:
            parts.append(f"{hpo_id} {hpo_name}")
        else:
            parts.append(hpo_id or hpo_name)

    if frequency and frequency != "-":
        parts.append(f"frequency={frequency}")

    return "; ".join(parts) if parts else None


def row_to_record(row: Dict[str, str]) -> EvidenceRecord:
    gene_id = normalize_gene_id(row.get("ncbi_gene_id"))
    gene_symbol = clean(row.get("gene_symbol"))

    hpo_id = clean(row.get("hpo_id"))
    hpo_name = clean(row.get("hpo_name"))
    frequency = clean(row.get("frequency"))
    disease_id = clean(row.get("disease_id"))

    source_record_id_parts = [
        value
        for value in (
            gene_id,
            disease_id,
            hpo_id,
        )
        if value
    ]

    return EvidenceRecord(
        source="HPO",
        evidence_type="gene_phenotype",
        gene_symbol=gene_symbol,
        gene_id=gene_id,
        disease_id=disease_id,
        disease_name=None,
        phenotype_ids=[hpo_id] if hpo_id else [],
        inheritance=[],
        classification=None,
        confidence=None,
        frequency=frequency,
        source_record_id="|".join(source_record_id_parts)
        if source_record_id_parts
        else None,
        source_version=None,
        source_release=None,
        citations=[],
        evidence_description=build_description(
            hpo_id,
            hpo_name,
            frequency,
        ),
        license=None,
    )


def iter_hpo_records(path: str) -> Iterable[EvidenceRecord]:
    with open(
        path,
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
                "HPO genes_to_phenotype file is missing required "
                f"columns: {sorted(missing)}"
            )

        for row in reader:
            yield row_to_record(row)


def import_hpo(path: str) -> Dict:

    records = [
        record.to_dict()
        for record in iter_hpo_records(path)
    ]

    genes = {
        record["gene_symbol"]
        for record in records
        if record.get("gene_symbol")
    }

    phenotypes = {
        phenotype
        for record in records
        for phenotype in record.get(
            "phenotype_ids",
            [],
        )
    }

    diseases = {
        record["disease_id"]
        for record in records
        if record.get("disease_id")
    }

    disease_prefixes = Counter()

    for disease_id in diseases:
        if ":" in disease_id:
            prefix = disease_id.split(":", 1)[0]
        else:
            prefix = disease_id

        disease_prefixes[prefix] += 1

    return {
        "schema_version": "1.0",
        "importer": {
            "name": "hpo_importer",
            "version": IMPORTER_VERSION,
        },
        "source": {
            "name": "Human Phenotype Ontology",
            "file": "genes_to_phenotype.txt",
        },
        "summary": {
            "record_count": len(records),
            "unique_genes": len(genes),
            "unique_phenotypes": len(phenotypes),
            "unique_diseases": len(diseases),
            "disease_prefixes": dict(
                sorted(
                    disease_prefixes.items()
                )
            ),
        },
        "records": records,
    }


def build_gene_index(data: Dict) -> Dict:
    by_gene = defaultdict(list)

    for record in data.get("records", []):
        gene_symbol = record.get("gene_symbol")

        if not gene_symbol:
            continue

        by_gene[
            gene_symbol.upper()
        ].append(record)

    return {
        "schema_version": "1.0",
        "source": "HPO",
        "index_type": "gene_to_phenotype",
        "gene_count": len(by_gene),
        "by_gene": dict(
            sorted(by_gene.items())
        ),
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Import HPO genes_to_phenotype.txt "
            "and optionally build a gene index."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
    )

    parser.add_argument(
        "--output-json",
        required=True,
    )

    parser.add_argument(
        "--index-json",
    )

    args = parser.parse_args()

    result = import_hpo(args.input)

    output_path = Path(args.output_json)
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            result,
            handle,
            indent=2,
        )
        handle.write("\n")

    print("HPO importer")
    print(
        f"  records:           "
        f"{result['summary']['record_count']}"
    )
    print(
        f"  unique genes:      "
        f"{result['summary']['unique_genes']}"
    )
    print(
        f"  unique phenotypes: "
        f"{result['summary']['unique_phenotypes']}"
    )
    print(
        f"  unique diseases:   "
        f"{result['summary']['unique_diseases']}"
    )

    if args.index_json:
        index = build_gene_index(result)

        index_path = Path(args.index_json)
        index_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with index_path.open(
            "w",
            encoding="utf-8",
        ) as handle:
            json.dump(
                index,
                handle,
                indent=2,
            )
            handle.write("\n")

        print(
            f"  indexed genes:     "
            f"{index['gene_count']}"
        )


if __name__ == "__main__":
    main()
