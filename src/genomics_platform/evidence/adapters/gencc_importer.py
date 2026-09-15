#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from src.genomics_platform.evidence.evidence_record import EvidenceRecord


IMPORTER_VERSION = "0.1.0"


def first_present(row: Dict[str, str], *keys: str) -> Optional[str]:
    for key in keys:
        value = row.get(key)
        if value is not None:
            value = value.strip()
            if value:
                return value
    return None


def split_multi(value: Optional[str]) -> List[str]:
    if not value:
        return []

    separators = [";", "|"]
    values = [value]

    for sep in separators:
        next_values = []
        for item in values:
            next_values.extend(item.split(sep))
        values = next_values

    return [item.strip() for item in values if item.strip()]


def normalize_gene_symbol(row: Dict[str, str]) -> Optional[str]:
    return first_present(
        row,
        "gene_symbol",
        "gene_symbol_text",
        "gene_label",
    )


def normalize_gene_id(row: Dict[str, str]) -> Optional[str]:
    return first_present(
        row,
        "gene_curie",
        "gene_id",
    )


def normalize_disease_id(row: Dict[str, str]) -> Optional[str]:
    return first_present(
        row,
        "disease_curie",
        "disease_id",
    )


def normalize_disease_name(row: Dict[str, str]) -> Optional[str]:
    return first_present(
        row,
        "disease_title",
        "disease_name",
        "disease_label",
    )


def normalize_inheritance(row: Dict[str, str]) -> List[str]:
    value = first_present(
        row,
        "moi_title",
        "mode_of_inheritance",
        "inheritance",
    )
    return split_multi(value)


def normalize_classification(row: Dict[str, str]) -> Optional[str]:
    return first_present(
        row,
        "classification_title",
        "classification",
    )


def normalize_citations(row: Dict[str, str]) -> List[str]:
    citations = []

    for key in (
        "pmids",
        "publication_pmids",
        "evidence_pmids",
    ):
        value = row.get(key)
        if value:
            citations.extend(split_multi(value))

    return sorted(set(citations))


def row_to_evidence_record(row: Dict[str, str]) -> EvidenceRecord:
    sgc_id = first_present(row, "sgc_id")
    version_number = first_present(row, "version_number")

    submitter = first_present(
        row,
        "submitter_title",
        "submitter",
    )

    classification = normalize_classification(row)

    description_parts = []

    if submitter:
        description_parts.append(f"Submitted by {submitter}")

    if classification:
        description_parts.append(
            f"Gene-disease validity classification: {classification}"
        )

    evidence_description = (
        ". ".join(description_parts)
        if description_parts
        else None
    )

    return EvidenceRecord(
        source="GenCC",
        evidence_type="gene_disease_validity",
        gene_symbol=normalize_gene_symbol(row),
        gene_id=normalize_gene_id(row),
        disease_id=normalize_disease_id(row),
        disease_name=normalize_disease_name(row),
        inheritance=normalize_inheritance(row),
        classification=classification,
        confidence=classification,
        source_record_id=sgc_id,
        source_version=version_number,
        source_release=first_present(
            row,
            "submitted_as_date",
            "submitted_date",
        ),
        citations=normalize_citations(row),
        evidence_description=evidence_description,
        license="CC0-1.0",
    )


def iter_gencc_csv(csv_path: str) -> Iterable[EvidenceRecord]:
    with open(csv_path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)

        if not reader.fieldnames:
            raise ValueError("GenCC CSV has no header.")

        if "sgc_id" not in reader.fieldnames:
            raise ValueError(
                "This does not appear to be the GenCC new-format export: "
                "required column 'sgc_id' is missing."
            )

        if "version_number" not in reader.fieldnames:
            raise ValueError(
                "This does not appear to be the GenCC new-format export: "
                "required column 'version_number' is missing."
            )

        for row in reader:
            yield row_to_evidence_record(row)


def import_gencc(csv_path: str):
    records = [record.to_dict() for record in iter_gencc_csv(csv_path)]

    genes = sorted(
        {
            record["gene_symbol"]
            for record in records
            if record["gene_symbol"]
        }
    )

    diseases = sorted(
        {
            record["disease_id"]
            for record in records
            if record["disease_id"]
        }
    )

    classifications = {}

    for record in records:
        classification = record.get("classification")
        if classification:
            classifications[classification] = (
                classifications.get(classification, 0) + 1
            )

    return {
        "schema_version": "1.0",
        "importer": {
            "name": "gencc_importer",
            "version": IMPORTER_VERSION,
        },
        "source": {
            "name": "GenCC",
            "license": "CC0-1.0",
            "format": "new",
        },
        "summary": {
            "record_count": len(records),
            "unique_genes": len(genes),
            "unique_diseases": len(diseases),
            "classifications": dict(sorted(classifications.items())),
        },
        "records": records,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Import GenCC new-format CSV into normalized evidence JSON."
    )

    parser.add_argument(
        "--input-csv",
        required=True,
    )

    parser.add_argument(
        "--output-json",
        required=True,
    )

    args = parser.parse_args()

    result = import_gencc(args.input_csv)

    output_path = Path(args.output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")

    print("GenCC importer")
    print(f"  records:         {result['summary']['record_count']}")
    print(f"  unique genes:    {result['summary']['unique_genes']}")
    print(f"  unique diseases: {result['summary']['unique_diseases']}")


if __name__ == "__main__":
    main()
