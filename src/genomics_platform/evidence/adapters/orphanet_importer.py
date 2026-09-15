#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from src.genomics_platform.evidence.evidence_record import EvidenceRecord


IMPORTER_VERSION = "0.1.0"
SOURCE_NAME = "Orphanet"
LICENSE = "CC-BY-4.0"


def clean(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None

    value = value.strip()
    return value if value else None


def child_text(
    element: Optional[ET.Element],
    path: str,
) -> Optional[str]:
    if element is None:
        return None

    child = element.find(path)

    if child is None:
        return None

    return clean(child.text)


def normalize_orpha_code(
    value: Optional[str],
) -> Optional[str]:
    value = clean(value)

    if not value:
        return None

    if value.upper().startswith("ORPHA:"):
        return f"ORPHA:{value.split(':', 1)[1]}"

    return f"ORPHA:{value}"


def extract_external_references(
    gene: ET.Element,
) -> Dict[str, List[str]]:
    references: Dict[str, List[str]] = {}

    for external in gene.findall(
        "./ExternalReferenceList/ExternalReference"
    ):
        source = child_text(external, "Source")
        reference = child_text(external, "Reference")

        if not source or not reference:
            continue

        references.setdefault(source, [])

        if reference not in references[source]:
            references[source].append(reference)

    return references


def extract_citations(
    value: Optional[str],
) -> List[str]:
    value = clean(value)

    if not value:
        return []

    citations: List[str] = []

    for pmid in re.findall(
        r"(\d+)\s*\[PMID\]",
        value,
        flags=re.IGNORECASE,
    ):
        citation = f"PMID:{pmid}"

        if citation not in citations:
            citations.append(citation)

    if citations:
        return citations

    return [value]


def preferred_gene_id(
    references: Dict[str, List[str]],
) -> Optional[str]:
    for source, prefix in (
        ("HGNC", "HGNC"),
        ("Ensembl", "ENSEMBL"),
    ):
        values = references.get(source) or []

        if values:
            return f"{prefix}:{values[0]}"

    return None


def build_description(
    association_type: Optional[str],
    association_status: Optional[str],
    gene_type: Optional[str],
    external_references: Dict[str, List[str]],
) -> Optional[str]:
    parts: List[str] = []

    if association_type:
        parts.append(
            f"association_type={association_type}"
        )

    if association_status:
        parts.append(
            f"association_status={association_status}"
        )

    if gene_type:
        parts.append(
            f"gene_type={gene_type}"
        )

    if external_references:
        refs = []

        for source in sorted(external_references):
            for value in external_references[source]:
                refs.append(f"{source}:{value}")

        if refs:
            parts.append(
                "external_references=" + ",".join(refs)
            )

    return "; ".join(parts) if parts else None


def association_to_record(
    disorder: ET.Element,
    association: ET.Element,
    source_release: Optional[str],
) -> Optional[EvidenceRecord]:
    orpha_code = normalize_orpha_code(
        child_text(disorder, "OrphaCode")
    )

    disease_name = child_text(
        disorder,
        "Name",
    )

    gene = association.find("Gene")

    if gene is None:
        return None

    gene_symbol = child_text(
        gene,
        "Symbol",
    )

    if not gene_symbol:
        return None

    gene_type = child_text(
        gene,
        "./GeneType/Name",
    )

    association_type = child_text(
        association,
        "./DisorderGeneAssociationType/Name",
    )

    association_status = child_text(
        association,
        "./DisorderGeneAssociationStatus/Name",
    )

    source_validation = child_text(
        association,
        "SourceOfValidation",
    )

    external_references = extract_external_references(
        gene
    )

    citations = extract_citations(
        source_validation
    )

    source_record_id = "|".join(
        value
        for value in (
            orpha_code,
            gene_symbol.upper(),
            association_type,
        )
        if value
    )

    return EvidenceRecord(
        source=SOURCE_NAME,
        evidence_type="gene_disease",
        gene_symbol=gene_symbol.upper(),
        gene_id=preferred_gene_id(
            external_references
        ),
        disease_id=orpha_code,
        disease_name=disease_name,
        phenotype_ids=[],
        inheritance=[],
        classification=None,
        confidence=association_status,
        frequency=None,
        source_record_id=source_record_id or None,
        source_version=IMPORTER_VERSION,
        source_release=source_release,
        citations=citations,
        evidence_description=build_description(
            association_type,
            association_status,
            gene_type,
            external_references,
        ),
        license=LICENSE,
    )


def iter_orphanet_records(
    path: str,
) -> Iterable[EvidenceRecord]:
    source_release: Optional[str] = None

    context = ET.iterparse(
        path,
        events=("start", "end"),
    )

    for event, elem in context:
        if event == "start" and elem.tag == "JDBOR":
            source_release = clean(
                elem.attrib.get("date")
            )
            continue

        if event != "end":
            continue

        if elem.tag != "Disorder":
            continue

        associations = elem.findall(
            "./DisorderGeneAssociationList/"
            "DisorderGeneAssociation"
        )

        for association in associations:
            record = association_to_record(
                elem,
                association,
                source_release,
            )

            if record is not None:
                yield record

        elem.clear()


def import_orphanet(
    path: str,
) -> Dict:
    records = [
        record.to_dict()
        for record in iter_orphanet_records(path)
    ]

    genes = {
        record["gene_symbol"]
        for record in records
        if record.get("gene_symbol")
    }

    diseases = {
        record["disease_id"]
        for record in records
        if record.get("disease_id")
    }

    statuses = Counter(
        record.get("confidence")
        for record in records
        if record.get("confidence")
    )

    return {
        "schema_version": "1.0",
        "importer": {
            "name": "orphanet_importer",
            "version": IMPORTER_VERSION,
        },
        "source": {
            "name": SOURCE_NAME,
            "product": "Product 6",
            "license": LICENSE,
        },
        "scope": {
            "clinical_classification": False,
            "variant_pathogenicity": False,
            "description": (
                "Orphanet gene-disease associations. "
                "Association status is source metadata "
                "and must not be interpreted as variant "
                "pathogenicity."
            ),
        },
        "summary": {
            "record_count": len(records),
            "unique_genes": len(genes),
            "unique_diseases": len(diseases),
            "association_statuses": dict(
                sorted(statuses.items())
            ),
        },
        "records": records,
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Import Orphanet Product 6 XML into "
            "normalized evidence JSON."
        )
    )

    parser.add_argument(
        "--input-xml",
        required=True,
    )

    parser.add_argument(
        "--output-json",
        required=True,
    )

    args = parser.parse_args()

    result = import_orphanet(
        args.input_xml
    )

    output = Path(args.output_json)
    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output.open(
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            result,
            handle,
            indent=2,
        )
        handle.write("\n")

    print("Orphanet importer")
    print(
        "  records:",
        result["summary"]["record_count"],
    )
    print(
        "  unique genes:",
        result["summary"]["unique_genes"],
    )
    print(
        "  unique diseases:",
        result["summary"]["unique_diseases"],
    )
    print(
        "  statuses:",
        result["summary"]["association_statuses"],
    )


if __name__ == "__main__":
    main()
