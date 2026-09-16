#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import urllib.parse
import urllib.request
from typing import Any, Dict, Iterable, List, Optional

from genomics_platform.evidence.evidence_record import EvidenceRecord


ADAPTER_VERSION = "0.1.0"
BASE_URL = "https://panelapp.genomicsengland.co.uk/api/v1/entities/"


CONFIDENCE_MAP = {
    "3": "green",
    "2": "amber",
    "1": "red",
    "0": "red",
    3: "green",
    2: "amber",
    1: "red",
    0: "red",
}


def first_present(obj: Dict[str, Any], *keys: str) -> Optional[Any]:
    for key in keys:
        if key in obj:
            value = obj[key]
            if value not in (None, "", [], {}):
                return value
    return None


def ensure_list(value: Any) -> List[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def unique_strings(values: Iterable[Any]) -> List[str]:
    output = []

    for value in values:
        if value is None:
            continue

        text = str(value).strip()

        if text and text not in output:
            output.append(text)

    return output


def extract_gene_symbol(entity: Dict[str, Any]) -> Optional[str]:
    gene_data = entity.get("gene_data") or {}

    return first_present(
        gene_data,
        "gene_symbol",
        "symbol",
    ) or first_present(
        entity,
        "entity_name",
        "gene_symbol",
        "name",
    )


def extract_hgnc_id(entity: Dict[str, Any]) -> Optional[str]:
    gene_data = entity.get("gene_data") or {}

    value = first_present(
        gene_data,
        "hgnc_id",
        "hgnc",
    )

    if not value:
        return None

    value = str(value)

    if value.upper().startswith("HGNC:"):
        return value

    return f"HGNC:{value}"


def extract_panel(entity: Dict[str, Any]) -> Dict[str, Any]:
    panel = entity.get("panel") or {}

    return {
        "id": first_present(
            panel,
            "id",
            "pk",
        ),
        "name": first_present(
            panel,
            "name",
            "title",
        ),
        "version": first_present(
            panel,
            "version",
            "version_number",
        ),
    }


def extract_confidence(entity: Dict[str, Any]) -> Optional[str]:
    level = first_present(
        entity,
        "confidence_level",
        "level_of_confidence",
        "gel_status",
    )

    if isinstance(level, dict):
        level = first_present(
            level,
            "value",
            "id",
            "name",
            "label",
        )

    if level is None:
        return None

    text = str(level).strip()

    mapped = CONFIDENCE_MAP.get(level)
    if mapped:
        return mapped

    mapped = CONFIDENCE_MAP.get(text)
    if mapped:
        return mapped

    lower = text.lower()

    if "green" in lower or "high" in lower:
        return "green"

    if "amber" in lower or "moderate" in lower:
        return "amber"

    if "red" in lower or "low" in lower:
        return "red"

    return text


def extract_inheritance(entity: Dict[str, Any]) -> List[str]:
    value = first_present(
        entity,
        "mode_of_inheritance",
        "moi",
    )

    if isinstance(value, str):
        return [value.strip()] if value.strip() else []

    return unique_strings(ensure_list(value))


def extract_phenotypes(entity: Dict[str, Any]) -> List[str]:
    phenotypes = first_present(
        entity,
        "phenotypes",
        "phenotype",
    )

    output = []

    for phenotype in ensure_list(phenotypes):
        if isinstance(phenotype, dict):
            value = first_present(
                phenotype,
                "name",
                "label",
                "phenotype",
                "id",
            )
        else:
            value = phenotype

        if value:
            output.append(str(value).strip())

    return unique_strings(output)


def extract_publications(entity: Dict[str, Any]) -> List[str]:
    publications = first_present(
        entity,
        "publications",
        "publication",
    )

    output = []

    for publication in ensure_list(publications):
        if isinstance(publication, dict):
            value = first_present(
                publication,
                "pmid",
                "id",
                "title",
            )
        else:
            value = publication

        if value:
            output.append(str(value).strip())

    return unique_strings(output)


def panelapp_entity_to_record(
    entity: Dict[str, Any],
) -> EvidenceRecord:

    panel = extract_panel(entity)
    gene_symbol = extract_gene_symbol(entity)

    source_id_parts = []

    if panel["id"] is not None:
        source_id_parts.append(f"panel:{panel['id']}")

    if gene_symbol:
        source_id_parts.append(f"gene:{gene_symbol}")

    source_record_id = (
        "|".join(source_id_parts)
        if source_id_parts
        else None
    )

    description_parts = []

    if panel["name"]:
        description_parts.append(
            f"PanelApp panel: {panel['name']}"
        )

    phenotypes = extract_phenotypes(entity)

    if phenotypes:
        description_parts.append(
            "Phenotypes: " + "; ".join(phenotypes)
        )

    return EvidenceRecord(
        source="PanelApp",
        evidence_type="expert_gene_panel",
        gene_symbol=gene_symbol,
        gene_id=extract_hgnc_id(entity),
        disease_id=None,
        disease_name="; ".join(phenotypes) if phenotypes else None,
        inheritance=extract_inheritance(entity),
        classification=None,
        confidence=extract_confidence(entity),
        source_record_id=source_record_id,
        source_version=(
            str(panel["version"])
            if panel["version"] is not None
            else None
        ),
        source_release=None,
        citations=extract_publications(entity),
        evidence_description=(
            ". ".join(description_parts)
            if description_parts
            else None
        ),
        license=None,
    )


def parse_api_payload(
    payload: Dict[str, Any],
) -> List[EvidenceRecord]:

    results = payload.get("results")

    if results is None:
        if isinstance(payload, list):
            results = payload
        else:
            results = [payload]

    records = []

    for entity in results:
        if not isinstance(entity, dict):
            continue

        entity_type = str(
            entity.get("entity_type", "gene")
        ).lower()

        if entity_type not in ("gene", ""):
            continue

        records.append(
            panelapp_entity_to_record(entity)
        )

    return records


def fetch_gene(
    gene_symbol: str,
) -> Dict[str, Any]:

    params = urllib.parse.urlencode(
        {
            "entity_name": gene_symbol,
            "format": "json",
        }
    )

    url = f"{BASE_URL}?{params}"

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent":
                "modular-human-genomics-pipeline/0.1"
        },
    )

    with urllib.request.urlopen(
        request,
        timeout=60,
    ) as response:
        return json.load(response)


def build_output(
    gene_symbol: str,
    records: List[EvidenceRecord],
) -> Dict[str, Any]:

    normalized_records = [
        record.to_dict()
        for record in records
        if (
            record.gene_symbol
            and record.gene_symbol.upper()
            == gene_symbol.upper()
        )
    ]

    confidence_counts = {}

    for record in normalized_records:
        confidence = record.get("confidence")

        if confidence:
            confidence_counts[confidence] = (
                confidence_counts.get(
                    confidence,
                    0,
                ) + 1
            )

    return {
        "schema_version": "1.0",
        "adapter": {
            "name": "panelapp_adapter",
            "version": ADAPTER_VERSION,
        },
        "source": {
            "name": "PanelApp",
            "endpoint": BASE_URL,
        },
        "query": {
            "gene_symbol": gene_symbol.upper(),
        },
        "status": (
            "FOUND"
            if normalized_records
            else "NOT_FOUND"
        ),
        "summary": {
            "record_count":
                len(normalized_records),
            "confidence_counts":
                dict(
                    sorted(
                        confidence_counts.items()
                    )
                ),
        },
        "records": normalized_records,
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Query PanelApp gene evidence "
            "and normalize it."
        )
    )

    parser.add_argument(
        "--gene",
        required=True,
    )

    parser.add_argument(
        "--output-json",
        required=True,
    )

    parser.add_argument(
        "--input-json",
        help=(
            "Optional local PanelApp API JSON "
            "for offline/regression testing."
        ),
    )

    args = parser.parse_args()

    if args.input_json:
        with open(
            args.input_json,
            encoding="utf-8",
        ) as handle:
            payload = json.load(handle)
    else:
        payload = fetch_gene(args.gene)

    records = parse_api_payload(payload)

    result = build_output(
        args.gene,
        records,
    )

    with open(
        args.output_json,
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            result,
            handle,
            indent=2,
        )
        handle.write("\n")

    print("PanelApp adapter")
    print(
        f"  gene:       "
        f"{result['query']['gene_symbol']}"
    )
    print(
        f"  status:     "
        f"{result['status']}"
    )
    print(
        f"  records:    "
        f"{result['summary']['record_count']}"
    )
    print(
        f"  confidence: "
        f"{result['summary']['confidence_counts']}"
    )


if __name__ == "__main__":
    main()
