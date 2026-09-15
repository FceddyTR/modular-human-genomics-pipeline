#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Dict, List

from src.genomics_platform.evidence.disease_normalizer import (
    normalize_disease,
)
from src.genomics_platform.evidence.hpo_store import (
    query_gene as query_hpo_gene,
)
from src.genomics_platform.evidence.orphanet_store import (
    query_gene as query_orphanet_gene,
)


ENGINE_VERSION = "0.2.0"
SCHEMA_VERSION = "1.0"


def load_json(path: str) -> Dict:
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def load_gencc_records(
    index_path: str,
    gene_symbol: str,
) -> List[Dict]:
    data = load_json(index_path)

    return list(
        data.get("by_gene", {}).get(
            gene_symbol.upper(),
            [],
        )
    )


def load_panelapp_records(
    panelapp_path: str,
    gene_symbol: str,
) -> List[Dict]:
    path = Path(panelapp_path)

    if not path.exists():
        return []

    data = load_json(panelapp_path)

    records = []

    for record in data.get("records", []):
        symbol = (
            record.get("gene_symbol") or ""
        ).upper()

        if symbol == gene_symbol.upper():
            records.append(record)

    return records


def load_hpo_records(
    database_path: str,
    gene_symbol: str,
) -> List[Dict]:
    return [
        record.to_dict()
        for record in query_hpo_gene(
            database_path,
            gene_symbol,
        )
    ]


def load_orphanet_records(
    database_path: str,
    gene_symbol: str,
) -> List[Dict]:
    path = Path(database_path)

    if not path.exists():
        return []

    return query_orphanet_gene(
        database_path,
        gene_symbol,
    )


def enrich_record(
    record: Dict,
    mondo_database: str,
) -> Dict:
    disease_id = record.get("disease_id")

    normalization = normalize_disease(
        mondo_database,
        disease_id,
    )

    return {
        "evidence": record,
        "disease_normalization": normalization,
    }


def build_bundle(
    gene_symbol: str,
    gencc_index: str,
    panelapp_json: str,
    hpo_database: str,
    mondo_database: str,
    orphanet_database: str,
) -> Dict:
    gene_symbol = gene_symbol.upper()

    gencc = load_gencc_records(
        gencc_index,
        gene_symbol,
    )

    panelapp = load_panelapp_records(
        panelapp_json,
        gene_symbol,
    )

    hpo = load_hpo_records(
        hpo_database,
        gene_symbol,
    )

    orphanet = load_orphanet_records(
        orphanet_database,
        gene_symbol,
    )

    source_records = (
        gencc
        + panelapp
        + hpo
        + orphanet
    )

    enriched = [
        enrich_record(
            record,
            mondo_database,
        )
        for record in source_records
    ]

    source_counts = Counter(
        item["evidence"].get(
            "source",
            "UNKNOWN",
        )
        for item in enriched
    )

    normalization_counts = Counter(
        item["disease_normalization"]["status"]
        for item in enriched
    )

    mondo_diseases = {}

    for item in enriched:
        normalization = item[
            "disease_normalization"
        ]

        mondo_id = normalization.get(
            "mondo_id"
        )

        if not mondo_id:
            continue

        mondo_diseases[mondo_id] = (
            normalization.get(
                "mondo_label"
            )
        )

    return {
        "schema_version": SCHEMA_VERSION,
        "engine": {
            "name": "gene_evidence_bundle",
            "version": ENGINE_VERSION,
        },
        "query": {
            "gene_symbol": gene_symbol,
        },
        "scope": {
            "clinical_classification": False,
            "variant_pathogenicity": False,
            "candidate_ranking": False,
            "description": (
                "Evidence aggregation and disease identifier "
                "normalization only."
            ),
        },
        "summary": {
            "record_count": len(enriched),
            "source_counts": dict(
                sorted(source_counts.items())
            ),
            "disease_normalization": dict(
                sorted(
                    normalization_counts.items()
                )
            ),
            "normalized_disease_count": len(
                mondo_diseases
            ),
        },
        "normalized_diseases": [
            {
                "mondo_id": mondo_id,
                "mondo_label": mondo_diseases[
                    mondo_id
                ],
            }
            for mondo_id in sorted(
                mondo_diseases
            )
        ],
        "records": enriched,
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate gene-level evidence from "
            "GenCC, PanelApp, HPO and Orphanet "
            "and normalize disease identifiers "
            "with MONDO."
        )
    )

    parser.add_argument(
        "--gene",
        required=True,
    )

    parser.add_argument(
        "--gencc-index",
        required=True,
    )

    parser.add_argument(
        "--panelapp-json",
        required=True,
    )

    parser.add_argument(
        "--hpo-database",
        required=True,
    )

    parser.add_argument(
        "--mondo-database",
        required=True,
    )

    parser.add_argument(
        "--orphanet-database",
        required=True,
    )

    parser.add_argument(
        "--output-json",
        required=True,
    )

    args = parser.parse_args()

    bundle = build_bundle(
        gene_symbol=args.gene,
        gencc_index=args.gencc_index,
        panelapp_json=args.panelapp_json,
        hpo_database=args.hpo_database,
        mondo_database=args.mondo_database,
        orphanet_database=args.orphanet_database,
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
            bundle,
            handle,
            indent=2,
        )
        handle.write("\n")

    summary = bundle["summary"]

    print("Gene Evidence Bundle")
    print(
        "  gene:",
        bundle["query"]["gene_symbol"],
    )
    print(
        "  records:",
        summary["record_count"],
    )
    print(
        "  sources:",
        summary["source_counts"],
    )
    print(
        "  normalization:",
        summary["disease_normalization"],
    )
    print(
        "  normalized diseases:",
        summary["normalized_disease_count"],
    )


if __name__ == "__main__":
    main()
