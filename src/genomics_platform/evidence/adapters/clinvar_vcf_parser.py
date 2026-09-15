#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import json
from typing import Dict, Iterable, List, Optional

from src.genomics_platform.evidence.clinvar_record import (
    ClinVarCondition,
    ClinVarRecord,
    MolecularConsequence,
)


PARSER_VERSION = "0.1.0"


def split_pipe(value: Optional[str]) -> List[str]:
    if not value or value == ".":
        return []

    return [
        item
        for item in value.split("|")
        if item
    ]


def parse_info(info_string: str) -> Dict[str, str]:
    result: Dict[str, str] = {}

    for item in info_string.split(";"):
        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        else:
            result[item] = "true"

    return result


def parse_geneinfo(
    value: Optional[str],
) -> tuple[List[str], List[str]]:
    symbols: List[str] = []
    ids: List[str] = []

    if not value:
        return symbols, ids

    for item in value.split("|"):
        if ":" not in item:
            continue

        symbol, gene_id = item.rsplit(":", 1)

        if symbol and symbol not in symbols:
            symbols.append(symbol)

        if gene_id and gene_id not in ids:
            ids.append(gene_id)

    return symbols, ids


def review_stars(
    review_status: Optional[str],
) -> Optional[int]:
    if not review_status:
        return None

    status = review_status.lower()

    if "practice_guideline" in status:
        return 4

    if "reviewed_by_expert_panel" in status:
        return 3

    if (
        "criteria_provided,_multiple_submitters,"
        "_no_conflicts" in status
    ):
        return 2

    if "criteria_provided" in status:
        return 1

    if (
        "no_assertion_criteria_provided" in status
        or "no_assertion_provided" in status
    ):
        return 0

    return None


def parse_conditions(
    names_value: Optional[str],
    ids_value: Optional[str],
) -> List[ClinVarCondition]:
    names = split_pipe(names_value)
    identifier_groups = split_pipe(ids_value)

    count = max(
        len(names),
        len(identifier_groups),
    )

    conditions: List[ClinVarCondition] = []

    for index in range(count):
        name = (
            names[index]
            if index < len(names)
            else None
        )

        identifiers: List[str] = []

        if index < len(identifier_groups):
            identifiers = [
                item
                for item in identifier_groups[
                    index
                ].split(",")
                if item
            ]

        conditions.append(
            ClinVarCondition(
                name=name,
                identifiers=identifiers,
            )
        )

    return conditions



def parse_molecular_consequences(
    value: Optional[str],
) -> List[MolecularConsequence]:
    if not value or value == ".":
        return []

    results: List[MolecularConsequence] = []

    # ClinVar MC is comma-separated. Each entry is:
    # SequenceOntologyID|molecular_consequence
    for item in value.split(","):
        item = item.strip()

        if not item:
            continue

        if "|" in item:
            so_id, consequence = item.split("|", 1)
        else:
            so_id = item
            consequence = None

        results.append(
            MolecularConsequence(
                so_id=so_id or None,
                consequence=consequence or None,
            )
        )

    return results

def variation_to_vcv(
    variation_id: Optional[str],
) -> Optional[str]:
    if not variation_id:
        return None

    try:
        number = int(variation_id)
    except ValueError:
        return None

    return f"VCV{number:09d}"


def parse_record(
    line: str,
    source_release: Optional[str],
    assembly: str = "GRCh38",
) -> List[ClinVarRecord]:
    fields = line.rstrip("\n").split("\t")

    if len(fields) < 8:
        raise ValueError(
            "VCF record has fewer than 8 columns"
        )

    chrom, pos, variation_id, ref, alts = (
        fields[0],
        fields[1],
        fields[2],
        fields[3],
        fields[4],
    )

    info = parse_info(fields[7])

    alt_list = alts.split(",")

    if len(alt_list) != 1:
        raise ValueError(
            "ClinVar parser v0.1 expects "
            "biallelic normalized records"
        )

    gene_symbols, gene_ids = parse_geneinfo(
        info.get("GENEINFO")
    )

    record = ClinVarRecord(
        assembly=assembly,
        chrom=chrom,
        pos=int(pos),
        ref=ref,
        alt=alt_list[0],
        variation_id=(
            variation_id
            if variation_id != "."
            else None
        ),
        allele_id=info.get("ALLELEID"),
        vcv_accession=variation_to_vcv(
            variation_id
            if variation_id != "."
            else None
        ),
        gene_symbols=gene_symbols,
        gene_ids=gene_ids,
        clinical_significance=split_pipe(
            info.get("CLNSIG")
        ),
        conflicting_significance=split_pipe(
            info.get("CLNSIGCONF")
        ),
        review_status=info.get("CLNREVSTAT"),
        review_stars=review_stars(
            info.get("CLNREVSTAT")
        ),
        conditions=parse_conditions(
            info.get("CLNDN"),
            info.get("CLNDISDB"),
        ),
        hgvs=split_pipe(
            info.get("CLNHGVS")
        ),
        molecular_consequences=parse_molecular_consequences(
            info.get("MC")
        ),
        scv_accessions=split_pipe(
            info.get("CLNSIGSCV")
        ),
        variant_type=info.get("CLNVC"),
        variant_type_so=info.get("CLNVCSO"),
        origin=info.get("ORIGIN"),
        source_release=source_release,
        source_reference=assembly,
    )

    return [record]


def iter_clinvar_records(
    path: str,
) -> Iterable[ClinVarRecord]:
    source_release: Optional[str] = None
    assembly = "GRCh38"

    opener = (
        gzip.open
        if path.endswith(".gz")
        else open
    )

    with opener(
        path,
        "rt",
        encoding="utf-8",
    ) as handle:
        for line in handle:
            if line.startswith("##fileDate="):
                source_release = (
                    line.strip().split("=", 1)[1]
                )
                continue

            if line.startswith("##reference="):
                assembly = (
                    line.strip().split("=", 1)[1]
                )
                continue

            if line.startswith("#"):
                continue

            for record in parse_record(
                line,
                source_release=source_release,
                assembly=assembly,
            ):
                yield record


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input-vcf",
        required=True,
    )

    parser.add_argument(
        "--limit",
        type=int,
        default=5,
    )

    args = parser.parse_args()

    for index, record in enumerate(
        iter_clinvar_records(args.input_vcf)
    ):
        print(
            json.dumps(
                record.to_dict(),
                indent=2,
            )
        )

        if index + 1 >= args.limit:
            break


if __name__ == "__main__":
    main()
