#!/usr/bin/env python3

import argparse
import gzip
import json
import re
from collections import Counter
from pathlib import Path


SCHEMA_VERSION = "1.0"
PARSER_VERSION = "0.1.0"


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def parse_info(info_string):
    info = {}

    if info_string == ".":
        return info

    for item in info_string.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True

    return info


def extract_csq_fields(header_line):
    """
    Extract the dynamic VEP CSQ schema from a VCF header such as:

    ##INFO=<ID=CSQ,... Format: Allele|Consequence|IMPACT|SYMBOL|Gene|...">
    """
    match = re.search(r"Format:\s*([^\">]+)", header_line)

    if not match:
        raise ValueError(
            "Found CSQ header but could not determine its Format definition."
        )

    return [field.strip() for field in match.group(1).split("|")]


def normalize_csq_entry(raw_entry, csq_fields):
    values = raw_entry.split("|")

    if len(values) < len(csq_fields):
        values.extend([""] * (len(csq_fields) - len(values)))

    annotation = {
        field: values[index] if index < len(values) else ""
        for index, field in enumerate(csq_fields)
    }

    # Convenient normalized fields for downstream code.
    annotation["_normalized"] = {
        "allele": annotation.get("Allele", ""),
        "consequence": annotation.get("Consequence", ""),
        "impact": annotation.get("IMPACT", ""),
        "symbol": annotation.get("SYMBOL", ""),
        "gene_id": annotation.get("Gene", ""),
        "feature_type": annotation.get("Feature_type", ""),
        "transcript_id": annotation.get("Feature", ""),
        "biotype": annotation.get("BIOTYPE", ""),
        "protein_id": annotation.get("ENSP", ""),
        "canonical": annotation.get("CANONICAL", "") == "YES",
        "variant_class": annotation.get("VARIANT_CLASS", ""),
    }

    return annotation


def parse_vep_vcf(vcf_path):
    csq_fields = None
    variants = []

    consequence_counts = Counter()
    impact_counts = Counter()
    gene_counts = Counter()

    sample_names = []

    with open_text(vcf_path) as handle:
        for line in handle:
            line = line.rstrip("\n")

            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = extract_csq_fields(line)
                continue

            if line.startswith("#CHROM"):
                columns = line.split("\t")
                if len(columns) > 9:
                    sample_names = columns[9:]
                continue

            if line.startswith("#"):
                continue

            if not line:
                continue

            columns = line.split("\t")

            if len(columns) < 8:
                raise ValueError(
                    f"Malformed VCF record with fewer than 8 columns: {line}"
                )

            chrom, pos, variant_id, ref, alt, qual, filt, info_string = columns[:8]

            info = parse_info(info_string)

            transcripts = []

            if "CSQ" in info:
                if csq_fields is None:
                    raise ValueError(
                        "VCF contains CSQ annotations but no CSQ Format header."
                    )

                for raw_csq in info["CSQ"].split(","):
                    annotation = normalize_csq_entry(raw_csq, csq_fields)
                    transcripts.append(annotation)

                    normalized = annotation["_normalized"]

                    consequence = normalized["consequence"]
                    if consequence:
                        for term in consequence.split("&"):
                            consequence_counts[term] += 1

                    impact = normalized["impact"]
                    if impact:
                        impact_counts[impact] += 1

                    symbol = normalized["symbol"]
                    if symbol:
                        gene_counts[symbol] += 1

            variant = {
                "variant_key": f"{chrom}:{pos}:{ref}:{alt}",
                "chrom": chrom,
                "pos": int(pos),
                "id": None if variant_id == "." else variant_id,
                "ref": ref,
                "alt": alt.split(","),
                "qual": None if qual == "." else qual,
                "filter": filt.split(";") if filt != "." else [],
                "transcript_annotations": transcripts,
            }

            variants.append(variant)

    if csq_fields is None:
        raise ValueError(
            "No VEP CSQ header found. Input does not appear to be a VEP-annotated VCF."
        )

    canonical_transcripts = sum(
        1
        for variant in variants
        for transcript in variant["transcript_annotations"]
        if transcript["_normalized"]["canonical"]
    )

    unique_genes = sorted(
        {
            transcript["_normalized"]["symbol"]
            for variant in variants
            for transcript in variant["transcript_annotations"]
            if transcript["_normalized"]["symbol"]
        }
    )

    output = {
        "schema_version": SCHEMA_VERSION,
        "parser": {
            "name": "genomics_platform.vep_parser",
            "version": PARSER_VERSION,
        },
        "scope": {
            "type": "functional_variant_annotation",
            "clinical_classification": False,
            "candidate_ranking": False,
            "statement": (
                "VEP consequence annotation only. "
                "This output is not ACMG/AMP classification and does not "
                "establish pathogenicity."
            ),
        },
        "input": {
            "vcf": str(Path(vcf_path)),
            "samples": sample_names,
        },
        "vep": {
            "csq_fields": csq_fields,
        },
        "summary": {
            "variant_records": len(variants),
            "transcript_annotations": sum(
                len(v["transcript_annotations"]) for v in variants
            ),
            "canonical_transcripts": canonical_transcripts,
            "unique_genes": len(unique_genes),
            "genes": unique_genes,
            "impact_counts": dict(sorted(impact_counts.items())),
            "consequence_counts": dict(sorted(consequence_counts.items())),
            "gene_annotation_counts": dict(
                sorted(gene_counts.items(), key=lambda item: (-item[1], item[0]))
            ),
        },
        "variants": variants,
    }

    return output


def main():
    parser = argparse.ArgumentParser(
        description="Parse an Ensembl VEP annotated VCF into structured JSON."
    )

    parser.add_argument(
        "--input-vcf",
        required=True,
        help="VEP annotated VCF or VCF.GZ",
    )

    parser.add_argument(
        "--output-json",
        required=True,
        help="Output JSON path",
    )

    args = parser.parse_args()

    result = parse_vep_vcf(args.input_vcf)

    output_path = Path(args.output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=False)
        handle.write("\n")

    print("VEP annotation parser")
    print(f"  Input:                 {args.input_vcf}")
    print(f"  Output:                {args.output_json}")
    print(f"  Variant records:       {result['summary']['variant_records']}")
    print(
        f"  Transcript annotations:{result['summary']['transcript_annotations']:>8}"
    )
    print(f"  Unique genes:          {result['summary']['unique_genes']}")


if __name__ == "__main__":
    main()
