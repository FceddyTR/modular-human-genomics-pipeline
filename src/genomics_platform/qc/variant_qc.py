#!/usr/bin/env python3

import argparse
import gzip
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


TRANSITIONS = {
    ("A", "G"),
    ("G", "A"),
    ("C", "T"),
    ("T", "C"),
}


def open_text(path):
    path = str(path)

    if path.endswith(".gz"):
        return gzip.open(path, "rt")

    return open(path, "r", encoding="utf-8")


def safe_float(value):
    if value in (None, "", "."):
        return None

    try:
        return float(value)
    except ValueError:
        return None


def safe_int(value):
    if value in (None, "", "."):
        return None

    try:
        return int(value)
    except ValueError:
        return None


def percentile(values, pct):
    if not values:
        return None

    values = sorted(values)

    if len(values) == 1:
        return values[0]

    rank = (len(values) - 1) * pct

    lo = math.floor(rank)
    hi = math.ceil(rank)

    if lo == hi:
        return values[lo]

    weight = rank - lo

    return (
        values[lo] * (1 - weight)
        + values[hi] * weight
    )


def variant_type(ref, alt):
    if alt.startswith("<"):
        return "SYMBOLIC"

    if len(ref) == 1 and len(alt) == 1:
        return "SNP"

    if len(ref) != len(alt):
        return "INDEL"

    return "MNV"


def parse_format(format_field, sample_field):
    if not format_field or not sample_field:
        return {}

    keys = format_field.split(":")
    values = sample_field.split(":")

    return dict(zip(keys, values))


def classify_genotype(gt):
    if not gt or gt in (".", "./.", ".|."):
        return "MISSING"

    alleles = gt.replace("|", "/").split("/")

    if any(a == "." for a in alleles):
        return "MISSING"

    if len(alleles) == 1:
        if alleles[0] == "0":
            return "HOM_REF"

        return "HEMI_ALT"

    if len(set(alleles)) == 1:
        if alleles[0] == "0":
            return "HOM_REF"

        return "HOM_ALT"

    return "HET"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate technical germline VCF QC metrics. "
            "This is not clinical interpretation."
        )
    )

    parser.add_argument(
        "--vcf",
        required=True,
    )

    parser.add_argument(
        "--output",
        required=True,
    )

    parser.add_argument(
        "--sample-id",
    )

    parser.add_argument(
        "--profile",
        choices=[
            "smoke",
            "wgs35x",
        ],
        default="smoke",
    )

    parser.add_argument(
        "--pipeline-version",
        default="development",
    )

    args = parser.parse_args()

    vcf_path = Path(args.vcf)

    if not vcf_path.exists():
        raise FileNotFoundError(
            f"VCF not found: {vcf_path}"
        )

    sample_name = None
    record_count = 0
    allele_count = 0

    type_counts = Counter()
    genotype_counts = Counter()
    filter_counts = Counter()
    chromosome_counts = Counter()

    qual_values = []
    dp_values = []
    gq_values = []

    transitions = 0
    transversions = 0

    multi_allelic_records = 0

    with open_text(vcf_path) as handle:
        for line in handle:

            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")

                if len(fields) >= 10:
                    sample_name = fields[9]

                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                continue

            record_count += 1

            chrom = fields[0]
            ref = fields[3]
            alts = fields[4].split(",")
            qual = safe_float(fields[5])
            filt = fields[6]

            chromosome_counts[chrom] += 1
            filter_counts[filt] += 1

            if len(alts) > 1:
                multi_allelic_records += 1

            if qual is not None:
                qual_values.append(qual)

            format_field = (
                fields[8]
                if len(fields) >= 9
                else ""
            )

            sample_field = (
                fields[9]
                if len(fields) >= 10
                else ""
            )

            sample_data = parse_format(
                format_field,
                sample_field,
            )

            gt = sample_data.get("GT")
            dp = safe_int(
                sample_data.get("DP")
            )
            gq = safe_float(
                sample_data.get("GQ")
            )

            genotype_counts[
                classify_genotype(gt)
            ] += 1

            if dp is not None:
                dp_values.append(dp)

            if gq is not None:
                gq_values.append(gq)

            for alt in alts:

                # DeepVariant gVCF placeholder should
                # never count as a biological variant.
                if alt == "<NON_REF>":
                    continue

                allele_count += 1

                vtype = variant_type(
                    ref,
                    alt,
                )

                type_counts[vtype] += 1

                if vtype == "SNP":
                    pair = (
                        ref.upper(),
                        alt.upper(),
                    )

                    if pair in TRANSITIONS:
                        transitions += 1
                    else:
                        transversions += 1

    pass_records = sum(
        count
        for key, count in filter_counts.items()
        if key in ("PASS", ".")
    )

    pass_rate = (
        pass_records / record_count
        if record_count
        else None
    )

    titv = (
        transitions / transversions
        if transversions
        else None
    )

    het = genotype_counts["HET"]
    hom_alt = genotype_counts["HOM_ALT"]

    het_hom_alt_ratio = (
        het / hom_alt
        if hom_alt
        else None
    )

    warnings = []
    failures = []

    if args.profile == "smoke":

        # Smoke profile tests technical validity only.
        # A tiny region / low-coverage BAM may correctly
        # produce zero variants.
        overall_status = "PASS"

        if record_count == 0:
            warnings.append(
                "VCF contains zero variant records. "
                "This is acceptable for the smoke profile "
                "if the tested region has insufficient "
                "coverage or no callable variants."
            )

    else:

        # Deliberately conservative development checks.
        # These are NOT clinically validated thresholds.
        overall_status = "PASS"

        if record_count == 0:
            failures.append(
                "No variant records detected in WGS profile."
            )

        if pass_rate is not None and pass_rate < 0.90:
            warnings.append(
                "PASS record fraction is below 90%."
            )

        if titv is not None and not (1.8 <= titv <= 2.3):
            warnings.append(
                "Genome-wide SNP Ti/Tv is outside the "
                "development expectation of 1.8-2.3."
            )

        if failures:
            overall_status = "FAIL"
        elif warnings:
            overall_status = "WARN"

    result = {
        "schema_version": "1.0",
        "module": "variant_qc",
        "created_at": datetime.now(
            timezone.utc
        ).isoformat(),
        "sample_id": (
            args.sample_id
            or sample_name
            or "unknown"
        ),
        "profile": args.profile,
        "pipeline_version": args.pipeline_version,
        "input": {
            "vcf": str(
                vcf_path.resolve()
            ),
        },
        "overall_status": overall_status,
        "technical_scope": (
            "Germline SNV/small-indel VCF QC. "
            "Not clinical interpretation."
        ),
        "metrics": {
            "records": record_count,
            "variant_alleles": allele_count,
            "multi_allelic_records": (
                multi_allelic_records
            ),
            "variant_types": dict(
                type_counts
            ),
            "filters": dict(
                filter_counts
            ),
            "pass_records": pass_records,
            "pass_rate": pass_rate,
            "genotypes": dict(
                genotype_counts
            ),
            "het_hom_alt_ratio": (
                het_hom_alt_ratio
            ),
            "transitions": transitions,
            "transversions": transversions,
            "ti_tv": titv,
            "quality": {
                "qual": {
                    "n": len(
                        qual_values
                    ),
                    "mean": (
                        sum(qual_values)
                        / len(qual_values)
                        if qual_values
                        else None
                    ),
                    "median": percentile(
                        qual_values,
                        0.5,
                    ),
                    "p10": percentile(
                        qual_values,
                        0.1,
                    ),
                    "p90": percentile(
                        qual_values,
                        0.9,
                    ),
                },
                "depth": {
                    "n": len(
                        dp_values
                    ),
                    "mean": (
                        sum(dp_values)
                        / len(dp_values)
                        if dp_values
                        else None
                    ),
                    "median": percentile(
                        dp_values,
                        0.5,
                    ),
                    "p10": percentile(
                        dp_values,
                        0.1,
                    ),
                    "p90": percentile(
                        dp_values,
                        0.9,
                    ),
                },
                "genotype_quality": {
                    "n": len(
                        gq_values
                    ),
                    "mean": (
                        sum(gq_values)
                        / len(gq_values)
                        if gq_values
                        else None
                    ),
                    "median": percentile(
                        gq_values,
                        0.5,
                    ),
                    "p10": percentile(
                        gq_values,
                        0.1,
                    ),
                    "p90": percentile(
                        gq_values,
                        0.9,
                    ),
                },
            },
            "records_by_chromosome": dict(
                chromosome_counts
            ),
        },
        "warnings": warnings,
        "failures": failures,
        "threshold_note": (
            "wgs35x thresholds are development "
            "heuristics and are not clinically validated."
        ),
    }

    output = Path(args.output)

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output.write_text(
        json.dumps(
            result,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "sample_id": result["sample_id"],
                "status": overall_status,
                "records": record_count,
                "snps": type_counts["SNP"],
                "indels": type_counts["INDEL"],
                "ti_tv": titv,
                "output": str(output),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
