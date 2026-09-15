from __future__ import annotations

import json
import urllib.error
import urllib.request
from typing import Any, Optional

from genomics_platform.evidence.gnomad_record import (
    DatasetFrequency,
    GnomADRecord,
    PopulationFrequency,
)


GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
ADAPTER_VERSION = "0.1.0"

QUERY = """
query Variant($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variant_id
    chrom
    pos
    ref
    alt

    exome {
      ac
      an
      homozygote_count
      populations {
        id
        ac
        an
        homozygote_count
      }
    }

    genome {
      ac
      an
      homozygote_count
      populations {
        id
        ac
        an
        homozygote_count
      }
    }
  }
}
"""


def canonical_chrom(chrom: str) -> str:
    value = str(chrom).strip()

    if value.lower().startswith("chr"):
        value = value[3:]

    if value.upper() == "M":
        return "MT"

    return value


def _population(row: dict[str, Any]) -> PopulationFrequency:
    return PopulationFrequency(
        id=str(row["id"]),
        ac=row.get("ac"),
        an=row.get("an"),
        homozygote_count=row.get("homozygote_count"),
    )


def _dataset(row: Optional[dict[str, Any]]) -> Optional[DatasetFrequency]:
    if row is None:
        return None

    populations = tuple(
        _population(item)
        for item in row.get("populations", [])
        if item.get("id") is not None
    )

    return DatasetFrequency(
        ac=row.get("ac"),
        an=row.get("an"),
        homozygote_count=row.get("homozygote_count"),
        populations=populations,
    )


def query_variant(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    *,
    dataset_id: str = "gnomad_r4",
    assembly: str = "GRCh38",
    timeout: int = 30,
    api_url: str = GNOMAD_API_URL,
) -> GnomADRecord:

    chrom = canonical_chrom(chrom)
    pos = int(pos)
    ref = ref.upper()
    alt = alt.upper()

    variant_id = f"{chrom}-{pos}-{ref}-{alt}"

    payload = {
        "query": QUERY,
        "variables": {
            "variantId": variant_id,
            "datasetId": dataset_id,
        },
    }

    request = urllib.request.Request(
        api_url,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "User-Agent": "modular-human-genomics-pipeline/0.1",
        },
        method="POST",
    )

    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            result = json.load(response)

    except (
        urllib.error.URLError,
        urllib.error.HTTPError,
        TimeoutError,
        OSError,
        json.JSONDecodeError,
    ) as exc:
        return GnomADRecord(
            assembly=assembly,
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            lookup_status="ERROR",
            source_dataset=dataset_id,
            error=f"{type(exc).__name__}: {exc}",
        )

    if result.get("errors"):
        messages = "; ".join(
            str(item.get("message", item))
            for item in result["errors"]
        )

        return GnomADRecord(
            assembly=assembly,
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            lookup_status="ERROR",
            source_dataset=dataset_id,
            error=f"GraphQL: {messages}",
        )

    variant = result.get("data", {}).get("variant")

    if variant is None:
        return GnomADRecord(
            assembly=assembly,
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            lookup_status="NOT_FOUND",
            source_dataset=dataset_id,
        )

    return GnomADRecord(
        assembly=assembly,
        chrom=str(variant.get("chrom", chrom)),
        pos=int(variant.get("pos", pos)),
        ref=str(variant.get("ref", ref)),
        alt=str(variant.get("alt", alt)),
        lookup_status="FOUND",
        exome=_dataset(variant.get("exome")),
        genome=_dataset(variant.get("genome")),
        source_dataset=dataset_id,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Query one normalized variant against gnomAD."
    )

    parser.add_argument("chrom")
    parser.add_argument("pos", type=int)
    parser.add_argument("ref")
    parser.add_argument("alt")
    parser.add_argument("--dataset", default="gnomad_r4")

    args = parser.parse_args()

    record = query_variant(
        args.chrom,
        args.pos,
        args.ref,
        args.alt,
        dataset_id=args.dataset,
    )

    print(json.dumps(record.to_dict(), indent=2))
