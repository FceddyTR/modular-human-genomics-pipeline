"""Build the frozen public benchmark cohort manifest.

This script performs cohort construction only. It does not run ranking
and does not create truth-blind real-case ranking manifests.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from genomics_platform.ranking.benchmark.public_cases import (
    build_public_cohort_manifest,
    load_eligible_source_records,
    select_diverse_cohort,
    split_public_benchmark,
)


DEFAULT_SOURCE = Path(
    "data/benchmark_sources/"
    "phenopacket_store/release/0.1.27"
)

DEFAULT_OUTPUT = Path(
    "results/public_benchmark/"
    "v0.1/cohort_manifest.json"
)


def canonical_json_bytes(
    payload: dict,
) -> bytes:
    return (
        json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        )
        + "\n"
    ).encode("utf-8")


def manifest_sha256(
    payload: dict,
) -> str:
    return hashlib.sha256(
        canonical_json_bytes(payload)
    ).hexdigest()


def build(
    source: Path,
) -> dict:
    records = load_eligible_source_records(
        source
    )

    cohort = select_diverse_cohort(
        records,
        target_size=500,
    )

    split = split_public_benchmark(
        cohort,
        development_size=300,
        locked_size=200,
    )

    manifest = (
        build_public_cohort_manifest(
            split
        )
    )

    return manifest.to_dict()


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--source",
        type=Path,
        default=DEFAULT_SOURCE,
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
    )

    args = parser.parse_args()

    payload = build(args.source)

    digest = manifest_sha256(
        payload
    )

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.output.write_bytes(
        canonical_json_bytes(payload)
    )

    fingerprint_path = (
        args.output.with_suffix(
            args.output.suffix + ".sha256"
        )
    )

    fingerprint_path.write_text(
        f"{digest}  {args.output.name}\n",
        encoding="utf-8",
    )

    print("=" * 78)
    print("PUBLIC BENCHMARK COHORT BUILD")
    print("=" * 78)
    print(
        "Development:",
        len(payload["development"]),
    )
    print(
        "Locked evaluation:",
        len(
            payload[
                "locked_evaluation"
            ]
        ),
    )
    print(
        "Total:",
        payload["cohort"][
            "total_cases"
        ],
    )
    print("Manifest:", args.output)
    print("SHA-256:", digest)
    print(
        "Fingerprint:",
        fingerprint_path,
    )
    print()
    print(
        "WARNING: manifest contains "
        "benchmark truth and is not "
        "ranking input."
    )


if __name__ == "__main__":
    main()
