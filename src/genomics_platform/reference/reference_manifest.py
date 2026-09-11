#!/usr/bin/env python3

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


VERSION = "0.1.0"
AUTHOR = "Efe Mıhcı"


def sha256_file(path, chunk_size=16 * 1024 * 1024):
    digest = hashlib.sha256()

    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)

            if not chunk:
                break

            digest.update(chunk)

    return digest.hexdigest()


def parse_fai(path):
    contigs = []

    with open(path) as handle:
        for line in handle:
            fields = line.rstrip().split("\t")

            if len(fields) < 2:
                continue

            contigs.append(
                {
                    "name": fields[0],
                    "length": int(fields[1]),
                }
            )

    return contigs


def main():
    parser = argparse.ArgumentParser(
        description="Generate reproducible reference manifest."
    )

    parser.add_argument(
        "--fasta",
        required=True,
    )

    parser.add_argument(
        "--output",
        required=True,
    )

    parser.add_argument(
        "--hash-indexes",
        action="store_true",
    )

    args = parser.parse_args()

    fasta = Path(args.fasta).resolve()
    fai = Path(str(fasta) + ".fai")

    if not fasta.exists():
        raise FileNotFoundError(fasta)

    if not fai.exists():
        raise FileNotFoundError(fai)

    bwa_extensions = [
        "amb",
        "ann",
        "bwt",
        "pac",
        "sa",
    ]

    bwa_indexes = []

    for extension in bwa_extensions:
        path = Path(
            str(fasta) + "." + extension
        )

        entry = {
            "file": path.name,
            "exists": path.exists(),
            "size_bytes": (
                path.stat().st_size
                if path.exists()
                else None
            ),
        }

        if args.hash_indexes and path.exists():
            entry["sha256"] = sha256_file(path)

        bwa_indexes.append(entry)

    contigs = parse_fai(fai)

    observed = {item["name"] for item in contigs}

    expected_primary = {
        *(f"chr{i}" for i in range(1, 23)),
        "chrX",
        "chrY",
        "chrM",
    }

    missing_primary = sorted(
        expected_primary - observed
    )

    manifest = {
        "manifest_version": VERSION,
        "author": AUTHOR,
        "generated_at": datetime.now(
            timezone.utc
        ).isoformat(),
        "reference": {
            "file": fasta.name,
            "size_bytes": fasta.stat().st_size,
            "sha256": sha256_file(fasta),
            "fai_file": fai.name,
            "fai_sha256": sha256_file(fai),
            "contig_count": len(contigs),
            "total_reference_bases": sum(
                item["length"]
                for item in contigs
            ),
            "missing_primary_contigs": missing_primary,
        },
        "bwa_index": bwa_indexes,
        "status": (
            "PASS"
            if not missing_primary
            and all(
                item["exists"]
                for item in bwa_indexes
            )
            else "FAIL"
        ),
    }

    output = Path(args.output)
    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output.write_text(
        json.dumps(
            manifest,
            indent=2,
        )
    )

    print(
        f"Reference manifest: {output}"
    )

    print(
        f"Status: {manifest['status']}"
    )

    print(
        f"SHA256: "
        f"{manifest['reference']['sha256']}"
    )


if __name__ == "__main__":
    main()
