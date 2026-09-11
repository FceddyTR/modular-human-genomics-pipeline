#!/usr/bin/env python3

import argparse
import json
import sys
from pathlib import Path


def detect_status(data: dict) -> str:
    candidates = [
        data.get("gate"),
        data.get("overall_status"),
        data.get("status"),
        data.get("qc_status"),
    ]

    for value in candidates:
        if isinstance(value, str):
            value = value.upper().strip()
            if value in {"PASS", "WARN", "FAIL"}:
                return value

    raise ValueError(
        "Could not find PASS/WARN/FAIL status in QC JSON. "
        "Checked: gate, overall_status, status, qc_status."
    )


def main():
    parser = argparse.ArgumentParser(
        description="Gate downstream analysis based on QC result."
    )
    parser.add_argument("--qc-json", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    qc_path = Path(args.qc_json)

    if not qc_path.exists():
        raise FileNotFoundError(qc_path)

    with qc_path.open() as handle:
        data = json.load(handle)

    status = detect_status(data)

    if status == "FAIL":
        print(
            f"[QC_GATE] {args.sample_id}: RAW QC FAIL. "
            "Downstream analysis blocked.",
            file=sys.stderr,
        )
        sys.exit(2)

    output = Path(args.output)

    output.write_text(
        f"sample_id\t{args.sample_id}\n"
        f"raw_qc_status\t{status}\n"
        f"downstream_allowed\ttrue\n"
    )

    print(
        f"[QC_GATE] {args.sample_id}: {status}. "
        "Downstream analysis allowed."
    )


if __name__ == "__main__":
    main()
