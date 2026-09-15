#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INPUT="${ROOT}/tests/fixtures/panelapp/panelapp_fixture.json"
OUTPUT="${ROOT}/tests/fixtures/panelapp/observed.panelapp.json"

PYTHONPATH="${ROOT}" python \
  "${ROOT}/src/genomics_platform/evidence/adapters/panelapp_adapter.py" \
  --gene TESTGENE \
  --input-json "${INPUT}" \
  --output-json "${OUTPUT}"

python - "${OUTPUT}" <<'PY'
import json
import sys

path = sys.argv[1]

with open(path, encoding="utf-8") as handle:
    data = json.load(handle)

print()
print("PanelApp assertions:")
print("  status      =", data["status"])
print("  records     =", data["summary"]["record_count"])
print("  confidence  =", data["summary"]["confidence_counts"])

assert data["status"] == "FOUND"
assert data["summary"]["record_count"] == 1
assert data["summary"]["confidence_counts"]["green"] == 1

record = data["records"][0]

assert record["source"] == "PanelApp"
assert record["gene_symbol"] == "TESTGENE"
assert record["gene_id"] == "HGNC:999999"

assert record["confidence"] == "green"

assert record["inheritance"] == [
    "MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted"
]

assert record["source_version"] == "2.1"

assert "Synthetic inflammatory disease" in (
    record["disease_name"] or ""
)

assert data["schema_version"] == "1.0"

print()
print("==============================================")
print(" PASS: PANELAPP ADAPTER")
print("==============================================")
PY
