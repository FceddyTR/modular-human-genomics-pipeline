#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INPUT="${ROOT}/tests/fixtures/gencc/gencc_new_format.csv"
OUTPUT="${ROOT}/tests/fixtures/gencc/observed.gencc.json"

PYTHONPATH="${ROOT}" python \
  "${ROOT}/src/genomics_platform/evidence/adapters/gencc_importer.py" \
  --input-csv "${INPUT}" \
  --output-json "${OUTPUT}"

python - "${OUTPUT}" <<'PY'
import json
import sys

path = sys.argv[1]

with open(path, encoding="utf-8") as handle:
    data = json.load(handle)

summary = data["summary"]

print()
print("GenCC assertions:")
print(f"  records         = {summary['record_count']}")
print(f"  unique genes    = {summary['unique_genes']}")
print(f"  unique diseases = {summary['unique_diseases']}")
print(f"  classifications = {summary['classifications']}")

assert summary["record_count"] == 2
assert summary["unique_genes"] == 2
assert summary["unique_diseases"] == 2

assert summary["classifications"]["Definitive"] == 1
assert summary["classifications"]["Strong"] == 1

r1 = data["records"][0]

assert r1["source"] == "GenCC"
assert r1["source_record_id"] == "SGC-TEST001"
assert r1["source_version"] == "2"
assert r1["gene_symbol"] == "TESTGENE1"
assert r1["gene_id"] == "HGNC:11111"
assert r1["disease_id"] == "MONDO:0000001"
assert r1["inheritance"] == ["Autosomal dominant"]
assert r1["classification"] == "Definitive"
assert r1["license"] == "CC0-1.0"

print()
print("==============================================")
print(" PASS: GENCC IMPORTER")
print("==============================================")
PY
