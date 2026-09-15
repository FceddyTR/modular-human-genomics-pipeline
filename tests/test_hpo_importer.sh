#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INPUT="${ROOT}/tests/fixtures/hpo/genes_to_phenotype.txt"
OUTPUT="${ROOT}/tests/fixtures/hpo/observed.hpo.json"
INDEX="${ROOT}/tests/fixtures/hpo/observed.hpo.index.json"

PYTHONPATH="${ROOT}" python \
  "${ROOT}/src/genomics_platform/evidence/adapters/hpo_importer.py" \
  --input "${INPUT}" \
  --output-json "${OUTPUT}" \
  --index-json "${INDEX}"

python - "${OUTPUT}" "${INDEX}" <<'PY'
import json
import sys

output_path = sys.argv[1]
index_path = sys.argv[2]

with open(output_path, encoding="utf-8") as f:
    data = json.load(f)

with open(index_path, encoding="utf-8") as f:
    index = json.load(f)

s = data["summary"]

print()
print("HPO assertions:")
print("  records           =", s["record_count"])
print("  unique genes      =", s["unique_genes"])
print("  unique phenotypes =", s["unique_phenotypes"])
print("  unique diseases   =", s["unique_diseases"])
print("  indexed genes     =", index["gene_count"])

assert s["record_count"] == 3
assert s["unique_genes"] == 2
assert s["unique_phenotypes"] == 3
assert s["unique_diseases"] == 3

assert s["disease_prefixes"]["OMIM"] == 2
assert s["disease_prefixes"]["ORPHA"] == 1

assert index["gene_count"] == 2
assert len(index["by_gene"]["TESTGENE"]) == 2

r = index["by_gene"]["TESTGENE"][0]

assert r["schema_version"] == "1.1"
assert r["source"] == "HPO"
assert r["evidence_type"] == "gene_phenotype"

assert r["gene_symbol"] == "TESTGENE"
assert r["gene_id"] == "NCBIGene:11111"

assert r["phenotype_ids"] == ["HP:0000001"]
assert r["disease_id"] == "OMIM:999001"

assert r["confidence"] is None
assert r["frequency"] == "2/3"

print()
print("==============================================")
print(" PASS: HPO IMPORTER v1.1")
print("==============================================")
PY
