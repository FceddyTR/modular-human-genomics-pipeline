#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INPUT="${ROOT}/tests/fixtures/vep/synthetic.vep.vcf"
OUTPUT="${ROOT}/tests/fixtures/vep/observed.vep.json"

python \
  "${ROOT}/src/genomics_platform/annotation/vep_parser.py" \
  --input-vcf "${INPUT}" \
  --output-json "${OUTPUT}"

python - "${OUTPUT}" <<'PY'
import json
import sys

path = sys.argv[1]

with open(path, encoding="utf-8") as handle:
    data = json.load(handle)

summary = data["summary"]

print()
print("VEP parser assertions:")
print(f"  variants       = {summary['variant_records']}")
print(f"  transcripts    = {summary['transcript_annotations']}")
print(f"  canonical      = {summary['canonical_transcripts']}")
print(f"  unique genes   = {summary['unique_genes']}")
print(f"  MODERATE       = {summary['impact_counts'].get('MODERATE', 0)}")
print(f"  LOW            = {summary['impact_counts'].get('LOW', 0)}")
print(f"  MODIFIER       = {summary['impact_counts'].get('MODIFIER', 0)}")

assert data["schema_version"] == "1.0"

assert summary["variant_records"] == 2
assert summary["transcript_annotations"] == 3
assert summary["canonical_transcripts"] == 2
assert summary["unique_genes"] == 3

assert summary["impact_counts"]["MODERATE"] == 1
assert summary["impact_counts"]["LOW"] == 1
assert summary["impact_counts"]["MODIFIER"] == 1

assert summary["consequence_counts"]["missense_variant"] == 1
assert summary["consequence_counts"]["synonymous_variant"] == 1
assert summary["consequence_counts"]["downstream_gene_variant"] == 1

variant1 = data["variants"][0]

assert variant1["variant_key"] == "chr20:10000000:A:G"

tx = variant1["transcript_annotations"][0]["_normalized"]

assert tx["symbol"] == "GENE1"
assert tx["impact"] == "MODERATE"
assert tx["canonical"] is True
assert tx["protein_id"] == "ENSP000001"

assert data["scope"]["clinical_classification"] is False
assert data["scope"]["candidate_ranking"] is False

print()
print("==============================================")
print(" PASS: VEP PARSER")
print("==============================================")
PY
