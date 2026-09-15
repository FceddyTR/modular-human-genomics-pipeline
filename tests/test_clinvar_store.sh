#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

cat > "$TMP/clinvar.vcf" <<'VCF'
##fileformat=VCFv4.2
##fileDate=2026-09-13
##source=ClinVar
##reference=GRCh38
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	66926	3385321	AG	A	.	.	ALLELEID=3544463;CLNDISDB=MONDO:MONDO:0019200,OMIM:268000,Orphanet:791;CLNDN=Retinitis_pigmentosa;CLNHGVS=NC_000001.11:g.66927del;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNSIGSCV=SCV005419006;CLNVC=Deletion;CLNVCSO=SO:0000159;GENEINFO=OR4F5:79501;MC=SO:0001627|intron_variant;ORIGIN=0
VCF

python -m src.genomics_platform.evidence.clinvar_store build \
  --input-vcf "$TMP/clinvar.vcf" \
  --database "$TMP/clinvar.sqlite"

python - "$TMP/clinvar.sqlite" <<'PY'
import sys

from src.genomics_platform.evidence.clinvar_store import (
    query_variant,
    query_variation_id,
)

db = sys.argv[1]

found = query_variant(
    db,
    "1",
    66926,
    "AG",
    "A",
)

assert found["status"] == "FOUND"
assert found["interpretation"] == "CLINVAR_RECORD_PRESENT"
assert len(found["records"]) == 1

r = found["records"][0]

assert r["schema_version"] == "1.1"
assert r["variant_key"] == "1:66926:AG:A"
assert r["variation_id"] == "3385321"
assert r["vcv_accession"] == "VCV003385321"
assert r["allele_id"] == "3544463"

assert r["gene_symbols"] == ["OR4F5"]
assert r["gene_ids"] == ["79501"]

assert r["clinical_significance"] == [
    "Uncertain_significance"
]

assert r["conflicting_significance"] == []

assert r["review_status"] == (
    "criteria_provided,_single_submitter"
)

assert r["review_stars"] == 1

assert r["scv_accessions"] == [
    "SCV005419006"
]

assert r["molecular_consequences"] == [
    {
        "so_id": "SO:0001627",
        "consequence": "intron_variant",
    }
]

by_id = query_variation_id(
    db,
    "3385321",
)

assert len(by_id) == 1
assert by_id[0]["variant_key"] == "1:66926:AG:A"

# chr-prefix normalization
with_chr = query_variant(
    db,
    "chr1",
    66926,
    "AG",
    "A",
)

assert with_chr["status"] == "FOUND"

# Critical semantic regression:
# absence from ClinVar must NOT become benign/rejected.
missing = query_variant(
    db,
    "1",
    999999,
    "A",
    "G",
)

assert missing["status"] == "NOT_FOUND"
assert missing["interpretation"] == "NO_CLINVAR_RECORD"
assert missing["records"] == []

assert "benign" not in str(missing).lower()
assert "reject" not in str(missing).lower()

print("")
print("==============================================")
print(" PASS: CLINVAR STORE v0.1")
print(" FOUND / NOT_FOUND semantics preserved")
print("==============================================")
PY
