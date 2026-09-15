#!/usr/bin/env bash
set -euo pipefail

export PYTHONPATH="${PYTHONPATH:-}:src"

CLINVAR_DB="results/evidence/clinvar/clinvar.GRCh38.sqlite"

if [[ ! -f "$CLINVAR_DB" ]]; then
    echo "ERROR: ClinVar database not found: $CLINVAR_DB"
    exit 1
fi

python - <<'PY'
import json

from genomics_platform.evidence.clinvar_store import query_variant
from genomics_platform.evidence.adapters.gnomad_api_adapter import query_variant as query_gnomad
from genomics_platform.evidence.variant_evidence_builder import (
    build_variant_evidence_bundle,
)
from genomics_platform.evidence.variant_evidence_contract_builder import (
    build_variant_evidence_contract,
)


CLINVAR_DB = "results/evidence/clinvar/clinvar.GRCh38.sqlite"

assembly = "GRCh38"
chrom = "1"
pos = 66926
ref = "AG"
alt = "A"


# ---------------------------------------------------------
# 1. REAL ClinVar evidence
# ---------------------------------------------------------

clinvar_result = query_variant(
    CLINVAR_DB,
    chrom,
    pos,
    ref,
    alt,
)

assert clinvar_result["status"] == "FOUND"
assert len(clinvar_result["records"]) >= 1

variation_ids = {
    record["variation_id"]
    for record in clinvar_result["records"]
}

assert "3385321" in variation_ids


# ---------------------------------------------------------
# 2. LIVE gnomAD lookup
#
# FOUND / NOT_FOUND are both valid completed lookups.
# ERROR is preserved as ERROR and must not be converted
# into NOT_FOUND.
# ---------------------------------------------------------

gnomad_result = query_gnomad(
    chrom=chrom,
    pos=pos,
    ref=ref,
    alt=alt,
    assembly=assembly,
)


# ---------------------------------------------------------
# 3. Synthetic VEP fixture
#
# This tests plumbing only.
# No biological interpretation is allowed from this fixture.
# ---------------------------------------------------------

vep_annotations = [
    {
        "Allele": "-",
        "Consequence": "intron_variant",
        "IMPACT": "MODIFIER",
        "SYMBOL": "OR4F5",
        "Gene": "ENSG00000186092",
        "Feature_type": "Transcript",
        "Feature": "ENST00000641515",
        "BIOTYPE": "protein_coding",
        "ENSP": "ENSP00000493376",
        "CANONICAL": "YES",
        "VARIANT_CLASS": "deletion",
        "_normalized": {
            "allele": "-",
            "consequence": "intron_variant",
            "impact": "MODIFIER",
            "symbol": "OR4F5",
            "gene_id": "ENSG00000186092",
            "feature_type": "Transcript",
            "transcript_id": "ENST00000641515",
            "biotype": "protein_coding",
            "protein_id": "ENSP00000493376",
            "canonical": True,
            "variant_class": "deletion",
        },
    }
]


# ---------------------------------------------------------
# 4. Build Unified Variant Evidence Bundle
# ---------------------------------------------------------

bundle = build_variant_evidence_bundle(
    assembly=assembly,
    chrom=chrom,
    pos=pos,
    ref=ref,
    alt=alt,
    vep_annotations=vep_annotations,
    clinvar_records=clinvar_result,
    gnomad_record=gnomad_result,
    vep_release="116",
    vep_provenance={
        "mode": "synthetic_fixture",
        "biological_interpretation_allowed": False,
        "parser_version": "0.1.0",
    },
)


# ---------------------------------------------------------
# 5. Bundle -> Stable Evidence Contract
# ---------------------------------------------------------

contract = build_variant_evidence_contract(bundle)

data = contract.to_dict()


# ---------------------------------------------------------
# 6. Identity assertions
# ---------------------------------------------------------

assert data["variant_key"] == "1:66926:AG:A"
assert data["identity"]["assembly"] == "GRCh38"
assert data["identity"]["chrom"] == "1"
assert data["identity"]["pos"] == 66926
assert data["identity"]["ref"] == "AG"
assert data["identity"]["alt"] == "A"


# ---------------------------------------------------------
# 7. Functional evidence assertions
# ---------------------------------------------------------

assert data["availability"]["VEP"]["status"] == "FOUND"

assert "OR4F5" in data["functional"]["genes"]
assert "intron_variant" in data["functional"]["consequences"]

assert (
    data["functional"]["canonical_transcript_count"]
    == 1
)


# ---------------------------------------------------------
# 8. Clinical evidence assertions
# ---------------------------------------------------------

assert (
    data["availability"]["ClinVar"]["status"]
    == "FOUND"
)

contract_variation_ids = {
    record["variation_id"]
    for record in data["clinical"]["records"]
}

assert "3385321" in contract_variation_ids


# ---------------------------------------------------------
# 9. gnomAD semantic preservation
# ---------------------------------------------------------

gnomad_status = data["availability"]["gnomAD"]["status"]

assert gnomad_status in {
    "FOUND",
    "NOT_FOUND",
    "ERROR",
}

if gnomad_status == "NOT_FOUND":
    assert data["population"]["exome"] is None
    assert data["population"]["genome"] is None

if gnomad_status == "ERROR":
    assert data["evidence_complete"] is False
    assert "gnomAD" in data["error_sources"]


# ---------------------------------------------------------
# 10. Provenance assertions
# ---------------------------------------------------------

provenance = {
    item["source"]: item
    for item in data["provenance"]
}

assert "VEP" in provenance
assert "ClinVar" in provenance
assert "gnomAD" in provenance

assert (
    provenance["VEP"]["metadata"]["mode"]
    == "synthetic_fixture"
)

assert (
    provenance["VEP"]["metadata"][
        "biological_interpretation_allowed"
    ]
    is False
)


# ---------------------------------------------------------
# 11. Safety assertions
# ---------------------------------------------------------

serialized = json.dumps(data).lower()

for forbidden in (
    '"ranking_score"',
    '"pathogenicity_probability"',
    '"acmg_classification"',
    '"tier"',
):
    assert forbidden not in serialized


print()
print(
    "PASS: VARIANT EVIDENCE CONTRACT "
    "INTEGRATION SMOKE v1.0"
)
print(f"variant: {data['variant_key']}")
print(
    "ClinVar:",
    data["availability"]["ClinVar"]["status"],
    f"({data['clinical']['record_count']} record(s))",
)
print(
    "ClinVar Variation IDs:",
    sorted(contract_variation_ids),
)
print(
    "VEP:",
    data["availability"]["VEP"]["status"],
    "(synthetic fixture; biological interpretation disabled)",
)
print(
    "gnomAD:",
    gnomad_status,
)
print(
    "evidence_complete:",
    data["evidence_complete"],
)
print(
    "error_sources:",
    data["error_sources"],
)

print()
print(
    json.dumps(
        {
            "variant_key": data["variant_key"],
            "contract_schema_version": data[
                "schema_version"
            ],
            "evidence_complete": data[
                "evidence_complete"
            ],
            "error_sources": data[
                "error_sources"
            ],
            "statuses": {
                source: value["status"]
                for source, value
                in data["availability"].items()
            },
            "functional": {
                "genes": data["functional"]["genes"],
                "consequences": data[
                    "functional"
                ]["consequences"],
            },
            "clinical": {
                "record_count": data[
                    "clinical"
                ]["record_count"],
                "variation_ids": sorted(
                    contract_variation_ids
                ),
            },
        },
        indent=2,
    )
)
PY
