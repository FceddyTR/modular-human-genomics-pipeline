#!/usr/bin/env bash
set -euo pipefail

export PYTHONPATH="${PYTHONPATH:-}:$(pwd)/src"

CLINVAR_DB="results/evidence/clinvar/clinvar.GRCh38.sqlite"

if [[ ! -f "$CLINVAR_DB" ]]; then
    echo "ERROR: ClinVar database not found: $CLINVAR_DB" >&2
    exit 1
fi

python - <<'PY'
import json

from genomics_platform.evidence.clinvar_store import query_variant
from genomics_platform.evidence.adapters.gnomad_api_adapter import (
    query_variant as query_gnomad,
)
from genomics_platform.evidence.variant_evidence_builder import (
    build_variant_evidence_bundle,
)


CLINVAR_DB = "results/evidence/clinvar/clinvar.GRCh38.sqlite"

CHROM = "1"
POS = 66926
REF = "AG"
ALT = "A"


# ------------------------------------------------------------
# REAL ClinVar local evidence
# ------------------------------------------------------------

clinvar_result = query_variant(
    CLINVAR_DB,
    CHROM,
    POS,
    REF,
    ALT,
)

assert clinvar_result["status"] == "FOUND"
assert clinvar_result["records"], (
    "Expected known ClinVar test variant to be FOUND"
)


# ------------------------------------------------------------
# REAL gnomAD live API evidence
#
# FOUND / NOT_FOUND / ERROR are all valid transport outcomes.
# ERROR must remain ERROR and must never become NOT_FOUND.
# ------------------------------------------------------------

gnomad_record = query_gnomad(
    CHROM,
    POS,
    REF,
    ALT,
)

assert gnomad_record.lookup_status in {
    "FOUND",
    "NOT_FOUND",
    "ERROR",
}


# ------------------------------------------------------------
# CONTROLLED VEP fixture
#
# The real VEP cache is not yet active in this integration test.
# This MUST NOT be interpreted as biological annotation.
# ------------------------------------------------------------

vep_annotations = [
    {
        "fixture": True,
        "note": (
            "Synthetic VEP component used only for integration plumbing. "
            "Not a biological annotation for this variant."
        ),
    }
]


bundle = build_variant_evidence_bundle(
    assembly="GRCh38",
    chrom=CHROM,
    pos=POS,
    ref=REF,
    alt=ALT,
    vep_annotations=vep_annotations,
    clinvar_records=clinvar_result,
    gnomad_record=gnomad_record,
    vep_release="116",
    vep_provenance={
        "mode": "synthetic_fixture",
        "biological_interpretation_allowed": False,
    },
)

data = bundle.to_dict()


# ------------------------------------------------------------
# Variant identity
# ------------------------------------------------------------

assert data["variant_key"] == "1:66926:AG:A"


# ------------------------------------------------------------
# ClinVar
# ------------------------------------------------------------

clinical = data["evidence"]["clinical"]

assert clinical["status"] == "FOUND"
assert len(clinical["data"]) >= 1

variation_ids = {
    str(record.get("variation_id"))
    for record in clinical["data"]
}

assert "3385321" in variation_ids


# ------------------------------------------------------------
# VEP provenance safety
# ------------------------------------------------------------

functional = data["evidence"]["functional"]

assert functional["status"] == "FOUND"
assert functional["provenance"]["mode"] == "synthetic_fixture"
assert (
    functional["provenance"]["biological_interpretation_allowed"]
    is False
)


# ------------------------------------------------------------
# gnomAD semantics
# ------------------------------------------------------------

population = data["evidence"]["population"]

assert population["status"] == gnomad_record.lookup_status


if population["status"] == "FOUND":

    assert population["data"] is not None

    exome = population["data"]["exome"]
    genome = population["data"]["genome"]

    if exome is not None and exome["an"]:
        assert exome["af"] == exome["ac"] / exome["an"]

    if genome is not None and genome["an"]:
        assert genome["af"] == genome["ac"] / genome["an"]


elif population["status"] == "NOT_FOUND":

    assert population["data"] is not None

    # NOT_FOUND is absence of a matching gnomAD record,
    # never an inferred zero frequency.
    assert population["data"]["exome"] is None
    assert population["data"]["genome"] is None
    assert population["error"] is None


elif population["status"] == "ERROR":

    # Critical safety invariant:
    # external service failure must survive aggregation as ERROR.
    assert population["error"] is not None
    assert data["evidence_complete"] is False
    assert "gnomAD" in data["error_sources"]


# ------------------------------------------------------------
# Technical completeness
# ------------------------------------------------------------

if population["status"] == "ERROR":
    assert data["evidence_complete"] is False
else:
    assert data["evidence_complete"] is True


# ------------------------------------------------------------
# No interpretation leakage
# ------------------------------------------------------------

assert "classification" not in data
assert "tier" not in data
assert "score" not in data
assert "pathogenicity_probability" not in data


print("PASS: UNIFIED VARIANT EVIDENCE INTEGRATION SMOKE v0.1")
print("variant:", data["variant_key"])

print(
    "ClinVar:",
    clinical["status"],
    f"({len(clinical['data'])} record(s))",
)

print(
    "ClinVar Variation IDs:",
    sorted(variation_ids),
)

print(
    "VEP:",
    functional["status"],
    "(synthetic fixture; biological interpretation disabled)",
)

print(
    "gnomAD:",
    population["status"],
)


if population["status"] == "FOUND":

    exome = population["data"]["exome"]
    genome = population["data"]["genome"]

    print(
        "gnomAD exome AF:",
        None if exome is None else exome["af"],
    )

    print(
        "gnomAD genome AF:",
        None if genome is None else genome["af"],
    )


elif population["status"] == "ERROR":

    print(
        "gnomAD live-service warning:",
        population["error"],
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
            "evidence_complete": data["evidence_complete"],
            "error_sources": data["error_sources"],
            "statuses": {
                key: value["status"]
                for key, value in data["evidence"].items()
            },
        },
        indent=2,
    )
)
PY
