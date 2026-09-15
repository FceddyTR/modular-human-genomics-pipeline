#!/usr/bin/env bash
set -euo pipefail

export PYTHONPATH="${PYTHONPATH:-}:src"

python - <<'PY'
import json
from pathlib import Path

from genomics_platform.evidence.clinvar_store import (
    query_variant as query_clinvar,
)
from genomics_platform.evidence.adapters.gnomad_api_adapter import (
    query_variant as query_gnomad,
)
from genomics_platform.evidence.variant_evidence_builder import (
    build_variant_evidence_bundle,
)
from genomics_platform.evidence.variant_evidence_contract_builder import (
    build_variant_evidence_contract,
)
from genomics_platform.interpretation.interpretation_engine import (
    interpret_variant,
)


# ------------------------------------------------------------
# Test variant
# ------------------------------------------------------------

ASSEMBLY = "GRCh38"
CHROM = "1"
POS = 66926
REF = "AG"
ALT = "A"

VARIANT_KEY = f"{CHROM}:{POS}:{REF}:{ALT}"

CLINVAR_DB = Path(
    "results/evidence/clinvar/clinvar.GRCh38.sqlite"
)

if not CLINVAR_DB.exists():
    raise SystemExit(
        f"FAIL: ClinVar SQLite database not found: {CLINVAR_DB}"
    )


# ------------------------------------------------------------
# Real ClinVar lookup
# ------------------------------------------------------------

clinvar_result = query_clinvar(
    CLINVAR_DB,
    CHROM,
    POS,
    REF,
    ALT,
)

assert clinvar_result["status"] == "FOUND"
assert clinvar_result["records"]

variation_ids = [
    record["variation_id"]
    for record in clinvar_result["records"]
]

assert "3385321" in variation_ids


# ------------------------------------------------------------
# Live gnomAD lookup
#
# Adapter returns GnomADRecord.
# ------------------------------------------------------------

gnomad_record = query_gnomad(
    chrom=CHROM,
    pos=POS,
    ref=REF,
    alt=ALT,
)

gnomad_status = gnomad_record.lookup_status

assert gnomad_status in {
    "FOUND",
    "NOT_FOUND",
    "ERROR",
}


# ------------------------------------------------------------
# Synthetic VEP fixture
#
# This is intentionally synthetic and tests integration plumbing.
# It must not be interpreted as a biological annotation of the
# real ClinVar variant.
#
# The structure mirrors the already passing evidence integration
# smoke test and the normalized VEP contract.
# ------------------------------------------------------------

vep_annotations = [
    {
        "Allele": "A",
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
            "allele": "A",
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
        "fixture": True,
        "note": (
            "Synthetic VEP component used only for integration "
            "plumbing. Not a biological annotation for this variant."
        ),
    }
]


# ------------------------------------------------------------
# Unified Evidence Bundle
#
# IMPORTANT:
# Use the real stable builder API.
# ------------------------------------------------------------

bundle = build_variant_evidence_bundle(
    assembly=ASSEMBLY,
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

assert bundle.variant_key == VARIANT_KEY


# ------------------------------------------------------------
# Stable Variant Evidence Contract
# ------------------------------------------------------------

contract = build_variant_evidence_contract(bundle)

assert contract.identity.variant_key == VARIANT_KEY


# ------------------------------------------------------------
# Interpretation Engine v0.1
# ------------------------------------------------------------

profile = interpret_variant(contract)
data = profile.to_dict()

assert data["variant_key"] == VARIANT_KEY

population = data["dimensions"]["population"]
clinical = data["dimensions"]["clinical"]
functional = data["dimensions"]["functional"]


# ------------------------------------------------------------
# Clinical evidence
# ------------------------------------------------------------

assert clinical["status"] == "AVAILABLE"

clinvar_significance = [
    observation
    for observation in clinical["observations"]
    if observation["code"]
    == "CLINVAR_SUBMITTED_SIGNIFICANCE"
]

assert clinvar_significance

assert any(
    "Uncertain_significance" in observation["value"]
    for observation in clinvar_significance
)

assert all(
    observation["details"]["platform_classification"] is False
    for observation in clinvar_significance
)


# ------------------------------------------------------------
# Functional evidence
# ------------------------------------------------------------

assert functional["status"] == "AVAILABLE"

functional_codes = {
    observation["code"]
    for observation in functional["observations"]
}

assert "VEP_GENE_CONTEXT" in functional_codes
assert "VEP_CONSEQUENCES" in functional_codes
assert "VEP_TRANSCRIPT" in functional_codes

consequence = next(
    observation
    for observation in functional["observations"]
    if observation["code"] == "VEP_CONSEQUENCES"
)

assert "intron_variant" in consequence["value"]


# ------------------------------------------------------------
# Population evidence semantics
# ------------------------------------------------------------

if gnomad_status == "NOT_FOUND":
    assert population["status"] == "UNAVAILABLE"

    observation = next(
        observation
        for observation in population["observations"]
        if observation["code"] == "GNOMAD_NOT_FOUND"
    )

    does_not_mean = observation["details"]["does_not_mean"]

    assert "AF=0" in does_not_mean
    assert "rare" in does_not_mean
    assert "pathogenic" in does_not_mean

elif gnomad_status == "ERROR":
    assert population["status"] == "ERROR"
    assert "population" in data["error_dimensions"]

elif gnomad_status == "FOUND":
    assert population["status"] in {
        "AVAILABLE",
        "PARTIAL",
    }


# ------------------------------------------------------------
# Safety boundary
# ------------------------------------------------------------

forbidden_keys = {
    "ranking_score",
    "score",
    "tier",
    "pathogenicity_probability",
    "acmg_classification",
    "diagnosis",
}

assert forbidden_keys.isdisjoint(data.keys())

for dimension_name in (
    "population",
    "clinical",
    "functional",
):
    assert forbidden_keys.isdisjoint(
        data["dimensions"][dimension_name].keys()
    )


# ------------------------------------------------------------
# Smoke summary
# ------------------------------------------------------------

summary = {
    "variant_key": VARIANT_KEY,
    "interpretation_schema_version": data["schema_version"],
    "dimensions": {
        "population": population["status"],
        "clinical": clinical["status"],
        "functional": functional["status"],
    },
    "available_dimensions": data["available_dimensions"],
    "error_dimensions": data["error_dimensions"],
    "source_statuses": {
        "ClinVar": "FOUND",
        "VEP": "FOUND",
        "gnomAD": gnomad_status,
    },
    "clinical": {
        "variation_ids": variation_ids,
        "submitted_significance": [
            observation["value"]
            for observation in clinvar_significance
        ],
    },
    "functional": {
        "genes": list(contract.functional.genes),
        "consequences": list(
            contract.functional.consequences
        ),
    },
    "safety_semantics": {
        "gnomad_not_found_means_af_zero": False,
        "gnomad_not_found_means_rare": False,
        "ranking_performed": False,
        "acmg_classification_performed": False,
        "pathogenicity_probability_assigned": False,
        "synthetic_vep_biological_interpretation_allowed": False,
    },
}

print(
    "PASS: INTERPRETATION ENGINE INTEGRATION SMOKE v0.1"
)

print(
    json.dumps(
        summary,
        indent=2,
    )
)
PY
