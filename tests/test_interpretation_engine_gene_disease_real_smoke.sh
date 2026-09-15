#!/usr/bin/env bash
set -euo pipefail

export PYTHONPATH="${PYTHONPATH:-}:src:tests"

python - <<'PY'
import json
from pathlib import Path

from genomics_platform.evidence.gene_evidence_bundle import (
    build_bundle,
)
from genomics_platform.interpretation.interpretation_engine import (
    interpret_variant,
)

from test_interpretation_engine import make_contract


GENE = "NCSTN"

PATHS = {
    "gencc_index": "data/evidence/gencc/gencc.index.json",
    "panelapp_json": "data/evidence/panelapp/NCSTN.panelapp.json",
    "hpo_database": "data/evidence/hpo/hpo.sqlite",
    "mondo_database": "data/evidence/mondo/mondo.sqlite",
    "orphanet_database": (
        "results/evidence/orphanet/orphanet.sqlite"
    ),
}


# ------------------------------------------------------------
# Verify real evidence artifacts
# ------------------------------------------------------------

for name, path in PATHS.items():
    if not Path(path).exists():
        raise SystemExit(
            f"FAIL: required artifact missing: {name} -> {path}"
        )


# ------------------------------------------------------------
# Build real NCSTN gene evidence context
# ------------------------------------------------------------

gene_context = build_bundle(
    gene_symbol=GENE,
    **PATHS,
)

assert gene_context["query"]["gene_symbol"] == GENE
assert gene_context["summary"]["record_count"] == 11

assert gene_context["summary"]["source_counts"] == {
    "GenCC": 3,
    "HPO": 3,
    "Orphanet": 1,
    "PanelApp": 4,
}


# ------------------------------------------------------------
# Interpretation Engine
#
# make_contract() is a synthetic variant contract used only to
# exercise engine plumbing. The real biological evidence being
# tested here is the NCSTN gene context.
# ------------------------------------------------------------

contract = make_contract()

profile = interpret_variant(
    contract,
    gene_evidence_context=gene_context,
)

data = profile.to_dict()

gene_disease = data["dimensions"]["gene_disease"]

assert gene_disease is not None
assert gene_disease["name"] == "gene_disease"

# Four NCSTN records have NO_DISEASE_ID.
# Therefore PARTIAL is the expected semantic state.
assert gene_disease["status"] == "PARTIAL"

assert "gene_disease" in data["available_dimensions"]


# ------------------------------------------------------------
# Verify NCSTN survived interpretation
# ------------------------------------------------------------

gene_observation = next(
    observation
    for observation in gene_disease["observations"]
    if observation["code"] == "GENE_DISEASE_GENE"
)

assert gene_observation["value"] == GENE


# ------------------------------------------------------------
# Verify all real evidence sources survived
# ------------------------------------------------------------

source_records = [
    observation
    for observation in gene_disease["observations"]
    if observation["code"]
    == "GENE_DISEASE_SOURCE_RECORD"
]

observed_sources = {
    observation["value"]["source"]
    for observation in source_records
}

assert {
    "GenCC",
    "HPO",
    "Orphanet",
    "PanelApp",
}.issubset(observed_sources)


# ------------------------------------------------------------
# Verify MONDO-normalized diseases survived
# ------------------------------------------------------------

normalized_observation = next(
    observation
    for observation in gene_disease["observations"]
    if observation["code"]
    == "GENE_DISEASE_NORMALIZED_DISEASES"
)

assert len(normalized_observation["value"]) == 2

assert all(
    disease["mondo_id"]
    for disease in normalized_observation["value"]
)


# ------------------------------------------------------------
# Verify incomplete normalization remains visible
# ------------------------------------------------------------

normalization_statuses = {
    observation["value"]["normalization_status"]
    for observation in source_records
}

assert "FOUND" in normalization_statuses
assert "NO_DISEASE_ID" in normalization_statuses


# ------------------------------------------------------------
# Safety boundary
# ------------------------------------------------------------

forbidden = {
    "score",
    "ranking_score",
    "tier",
    "acmg_classification",
    "pathogenicity_probability",
    "diagnosis",
}

assert forbidden.isdisjoint(data.keys())
assert forbidden.isdisjoint(gene_disease.keys())

for observation in source_records:
    assert (
        observation["details"]["variant_pathogenicity"]
        is False
    )

    assert (
        observation["details"]["candidate_ranking"]
        is False
    )


# ------------------------------------------------------------
# Human-readable summary
# ------------------------------------------------------------

summary = {
    "engine_version": "0.2.0",
    "gene": GENE,
    "real_gene_evidence": {
        "record_count": (
            gene_context["summary"]["record_count"]
        ),
        "source_counts": (
            gene_context["summary"]["source_counts"]
        ),
        "disease_normalization": (
            gene_context["summary"]["disease_normalization"]
        ),
        "normalized_disease_count": (
            gene_context["summary"]["normalized_disease_count"]
        ),
    },
    "interpretation_profile": {
        "gene_disease_status": gene_disease["status"],
        "available_dimensions": data["available_dimensions"],
        "error_dimensions": data["error_dimensions"],
        "observation_count": len(
            gene_disease["observations"]
        ),
    },
    "semantic_boundaries": {
        "gene_disease_evidence_is_variant_classification": False,
        "gencc_classification_is_acmg_classification": False,
        "panelapp_confidence_is_variant_pathogenicity": False,
        "mondo_normalization_is_causality": False,
        "candidate_ranking_performed": False,
    },
}

print(
    "PASS: INTERPRETATION ENGINE + REAL "
    "GENE-DISEASE CONTEXT SMOKE v0.2"
)

print(
    json.dumps(
        summary,
        indent=2,
    )
)
PY
