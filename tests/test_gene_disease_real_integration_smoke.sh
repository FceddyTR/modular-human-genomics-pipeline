#!/usr/bin/env bash
set -euo pipefail

export PYTHONPATH="${PYTHONPATH:-}:src"

python - <<'PY'
import json
from pathlib import Path

from genomics_platform.evidence.gene_evidence_bundle import (
    build_bundle,
)
from genomics_platform.interpretation.gene_disease_interpreter import (
    interpret_gene_disease,
)


GENE = "NCSTN"

GENCC = "data/evidence/gencc/gencc.index.json"
PANELAPP = "data/evidence/panelapp/NCSTN.panelapp.json"
HPO = "data/evidence/hpo/hpo.sqlite"
MONDO = "data/evidence/mondo/mondo.sqlite"
ORPHANET = "results/evidence/orphanet/orphanet.sqlite"


# ------------------------------------------------------------
# Required local evidence artifacts
# ------------------------------------------------------------

paths = {
    "GenCC": GENCC,
    "PanelApp": PANELAPP,
    "HPO": HPO,
    "MONDO": MONDO,
    "Orphanet": ORPHANET,
}

for source, path in paths.items():
    if not Path(path).exists():
        raise SystemExit(
            f"FAIL: required {source} artifact missing: {path}"
        )


# ------------------------------------------------------------
# Real Gene Evidence Bundle
# ------------------------------------------------------------

bundle = build_bundle(
    gene_symbol=GENE,
    gencc_index=GENCC,
    panelapp_json=PANELAPP,
    hpo_database=HPO,
    mondo_database=MONDO,
    orphanet_database=ORPHANET,
)

assert bundle["query"]["gene_symbol"] == GENE

summary = bundle["summary"]

assert summary["record_count"] > 0
assert bundle["records"]


# ------------------------------------------------------------
# Verify real source aggregation
# ------------------------------------------------------------

source_counts = summary["source_counts"]

expected_sources = {
    "GenCC",
    "PanelApp",
    "HPO",
    "Orphanet",
}

missing_sources = (
    expected_sources - set(source_counts)
)

assert not missing_sources, (
    "Expected NCSTN evidence sources missing: "
    f"{sorted(missing_sources)}"
)


# ------------------------------------------------------------
# Verify bundle safety boundary
# ------------------------------------------------------------

scope = bundle["scope"]

assert scope["clinical_classification"] is False
assert scope["variant_pathogenicity"] is False
assert scope["candidate_ranking"] is False


# ------------------------------------------------------------
# Gene-Disease Interpreter
# ------------------------------------------------------------

dimension = interpret_gene_disease(bundle)
data = dimension.to_dict()

assert data["name"] == "gene_disease"

assert data["status"] in {
    "AVAILABLE",
    "PARTIAL",
}

assert data["observations"]


# ------------------------------------------------------------
# Gene identity preserved
# ------------------------------------------------------------

gene_observation = next(
    observation
    for observation in data["observations"]
    if observation["code"] == "GENE_DISEASE_GENE"
)

assert gene_observation["value"] == GENE


# ------------------------------------------------------------
# Source records preserved
# ------------------------------------------------------------

source_records = [
    observation
    for observation in data["observations"]
    if observation["code"]
    == "GENE_DISEASE_SOURCE_RECORD"
]

assert source_records

observed_sources = {
    observation["value"]["source"]
    for observation in source_records
}

assert expected_sources.issubset(
    observed_sources
)


# ------------------------------------------------------------
# MONDO normalization preserved
# ------------------------------------------------------------

normalized = bundle["normalized_diseases"]

assert normalized

for disease in normalized:
    assert disease["mondo_id"]
    assert disease["mondo_label"]


# ------------------------------------------------------------
# Critical semantic boundary
# ------------------------------------------------------------

for observation in source_records:
    details = observation["details"]

    assert (
        details["variant_pathogenicity"]
        is False
    )

    assert (
        details["candidate_ranking"]
        is False
    )


# ------------------------------------------------------------
# No ranking / classification outputs
# ------------------------------------------------------------

forbidden = {
    "score",
    "ranking_score",
    "tier",
    "acmg_classification",
    "pathogenicity_probability",
    "diagnosis",
}

assert forbidden.isdisjoint(
    data.keys()
)


# ------------------------------------------------------------
# Human-readable summary
# ------------------------------------------------------------

result = {
    "gene": GENE,
    "gene_evidence_bundle": {
        "schema_version": bundle["schema_version"],
        "engine": bundle["engine"],
        "record_count": summary["record_count"],
        "source_counts": source_counts,
        "disease_normalization": (
            summary["disease_normalization"]
        ),
        "normalized_disease_count": (
            summary["normalized_disease_count"]
        ),
    },
    "gene_disease_interpretation": {
        "status": data["status"],
        "observation_count": len(
            data["observations"]
        ),
    },
    "safety_semantics": {
        "variant_classification_performed": False,
        "candidate_ranking_performed": False,
        "pathogenicity_probability_assigned": False,
        "mondo_normalization_means_causality": False,
    },
}

print(
    "PASS: REAL GENE-DISEASE INTEGRATION SMOKE v0.1"
)

print(
    json.dumps(
        result,
        indent=2,
    )
)
PY
