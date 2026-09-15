#!/usr/bin/env bash
set -euo pipefail

export PYTHONPATH="${PYTHONPATH:-}:src:tests"

python - <<'PY'
import json
from pathlib import Path

from genomics_platform.evidence.gene_evidence_bundle import build_bundle
from genomics_platform.interpretation.interpretation_engine import interpret_variant
from test_interpretation_engine import make_contract


GENE = "NCSTN"

PATHS = {
    "gencc_index": "data/evidence/gencc/gencc.index.json",
    "panelapp_json": "data/evidence/panelapp/NCSTN.panelapp.json",
    "hpo_database": "data/evidence/hpo/hpo.sqlite",
    "mondo_database": "data/evidence/mondo/mondo.sqlite",
    "orphanet_database": "results/evidence/orphanet/orphanet.sqlite",
}

for name, path in PATHS.items():
    if not Path(path).exists():
        raise SystemExit(
            f"FAIL: missing artifact: {name} -> {path}"
        )


# ------------------------------------------------------------
# Real NCSTN evidence bundle
# ------------------------------------------------------------

gene_context = build_bundle(
    gene_symbol=GENE,
    **PATHS,
)

assert gene_context["summary"]["record_count"] == 11


# ------------------------------------------------------------
# Interpretation Engine v0.4
#
# Synthetic variant contract + real NCSTN gene evidence.
# Heterozygous genotype is used only to exercise inheritance
# compatibility plumbing.
# ------------------------------------------------------------

profile = interpret_variant(
    make_contract(),
    gene_evidence_context=gene_context,
    present_hpo_terms=[
        "HP:0040154",
        "HP:0000987",
    ],
    genotype_state="heterozygous",
)

data = profile.to_dict()

gene_disease = data["dimensions"]["gene_disease"]
phenotype = data["dimensions"]["phenotype"]
inheritance = data["dimensions"]["inheritance"]

assert gene_disease["status"] == "PARTIAL"
assert phenotype["status"] == "AVAILABLE"
assert inheritance["status"] == "AVAILABLE"


# ------------------------------------------------------------
# Extract inheritance source observations
# ------------------------------------------------------------

records = [
    observation
    for observation in inheritance["observations"]
    if observation["code"]
    == "INHERITANCE_SOURCE_EVIDENCE"
]

assert records


# ------------------------------------------------------------
# Real source coverage
# ------------------------------------------------------------

sources = {
    record["value"]["source"]
    for record in records
}

assert "GenCC" in sources
assert "PanelApp" in sources
assert "HPO" in sources


# ------------------------------------------------------------
# Real inheritance models
# ------------------------------------------------------------

models = {
    record["value"]["normalized_inheritance"]
    for record in records
}

assert "AUTOSOMAL_DOMINANT" in models

assert (
    "AUTOSOMAL_DOMINANT_OR_MONOALLELIC"
    in models
)


# ------------------------------------------------------------
# HPO HP:0000006 must survive as inheritance evidence
# ------------------------------------------------------------

hpo_records = [
    record
    for record in records
    if record["value"]["source"] == "HPO"
]

assert any(
    record["value"]["raw_inheritance"]
    == "Autosomal dominant inheritance"
    for record in hpo_records
)


# ------------------------------------------------------------
# Heterozygous compatibility
# ------------------------------------------------------------

compatible = [
    record
    for record in records
    if record["value"]["genotype_compatibility"]
    == "COMPATIBLE"
]

assert compatible

assert any(
    record["value"]["normalized_inheritance"]
    == "AUTOSOMAL_DOMINANT"
    for record in compatible
)


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary = next(
    observation
    for observation in inheritance["observations"]
    if observation["code"]
    == "INHERITANCE_SUMMARY"
)

assert "COMPATIBLE" in (
    summary["value"]["compatibility_states"]
)


# ------------------------------------------------------------
# All major interpretation dimensions coexist
# ------------------------------------------------------------

for dimension in (
    "population",
    "clinical",
    "functional",
    "gene_disease",
    "phenotype",
    "inheritance",
):
    assert dimension in data["available_dimensions"]

assert data["error_dimensions"] == []


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
assert forbidden.isdisjoint(inheritance.keys())


# ------------------------------------------------------------
# Human-readable result
# ------------------------------------------------------------

result = {
    "engine_version": "0.4.0",
    "gene": GENE,
    "case_context": {
        "genotype_state": "heterozygous",
        "present_hpo_terms": [
            "HP:0040154",
            "HP:0000987",
        ],
    },
    "real_gene_evidence": {
        "record_count": (
            gene_context["summary"]["record_count"]
        ),
        "source_counts": (
            gene_context["summary"]["source_counts"]
        ),
    },
    "interpretation_profile": {
        "gene_disease_status": (
            gene_disease["status"]
        ),
        "phenotype_status": (
            phenotype["status"]
        ),
        "inheritance_status": (
            inheritance["status"]
        ),
        "inheritance_source_count": len(records),
        "inheritance_sources": sorted(sources),
        "normalized_models": sorted(models),
        "compatible_evidence_count": len(
            compatible
        ),
        "available_dimensions": (
            data["available_dimensions"]
        ),
        "error_dimensions": (
            data["error_dimensions"]
        ),
    },
    "semantic_boundaries": {
        "compatibility_means_pathogenicity": False,
        "candidate_ranking_performed": False,
        "variant_classification_performed": False,
        "diagnosis_assigned": False,
        "phase_inferred": False,
        "compound_heterozygosity_inferred": False,
    },
}

print(
    "PASS: INTERPRETATION ENGINE + REAL "
    "INHERITANCE CONTEXT SMOKE v0.4"
)

print(
    json.dumps(
        result,
        indent=2,
    )
)
PY
