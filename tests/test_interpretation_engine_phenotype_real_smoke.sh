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

PRESENT_HPO = [
    "HP:0040154",  # Acne inversa
    "HP:0000987",  # Atypical scarring of skin
]

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
# Required real evidence artifacts
# ------------------------------------------------------------

for name, path in PATHS.items():
    if not Path(path).exists():
        raise SystemExit(
            f"FAIL: missing artifact: {name} -> {path}"
        )


# ------------------------------------------------------------
# Build real NCSTN gene evidence
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
# Run Interpretation Engine v0.3
#
# Variant contract remains synthetic. This smoke tests the
# integration plumbing of real gene/phenotype evidence.
# ------------------------------------------------------------

profile = interpret_variant(
    make_contract(),
    gene_evidence_context=gene_context,
    present_hpo_terms=PRESENT_HPO,
)

data = profile.to_dict()

gene_disease = data["dimensions"]["gene_disease"]
phenotype = data["dimensions"]["phenotype"]


# ------------------------------------------------------------
# Gene-disease dimension
# ------------------------------------------------------------

assert gene_disease is not None
assert gene_disease["status"] == "PARTIAL"


# ------------------------------------------------------------
# Phenotype dimension
# ------------------------------------------------------------

assert phenotype is not None
assert phenotype["status"] == "AVAILABLE"

summary = next(
    observation
    for observation in phenotype["observations"]
    if observation["code"]
    == "PHENOTYPE_EXACT_MATCH_SUMMARY"
)

assert summary["value"]["matched_present"] == PRESENT_HPO
assert summary["value"]["unmatched_present"] == []
assert summary["value"]["matched_absent"] == []
assert summary["value"]["unmatched_absent"] == []


# ------------------------------------------------------------
# Both dimensions exposed through profile
# ------------------------------------------------------------

assert "gene_disease" in data["available_dimensions"]
assert "phenotype" in data["available_dimensions"]

assert data["error_dimensions"] == []


# ------------------------------------------------------------
# Verify real HPO source evidence survived
# ------------------------------------------------------------

source_associations = [
    observation
    for observation in phenotype["observations"]
    if observation["code"]
    == "PHENOTYPE_SOURCE_ASSOCIATION"
]

assert len(source_associations) == 3

observed_hpo_terms = {
    term
    for observation in source_associations
    for term in observation["value"]["phenotype_ids"]
}

assert "HP:0040154" in observed_hpo_terms
assert "HP:0000987" in observed_hpo_terms

# Inheritance-coded HPO evidence must remain preserved.
assert "HP:0000006" in observed_hpo_terms


# ------------------------------------------------------------
# Safety boundary
# ------------------------------------------------------------

forbidden = {
    "score",
    "ranking_score",
    "phenotype_score",
    "tier",
    "acmg_classification",
    "pathogenicity_probability",
    "diagnosis",
}

assert forbidden.isdisjoint(data.keys())
assert forbidden.isdisjoint(phenotype.keys())
assert forbidden.isdisjoint(gene_disease.keys())


# ------------------------------------------------------------
# Human-readable result
# ------------------------------------------------------------

result = {
    "engine_version": "0.3.0",
    "gene": GENE,
    "patient_phenotype": {
        "present_hpo_terms": PRESENT_HPO,
    },
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
        "phenotype_status": phenotype["status"],
        "matched_present_hpo": (
            summary["value"]["matched_present"]
        ),
        "available_dimensions": (
            data["available_dimensions"]
        ),
        "error_dimensions": (
            data["error_dimensions"]
        ),
    },
    "semantic_boundaries": {
        "phenotype_similarity_score_assigned": False,
        "candidate_ranking_performed": False,
        "variant_classification_performed": False,
        "diagnosis_assigned": False,
        "inheritance_interpretation_performed": False,
    },
}

print(
    "PASS: INTERPRETATION ENGINE + REAL "
    "PHENOTYPE CONTEXT SMOKE v0.3"
)

print(
    json.dumps(
        result,
        indent=2,
    )
)
PY
