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
from genomics_platform.interpretation.case_model import (
    build_candidate_case,
)

from test_interpretation_engine import make_contract


GENE = "NCSTN"

PRESENT_HPO = [
    "HP:0040154",
    "HP:0000987",
]

PATHS = {
    "gencc_index":
        "data/evidence/gencc/gencc.index.json",
    "panelapp_json":
        "data/evidence/panelapp/NCSTN.panelapp.json",
    "hpo_database":
        "data/evidence/hpo/hpo.sqlite",
    "mondo_database":
        "data/evidence/mondo/mondo.sqlite",
    "orphanet_database":
        "results/evidence/orphanet/orphanet.sqlite",
}

for name, path in PATHS.items():
    if not Path(path).exists():
        raise SystemExit(
            f"FAIL: missing artifact: "
            f"{name} -> {path}"
        )


# ------------------------------------------------------------
# 1. Real NCSTN evidence
# ------------------------------------------------------------

gene_context = build_bundle(
    gene_symbol=GENE,
    **PATHS,
)

assert (
    gene_context["summary"]["record_count"]
    == 11
)

assert (
    gene_context["summary"]["source_counts"]
    == {
        "GenCC": 3,
        "HPO": 3,
        "Orphanet": 1,
        "PanelApp": 4,
    }
)


# ------------------------------------------------------------
# 2. Variant contract
#
# Synthetic contract is intentionally used here only for
# plumbing. Gene evidence itself is real.
# ------------------------------------------------------------

contract = make_contract()


# ------------------------------------------------------------
# 3. Interpretation Engine v0.4
# ------------------------------------------------------------

profile = interpret_variant(
    contract,
    gene_evidence_context=gene_context,
    present_hpo_terms=PRESENT_HPO,
    genotype_state="heterozygous",
    proband_sex="female",
)

assert profile.gene_disease is not None
assert profile.phenotype is not None
assert profile.inheritance is not None

assert (
    profile.gene_disease.status
    == "PARTIAL"
)

assert (
    profile.phenotype.status
    == "AVAILABLE"
)

assert (
    profile.inheritance.status
    == "AVAILABLE"
)


# ------------------------------------------------------------
# 4. Build CandidateCase
# ------------------------------------------------------------

candidate = build_candidate_case(
    candidate_id="NCSTN-demo-001",
    evidence_contract=contract,
    interpretation_profile=profile,
    gene_symbol=GENE,
    present_hpo_terms=PRESENT_HPO,
    genotype_state="heterozygous",
    proband_sex="female",
    family_context=(
        "Synthetic case context for "
        "integration testing only."
    ),
)

data = candidate.to_dict()


# ------------------------------------------------------------
# 5. Identity integrity
# ------------------------------------------------------------

assert (
    candidate.variant_key
    == contract.identity.variant_key
)

assert (
    data["variant_key"]
    == profile.variant_key
)

assert data["gene_symbol"] == GENE


# ------------------------------------------------------------
# 6. Case context integrity
# ------------------------------------------------------------

assert (
    data["case_context"][
        "present_hpo_terms"
    ]
    == PRESENT_HPO
)

assert (
    data["case_context"][
        "genotype_state"
    ]
    == "heterozygous"
)

assert (
    data["case_context"][
        "proband_sex"
    ]
    == "female"
)


# ------------------------------------------------------------
# 7. Six interpretation dimensions
# ------------------------------------------------------------

expected_dimensions = {
    "population",
    "clinical",
    "functional",
    "phenotype",
    "inheritance",
    "gene_disease",
}

available = set(
    data["summary"][
        "available_dimensions"
    ]
)

assert available == expected_dimensions

assert (
    data["summary"]["error_dimensions"]
    == []
)


# ------------------------------------------------------------
# 8. Serialized interpretation survives
# ------------------------------------------------------------

dimensions = (
    data["interpretation_profile"]
    ["dimensions"]
)

assert (
    dimensions["gene_disease"]["status"]
    == "PARTIAL"
)

assert (
    dimensions["phenotype"]["status"]
    == "AVAILABLE"
)

assert (
    dimensions["inheritance"]["status"]
    == "AVAILABLE"
)


# ------------------------------------------------------------
# 9. Semantic safety boundary
# ------------------------------------------------------------

boundaries = data[
    "semantic_boundaries"
]

assert (
    boundaries["candidate_ranking"]
    is False
)

assert (
    boundaries["candidate_tier"]
    is False
)

assert (
    boundaries["variant_classification"]
    is False
)

assert (
    boundaries[
        "acmg_amp_classification"
    ]
    is False
)

assert (
    boundaries[
        "pathogenicity_probability"
    ]
    is False
)

assert boundaries["diagnosis"] is False


# ------------------------------------------------------------
# 10. JSON serialization
# ------------------------------------------------------------

serialized = json.dumps(
    data,
    indent=2,
)

assert "NCSTN-demo-001" in serialized
assert "NCSTN" in serialized
assert "HP:0040154" in serialized


# ------------------------------------------------------------
# Result
# ------------------------------------------------------------

result = {
    "case_model_version": "0.1.0",
    "candidate_id": (
        candidate.candidate_id
    ),
    "variant_key": (
        candidate.variant_key
    ),
    "gene": candidate.gene_symbol,
    "real_gene_evidence": {
        "record_count": (
            gene_context[
                "summary"
            ]["record_count"]
        ),
        "source_counts": (
            gene_context[
                "summary"
            ]["source_counts"]
        ),
    },
    "case_context": (
        data["case_context"]
    ),
    "interpretation": {
        "available_dimensions": (
            data["summary"][
                "available_dimensions"
            ]
        ),
        "error_dimensions": (
            data["summary"][
                "error_dimensions"
            ]
        ),
        "gene_disease_status": (
            dimensions[
                "gene_disease"
            ]["status"]
        ),
        "phenotype_status": (
            dimensions[
                "phenotype"
            ]["status"]
        ),
        "inheritance_status": (
            dimensions[
                "inheritance"
            ]["status"]
        ),
    },
    "semantic_boundaries": boundaries,
}

print(
    "PASS: REAL NCSTN CASE MODEL "
    "SMOKE v0.1"
)

print(
    json.dumps(
        result,
        indent=2,
    )
)
PY
