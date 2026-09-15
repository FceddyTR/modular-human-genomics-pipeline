import pytest

from genomics_platform.interpretation.case_model import (
    build_candidate_case,
)
from genomics_platform.interpretation.interpretation_engine import (
    interpret_variant,
)
from genomics_platform.phenotype.hpo_semantic_similarity import (
    HPOSemanticSimilarity,
)
from genomics_platform.ranking.components.semantic_phenotype_component import (
    semantic_phenotype_component,
)
from test_interpretation_engine import (
    make_contract,
)
from test_phenotype_interpreter import (
    make_bundle,
)


ONTOLOGY = (
    "data/evidence/hpo_ontology/"
    "hpo_ontology.sqlite"
)

IC = (
    "data/evidence/hpo_ontology/"
    "hpo_ic.sqlite"
)


@pytest.fixture
def semantic():
    return HPOSemanticSimilarity(
        ontology_database=ONTOLOGY,
        information_content_database=IC,
    )


def test_gene_mismatch_blocks_semantic(
    semantic,
):
    contract = make_contract()
    profile = interpret_variant(
        contract
    )

    candidate = build_candidate_case(
        candidate_id="mismatch",
        evidence_contract=contract,
        interpretation_profile=profile,
        gene_symbol="NCSTN",
        present_hpo_terms=[
            "HP:0040154",
        ],
    )

    result = (
        semantic_phenotype_component(
            candidate,
            20.0,
            semantic,
        )
    )

    assert (
        candidate.gene_consistency
        == "MISMATCH"
    )

    assert result.status == "BLOCKED"
    assert result.contribution == 0.0


def test_zero_maximum_disables_component(
    semantic,
):
    contract = make_contract()
    profile = interpret_variant(
        contract
    )

    candidate = build_candidate_case(
        candidate_id="disabled",
        evidence_contract=contract,
        interpretation_profile=profile,
    )

    result = (
        semantic_phenotype_component(
            candidate,
            0.0,
            semantic,
        )
    )

    assert result.contribution == 0.0
    assert result.max_contribution == 0.0


def test_no_patient_hpo_is_unavailable(
    semantic,
):
    contract = make_contract()
    profile = interpret_variant(
        contract
    )

    candidate = build_candidate_case(
        candidate_id="no-hpo",
        evidence_contract=contract,
        interpretation_profile=profile,
        gene_symbol="OR4F5",
    )

    result = (
        semantic_phenotype_component(
            candidate,
            20.0,
            semantic,
        )
    )

    assert (
        candidate.gene_consistency
        == "MATCH"
    )

    assert result.status == "UNAVAILABLE"
    assert result.contribution == 0.0


def test_semantic_boundaries_in_limitations(
    semantic,
):
    contract = make_contract()

    profile = interpret_variant(
        contract,
        gene_evidence_context=make_bundle(),
        present_hpo_terms=[
            "HP:0040154",
        ],
    )

    candidate = build_candidate_case(
        candidate_id="boundaries",
        evidence_contract=contract,
        interpretation_profile=profile,
        gene_symbol="OR4F5",
        present_hpo_terms=[
            "HP:0040154",
        ],
    )

    result = (
        semantic_phenotype_component(
            candidate,
            20.0,
            semantic,
        )
    )

    assert (
        candidate.gene_consistency
        == "MATCH"
    )

    assert result.status == "SUPPORTING"

    assert (
        0.0
        < result.contribution
        <= 20.0
    )

    text = " ".join(
        result.limitations
    ).lower()

    assert "probability" in text
    assert "diagnosis" in text
    assert "pathogenicity" in text

