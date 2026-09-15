"""Controlled benchmark for ontology-aware phenotype ranking.

The cohort is constructed independently of the benchmark truth.

Three ranking modes are compared:

1. variant-first
2. exact phenotype
3. ontology-aware semantic phenotype

The causal candidate intentionally has no exact phenotype match,
but has an ontology-related HPO association.

Truth is introduced only after all ranking runs have completed.
"""

from dataclasses import replace

import pytest

from genomics_platform.interpretation.case_model import (
    build_candidate_case,
)
from genomics_platform.interpretation.interpretation_engine import (
    interpret_variant,
)
from genomics_platform.interpretation.interpretation_profile import (
    EvidenceDimension,
    InterpretationObservation,
)
from genomics_platform.phenotype.hpo_semantic_similarity import (
    HPOSemanticSimilarity,
)
from genomics_platform.ranking.ranking_engine import (
    rank_candidates,
)
from genomics_platform.ranking.ranking_modes import (
    RankingMode,
)
from genomics_platform.ranking.benchmark.truth_set import (
    RankingTruth,
)
from genomics_platform.ranking.benchmark.metrics import (
    evaluate_ranking,
)

from test_interpretation_engine import (
    make_contract,
)


PATIENT_HPO = "HP:0040154"

# Ontology-related to PATIENT_HPO.
CAUSAL_ASSOCIATED_HPO = "HP:0001061"

# Much weaker semantic relationship.
DECOY_ASSOCIATED_HPO = "HP:0000987"


@pytest.fixture
def semantic():
    return HPOSemanticSimilarity(
        ontology_database=(
            "data/evidence/hpo_ontology/"
            "hpo_ontology.sqlite"
        ),
        information_content_database=(
            "data/evidence/hpo_ontology/"
            "hpo_ic.sqlite"
        ),
    )


def phenotype_dimension(
    associated_hpo,
):
    """Create phenotype evidence with no exact patient match."""

    return EvidenceDimension(
        name="phenotype",
        status="PARTIAL",
        observations=(
            InterpretationObservation(
                code=(
                    "PHENOTYPE_EXACT_MATCH_SUMMARY"
                ),
                label=(
                    "Controlled exact phenotype "
                    "summary"
                ),
                value={
                    "matched_present": [],
                    "unmatched_present": [
                        PATIENT_HPO
                    ],
                    "matched_absent": [],
                    "unmatched_absent": [],
                },
                source="controlled_fixture",
            ),
            InterpretationObservation(
                code=(
                    "PHENOTYPE_SOURCE_ASSOCIATION"
                ),
                label=(
                    "Controlled HPO source "
                    "association"
                ),
                value={
                    "source": "HPO",
                    "gene_symbol": None,
                    "disease_id": (
                        "SYNTHETIC:DISEASE"
                    ),
                    "phenotype_ids": [
                        associated_hpo
                    ],
                    "present_exact_matches": [],
                    "absent_exact_matches": [],
                },
                source="HPO",
            ),
        ),
        limitations=(
            "Synthetic controlled semantic "
            "ranking fixture.",
        ),
    )


def make_contract_for_gene(
    gene,
    position,
    consequence,
    impact,
    significance,
    af,
):
    base = make_contract()

    identity = replace(
        base.identity,
        pos=position,
    )

    transcript = replace(
        base.functional.transcripts[0],
        gene_symbol=gene,
        gene_id=f"SYNTHETIC:{gene}",
        consequences=(consequence,),
        impact=impact,
    )

    functional = replace(
        base.functional,
        genes=(gene,),
        gene_ids=(f"SYNTHETIC:{gene}",),
        consequences=(consequence,),
        impacts=(impact,),
        transcripts=(transcript,),
    )

    clinical_record = replace(
        base.clinical.records[0],
        variation_id=f"SYNTHETIC-{gene}",
        clinical_significance=(
            significance,
        ),
    )

    clinical = replace(
        base.clinical,
        records=(clinical_record,),
    )

    exome = replace(
        base.population.exome,
        af=af,
        ac=1,
        an=100000,
    )

    population = replace(
        base.population,
        exome=exome,
        genome=None,
    )

    return replace(
        base,
        identity=identity,
        functional=functional,
        clinical=clinical,
        population=population,
    )


def make_candidate(
    *,
    candidate_id,
    gene,
    position,
    consequence,
    impact,
    significance,
    af,
    associated_hpo,
):
    contract = make_contract_for_gene(
        gene=gene,
        position=position,
        consequence=consequence,
        impact=impact,
        significance=significance,
        af=af,
    )

    profile = interpret_variant(
        contract
    )

    profile = replace(
        profile,
        phenotype=phenotype_dimension(
            associated_hpo
        ),
    )

    candidate = build_candidate_case(
        candidate_id=candidate_id,
        evidence_contract=contract,
        interpretation_profile=profile,
        gene_symbol=gene,
        present_hpo_terms=(
            PATIENT_HPO,
        ),
    )

    assert (
        candidate.gene_consistency
        == "MATCH"
    )

    return candidate


def build_cohort():
    # Strong variant-level decoy.
    strong_decoy = make_candidate(
        candidate_id="strong-decoy",
        gene="DECOY1",
        position=66927,
        consequence="frameshift_variant",
        impact="HIGH",
        significance="Pathogenic",
        af=1e-6,
        associated_hpo=(
            DECOY_ASSOCIATED_HPO
        ),
    )

    # Benchmark truth.
    #
    # Weaker variant evidence than strong-decoy,
    # but phenotype is ontology-near PATIENT_HPO.
    causal = make_candidate(
        candidate_id="causal",
        gene="CAUSAL1",
        position=66928,
        consequence="missense_variant",
        impact="MODERATE",
        significance="Uncertain_significance",
        af=1e-5,
        associated_hpo=(
            CAUSAL_ASSOCIATED_HPO
        ),
    )

    # Weak variant and weak phenotype decoy.
    weak_decoy = make_candidate(
        candidate_id="weak-decoy",
        gene="DECOY2",
        position=66929,
        consequence="intron_variant",
        impact="MODIFIER",
        significance="Uncertain_significance",
        af=0.005,
        associated_hpo=(
            DECOY_ASSOCIATED_HPO
        ),
    )

    return (
        strong_decoy,
        causal,
        weak_decoy,
    )


def component(
    result_candidate,
    name,
):
    return next(
        item
        for item
        in result_candidate.components
        if item.name == name
    )


def ids(result):
    return [
        candidate.candidate_id
        for candidate
        in result.candidates
    ]


def test_exact_mode_cannot_use_ontology_relationship():
    result = rank_candidates(
        build_cohort(),
        mode=RankingMode.PHENOTYPE,
    )

    by_id = {
        candidate.candidate_id:
        candidate
        for candidate
        in result.candidates
    }

    causal_phenotype = component(
        by_id["causal"],
        "phenotype",
    )

    assert (
        causal_phenotype.contribution
        == 0.0
    )


def test_semantic_mode_detects_ontology_relationship(
    semantic,
):
    result = rank_candidates(
        build_cohort(),
        mode=(
            RankingMode.PHENOTYPE_SEMANTIC
        ),
        semantic_engine=semantic,
    )

    by_id = {
        candidate.candidate_id:
        candidate
        for candidate
        in result.candidates
    }

    causal = component(
        by_id["causal"],
        "phenotype",
    )

    strong_decoy = component(
        by_id["strong-decoy"],
        "phenotype",
    )

    assert causal.status == "SUPPORTING"

    assert (
        causal.contribution
        > strong_decoy.contribution
    )

    assert (
        0.0
        < causal.contribution
        <= 20.0
    )


def test_three_mode_benchmark(
    semantic,
):
    cohort = build_cohort()

    variant_first = rank_candidates(
        cohort,
        mode=RankingMode.VARIANT_FIRST,
    )

    exact = rank_candidates(
        cohort,
        mode=RankingMode.PHENOTYPE,
    )

    semantic_result = rank_candidates(
        cohort,
        mode=(
            RankingMode.PHENOTYPE_SEMANTIC
        ),
        semantic_engine=semantic,
    )

    # IMPORTANT:
    # Truth is constructed only after all
    # ranking runs have completed.
    truth = RankingTruth(
        case_id=(
            "semantic-controlled-case-001"
        ),
        causal_candidate_ids=(
            "causal",
        ),
    )

    variant_metrics = evaluate_ranking(
        variant_first,
        truth,
        k_values=(1, 3),
    )

    exact_metrics = evaluate_ranking(
        exact,
        truth,
        k_values=(1, 3),
    )

    semantic_metrics = evaluate_ranking(
        semantic_result,
        truth,
        k_values=(1, 3),
    )

    print()
    print("=== VARIANT FIRST ===")
    for candidate in variant_first.candidates:
        print(
            candidate.candidate_id,
            candidate.prioritization_score,
            component(
                candidate,
                "phenotype",
            ).contribution,
        )

    print()
    print("=== EXACT PHENOTYPE ===")
    for candidate in exact.candidates:
        print(
            candidate.candidate_id,
            candidate.prioritization_score,
            component(
                candidate,
                "phenotype",
            ).contribution,
        )

    print()
    print("=== SEMANTIC PHENOTYPE ===")
    for candidate in semantic_result.candidates:
        print(
            candidate.candidate_id,
            candidate.prioritization_score,
            component(
                candidate,
                "phenotype",
            ).contribution,
        )

    print()
    print("=== CAUSAL RANKS ===")
    print(
        "variant_first:",
        variant_metrics.best_causal_rank,
    )
    print(
        "exact:",
        exact_metrics.best_causal_rank,
    )
    print(
        "semantic:",
        semantic_metrics.best_causal_rank,
    )

    print()
    print("=== RECIPROCAL RANK ===")
    print(
        "variant_first:",
        variant_metrics.reciprocal_rank,
    )
    print(
        "exact:",
        exact_metrics.reciprocal_rank,
    )
    print(
        "semantic:",
        semantic_metrics.reciprocal_rank,
    )

    # Semantic phenotype must recover information
    # that exact-HPO matching cannot use.
    #
    # Rank improvement is NOT required here:
    # stronger non-phenotype evidence may still keep
    # a decoy above the causal candidate.

    exact_by_id = {
        candidate.candidate_id: candidate
        for candidate in exact.candidates
    }

    semantic_by_id = {
        candidate.candidate_id: candidate
        for candidate in semantic_result.candidates
    }

    exact_causal_phenotype = component(
        exact_by_id["causal"],
        "phenotype",
    )

    semantic_causal_phenotype = component(
        semantic_by_id["causal"],
        "phenotype",
    )

    assert (
        semantic_causal_phenotype.contribution
        > exact_causal_phenotype.contribution
    )

    exact_gap = (
        exact_by_id["strong-decoy"]
        .prioritization_score
        - exact_by_id["causal"]
        .prioritization_score
    )

    semantic_gap = (
        semantic_by_id["strong-decoy"]
        .prioritization_score
        - semantic_by_id["causal"]
        .prioritization_score
    )

    assert semantic_gap < exact_gap

    # Benchmark truth remains post-ranking only.
    assert (
        semantic_metrics.best_causal_rank
        is not None
    )

    assert (
        semantic_metrics.hit_at_k[3]
        is True
    )


def test_semantic_truth_not_exposed(
    semantic,
):
    result = rank_candidates(
        build_cohort(),
        mode="phenotype_semantic",
        semantic_engine=semantic,
    )

    data = result.to_dict()

    assert "truth" not in data

    assert (
        data["semantic_boundaries"][
            "pathogenicity_probability"
        ]
        is False
    )

    assert (
        data["semantic_boundaries"][
            "diagnosis"
        ]
        is False
    )


def test_semantic_rank_recovery_when_variant_evidence_is_balanced(
    semantic,
):
    """Semantic phenotype can recover causal rank.

    This controlled case differs from the information-recovery
    benchmark above.

    The decoy has a modest variant-evidence advantage:
    VUS + frameshift + very low AF.

    The causal candidate has:
    VUS + missense + very low AF,
    but a substantially stronger ontology relationship to the
    patient's phenotype.

    Truth is introduced only after both rankings are complete.
    """

    decoy = make_candidate(
        candidate_id="rank-recovery-decoy",
        gene="DECOY_RECOVERY",
        position=67101,
        consequence="frameshift_variant",
        impact="HIGH",
        significance="Uncertain_significance",
        af=1e-6,
        associated_hpo=(
            DECOY_ASSOCIATED_HPO
        ),
    )

    causal = make_candidate(
        candidate_id="rank-recovery-causal",
        gene="CAUSAL_RECOVERY",
        position=67102,
        consequence="missense_variant",
        impact="MODERATE",
        significance="Uncertain_significance",
        af=1e-5,
        associated_hpo=(
            CAUSAL_ASSOCIATED_HPO
        ),
    )

    cohort = (
        decoy,
        causal,
    )

    exact = rank_candidates(
        cohort,
        mode=RankingMode.PHENOTYPE,
    )

    semantic_result = rank_candidates(
        cohort,
        mode=(
            RankingMode.PHENOTYPE_SEMANTIC
        ),
        semantic_engine=semantic,
    )

    # Truth remains outside the ranking engine
    # and is constructed only after ranking.
    truth = RankingTruth(
        case_id=(
            "semantic-rank-recovery-001"
        ),
        causal_candidate_ids=(
            "rank-recovery-causal",
        ),
    )

    exact_metrics = evaluate_ranking(
        exact,
        truth,
        k_values=(1, 2),
    )

    semantic_metrics = evaluate_ranking(
        semantic_result,
        truth,
        k_values=(1, 2),
    )

    exact_by_id = {
        candidate.candidate_id: candidate
        for candidate in exact.candidates
    }

    semantic_by_id = {
        candidate.candidate_id: candidate
        for candidate
        in semantic_result.candidates
    }

    exact_causal_phenotype = component(
        exact_by_id[
            "rank-recovery-causal"
        ],
        "phenotype",
    )

    semantic_causal_phenotype = component(
        semantic_by_id[
            "rank-recovery-causal"
        ],
        "phenotype",
    )

    semantic_decoy_phenotype = component(
        semantic_by_id[
            "rank-recovery-decoy"
        ],
        "phenotype",
    )

    # Exact matching cannot use the ontology-near term.
    assert (
        exact_causal_phenotype.contribution
        == 0.0
    )

    # Semantic similarity recovers phenotype information.
    assert (
        semantic_causal_phenotype.contribution
        > semantic_decoy_phenotype.contribution
    )

    # Observed controlled rank recovery:
    # exact #2 -> semantic #1.
    assert (
        exact_metrics.best_causal_rank
        == 2
    )

    assert (
        semantic_metrics.best_causal_rank
        == 1
    )

    assert (
        semantic_metrics.reciprocal_rank
        > exact_metrics.reciprocal_rank
    )

    assert (
        exact_metrics.hit_at_k[1]
        is False
    )

    assert (
        semantic_metrics.hit_at_k[1]
        is True
    )

    # The benchmark truth must never become
    # part of ranking output semantics.
    data = semantic_result.to_dict()

    assert "truth" not in data

    assert (
        data["semantic_boundaries"][
            "pathogenicity_probability"
        ]
        is False
    )

    assert (
        data["semantic_boundaries"][
            "diagnosis"
        ]
        is False
    )
