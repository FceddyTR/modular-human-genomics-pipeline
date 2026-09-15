"""Controlled comparative ranking acceptance test.

The synthetic evidence profiles are defined independently of
the expected ranking result.

Truth is introduced only AFTER ranking has completed.
"""

from dataclasses import replace

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
from genomics_platform.ranking.ranking_engine import (
    rank_candidates,
)
from genomics_platform.ranking.benchmark.truth_set import (
    RankingTruth,
)
from genomics_platform.ranking.benchmark.metrics import (
    evaluate_ranking,
)
from genomics_platform.ranking.benchmark.comparison import (
    RankingComparison,
)

from test_interpretation_engine import make_contract


HPO_1 = "HP:0040154"
HPO_2 = "HP:0000987"


def phenotype_dimension(
    matched=(),
    unmatched=(),
):
    return EvidenceDimension(
        name="phenotype",
        status="AVAILABLE",
        observations=(
            InterpretationObservation(
                code="PHENOTYPE_EXACT_MATCH_SUMMARY",
                label="Controlled phenotype fixture",
                value={
                    "matched_present": list(matched),
                    "unmatched_present": list(unmatched),
                    "matched_absent": [],
                    "unmatched_absent": [],
                },
                source="controlled_fixture",
            ),
        ),
        limitations=(
            "Synthetic controlled ranking fixture.",
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
    candidate_id,
    gene,
    position,
    consequence,
    impact,
    significance,
    af,
    matched=(),
    unmatched=(),
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
            matched=matched,
            unmatched=unmatched,
        ),
    )

    candidate = build_candidate_case(
        candidate_id=candidate_id,
        evidence_contract=contract,
        interpretation_profile=profile,
        gene_symbol=gene,
        present_hpo_terms=(
            list(matched)
            + list(unmatched)
        ),
    )

    assert candidate.gene_consistency == "MATCH"

    return candidate


def build_cohort():
    # Strong variant-level decoy:
    # P/LP + very low AF + high-priority consequence,
    # but no phenotype match.
    strong_decoy = make_candidate(
        candidate_id="strong-decoy",
        gene="DECOY1",
        position=66927,
        consequence="frameshift_variant",
        impact="HIGH",
        significance="Pathogenic",
        af=1e-6,
        unmatched=(HPO_1, HPO_2),
    )

    # Known truth for benchmark:
    # VUS + very low AF + missense consequence,
    # but complete phenotype match.
    causal = make_candidate(
        candidate_id="causal",
        gene="CAUSAL1",
        position=66928,
        consequence="missense_variant",
        impact="MODERATE",
        significance="Uncertain_significance",
        af=1e-5,
        matched=(HPO_1, HPO_2),
    )

    # Weak decoy:
    # VUS + higher AF + intronic consequence
    # and no phenotype match.
    weak_decoy = make_candidate(
        candidate_id="weak-decoy",
        gene="DECOY2",
        position=66929,
        consequence="intron_variant",
        impact="MODIFIER",
        significance="Uncertain_significance",
        af=0.005,
        unmatched=(HPO_1, HPO_2),
    )

    return (
        strong_decoy,
        causal,
        weak_decoy,
    )


def component(result_candidate, name):
    return next(
        item
        for item in result_candidate.components
        if item.name == name
    )


def test_controlled_contracts_are_gene_consistent():
    for candidate in build_cohort():
        assert candidate.gene_consistency == "MATCH"


def test_variant_first_prefers_strong_variant_decoy():
    cohort = build_cohort()

    result = rank_candidates(
        cohort,
        mode="variant_first",
    )

    assert (
        result.candidates[0].candidate_id
        == "strong-decoy"
    )


def test_phenotype_is_disabled_in_variant_first():
    result = rank_candidates(
        build_cohort(),
        mode="variant_first",
    )

    for candidate in result.candidates:
        phenotype = component(
            candidate,
            "phenotype",
        )

        assert phenotype.max_contribution == 0.0
        assert phenotype.contribution == 0.0


def test_phenotype_mode_rewards_matching_candidate():
    result = rank_candidates(
        build_cohort(),
        mode="phenotype",
    )

    by_id = {
        candidate.candidate_id: candidate
        for candidate in result.candidates
    }

    causal = component(
        by_id["causal"],
        "phenotype",
    )

    strong_decoy = component(
        by_id["strong-decoy"],
        "phenotype",
    )

    assert causal.contribution == 20.0
    assert strong_decoy.contribution == 0.0


def test_benchmark_truth_is_post_ranking_only():
    cohort = build_cohort()

    baseline = rank_candidates(
        cohort,
        mode="variant_first",
    )

    phenotype = rank_candidates(
        cohort,
        mode="phenotype",
    )

    # Truth is intentionally constructed only
    # after both ranking runs are complete.
    truth = RankingTruth(
        case_id="controlled-case-001",
        causal_candidate_ids=(
            "causal",
        ),
    )

    baseline_metrics = evaluate_ranking(
        baseline,
        truth,
        k_values=(1, 3),
    )

    phenotype_metrics = evaluate_ranking(
        phenotype,
        truth,
        k_values=(1, 3),
    )

    comparison = RankingComparison(
        baseline=baseline_metrics,
        phenotype=phenotype_metrics,
    )

    assert (
        baseline_metrics.best_causal_rank
        is not None
    )

    assert (
        phenotype_metrics.best_causal_rank
        is not None
    )

    assert (
        comparison.causal_rank_improvement
        is not None
    )


def test_benchmark_reports_real_observed_delta():
    cohort = build_cohort()

    baseline = rank_candidates(
        cohort,
        mode="variant_first",
    )

    phenotype = rank_candidates(
        cohort,
        mode="phenotype",
    )

    truth = RankingTruth(
        case_id="controlled-case-001",
        causal_candidate_ids=(
            "causal",
        ),
    )

    baseline_metrics = evaluate_ranking(
        baseline,
        truth,
        k_values=(1, 3),
    )

    phenotype_metrics = evaluate_ranking(
        phenotype,
        truth,
        k_values=(1, 3),
    )

    comparison = RankingComparison(
        baseline=baseline_metrics,
        phenotype=phenotype_metrics,
    )

    # We do NOT hard-code that phenotype must
    # magically make the truth rank #1.
    #
    # We require the phenotype signal to improve
    # the known truth relative to baseline.
    assert (
        comparison.causal_rank_improvement
        > 0
    )

    assert (
        phenotype_metrics.reciprocal_rank
        > baseline_metrics.reciprocal_rank
    )


def test_scores_are_not_pathogenicity_probabilities():
    result = rank_candidates(
        build_cohort(),
        mode="phenotype",
    )

    boundaries = (
        result.to_dict()[
            "semantic_boundaries"
        ]
    )

    assert (
        boundaries[
            "pathogenicity_probability"
        ]
        is False
    )

    assert (
        boundaries[
            "acmg_amp_classification"
        ]
        is False
    )

    assert boundaries["diagnosis"] is False
