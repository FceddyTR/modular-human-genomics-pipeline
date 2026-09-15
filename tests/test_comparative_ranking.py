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
from genomics_platform.ranking.ranking_modes import (
    RankingMode,
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


def phenotype_dimension(
    matched,
    unmatched=(),
):
    return EvidenceDimension(
        name="phenotype",
        status="AVAILABLE",
        observations=(
            InterpretationObservation(
                code=(
                    "PHENOTYPE_EXACT_MATCH_SUMMARY"
                ),
                label=(
                    "Synthetic phenotype exact "
                    "match summary"
                ),
                value={
                    "matched_present": list(
                        matched
                    ),
                    "unmatched_present": list(
                        unmatched
                    ),
                    "matched_absent": [],
                    "unmatched_absent": [],
                },
                source=(
                    "synthetic_ranking_fixture"
                ),
            ),
        ),
        limitations=(
            "Synthetic ranking fixture.",
        ),
    )


def make_candidate(
    candidate_id,
    matched=(),
    unmatched=(),
):
    contract = make_contract()

    # Fixture contract is annotated to OR4F5.
    gene = "OR4F5"

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

    assert (
        candidate.gene_consistency
        == "MATCH"
    )

    return candidate


def ids(result):
    return [
        candidate.candidate_id
        for candidate
        in result.candidates
    ]


def test_variant_first_does_not_use_phenotype():
    causal = make_candidate(
        "z-causal",
        matched=(
            "HP:0040154",
            "HP:0000987",
        ),
    )

    decoy_a = make_candidate(
        "a-decoy",
        unmatched=(
            "HP:0040154",
            "HP:0000987",
        ),
    )

    result = rank_candidates(
        [
            causal,
            decoy_a,
        ],
        mode="variant_first",
    )

    # Same neutral baseline -> deterministic ID tie-break.
    assert ids(result) == [
        "a-decoy",
        "z-causal",
    ]

    assert (
        result.candidates[0]
        .prioritization_score
        == result.candidates[1]
        .prioritization_score
    )


def test_phenotype_mode_promotes_match():
    causal = make_candidate(
        "z-causal",
        matched=(
            "HP:0040154",
            "HP:0000987",
        ),
    )

    decoy_a = make_candidate(
        "a-decoy",
        unmatched=(
            "HP:0040154",
            "HP:0000987",
        ),
    )

    result = rank_candidates(
        [
            causal,
            decoy_a,
        ],
        mode="phenotype",
    )

    assert ids(result)[0] == "z-causal"

    assert (
        result.candidates[0]
        .prioritization_score
        > result.candidates[1]
        .prioritization_score
    )


def test_same_cohort_comparative_benchmark():
    causal = make_candidate(
        "z-causal",
        matched=(
            "HP:0040154",
            "HP:0000987",
        ),
    )

    decoys = [
        make_candidate(
            "a-decoy",
            unmatched=(
                "HP:0040154",
                "HP:0000987",
            ),
        ),
        make_candidate(
            "b-decoy",
            unmatched=(
                "HP:0040154",
                "HP:0000987",
            ),
        ),
        make_candidate(
            "c-decoy",
            unmatched=(
                "HP:0040154",
                "HP:0000987",
            ),
        ),
    ]

    cohort = [
        causal,
        *decoys,
    ]

    baseline = rank_candidates(
        cohort,
        mode=RankingMode.VARIANT_FIRST,
    )

    phenotype = rank_candidates(
        cohort,
        mode=RankingMode.PHENOTYPE,
    )

    truth = RankingTruth(
        case_id="synthetic-case-001",
        causal_candidate_ids=(
            "z-causal",
        ),
    )

    baseline_metrics = evaluate_ranking(
        baseline,
        truth,
        k_values=(
            1,
            3,
            5,
        ),
    )

    phenotype_metrics = evaluate_ranking(
        phenotype,
        truth,
        k_values=(
            1,
            3,
            5,
        ),
    )

    comparison = RankingComparison(
        baseline=baseline_metrics,
        phenotype=phenotype_metrics,
    )

    assert (
        baseline_metrics.best_causal_rank
        == 4
    )

    assert (
        phenotype_metrics.best_causal_rank
        == 1
    )

    assert (
        comparison.causal_rank_improvement
        == 3
    )

    assert (
        baseline_metrics.hit_at_k[1]
        is False
    )

    assert (
        phenotype_metrics.hit_at_k[1]
        is True
    )

    assert (
        phenotype_metrics.reciprocal_rank
        > baseline_metrics.reciprocal_rank
    )


def test_truth_never_passed_to_ranking():
    causal = make_candidate(
        "causal",
        matched=("HP:0040154",),
    )

    result = rank_candidates(
        [causal],
        mode="phenotype",
    )

    data = result.to_dict()

    assert "truth" not in data
    assert "causal" not in (
        data["semantic_boundaries"]
    )


def test_ranking_semantics():
    result = rank_candidates(
        [
            make_candidate(
                "candidate",
                matched=(
                    "HP:0040154",
                ),
            )
        ],
        mode="phenotype",
    )

    boundaries = (
        result.to_dict()[
            "semantic_boundaries"
        ]
    )

    assert (
        boundaries[
            "relative_prioritization_only"
        ]
        is True
    )

    assert (
        boundaries[
            "pathogenicity_probability"
        ]
        is False
    )

    assert boundaries["diagnosis"] is False


def test_gene_mismatch_cannot_receive_phenotype_boost():
    contract = make_contract()
    profile = interpret_variant(contract)

    profile = replace(
        profile,
        phenotype=phenotype_dimension(
            matched=(
                "HP:0040154",
                "HP:0000987",
            ),
        ),
    )

    candidate = build_candidate_case(
        candidate_id="mismatch-phenotype",
        evidence_contract=contract,
        interpretation_profile=profile,
        gene_symbol="NCSTN",
        present_hpo_terms=(
            "HP:0040154",
            "HP:0000987",
        ),
    )

    assert candidate.gene_consistency == "MISMATCH"

    result = rank_candidates(
        [candidate],
        mode="phenotype",
    )

    ranked = result.candidates[0]

    phenotype = next(
        component
        for component in ranked.components
        if component.name == "phenotype"
    )

    assert phenotype.status == "BLOCKED"
    assert phenotype.contribution == 0.0

    assert any(
        "DATA_INTEGRITY_WARNING" in warning
        for warning in ranked.warnings
    )
