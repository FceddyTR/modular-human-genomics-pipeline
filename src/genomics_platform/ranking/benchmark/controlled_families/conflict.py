"""Controlled conflicting-evidence ranking stress cases."""

from genomics_platform.ranking.benchmark.controlled_suite import (
    ControlledCase,
    ControlledCaseSpec,
)
from genomics_platform.ranking.benchmark.synthetic_factory import (
    DEFAULT_DISTANT_HPO,
    DEFAULT_NEAR_HPO,
    DEFAULT_PATIENT_HPO,
    make_gene_context,
    make_synthetic_candidate,
)


def _case(case_id, description, expected_behavior, candidates):
    return ControlledCase(
        spec=ControlledCaseSpec(
            case_id=case_id,
            family="conflict",
            description=description,
            causal_candidate_ids=("causal",),
            expected_behavior=expected_behavior,
            tags=("synthetic", "controlled", "stress", "conflict"),
        ),
        candidates=tuple(candidates),
    )


def phenotype_vs_variant():
    return _case(
        "conflict-phenotype-vs-variant-001",
        "Ontology-near phenotype competes with stronger variant evidence.",
        "Semantic mode should expose the phenotype/variant trade-off.",
        (
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="CON1D",
                position=78001,
                consequence="frameshift_variant",
                impact="HIGH",
                significance="Pathogenic",
                af=1e-6,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_DISTANT_HPO,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="CON1C",
                position=78002,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_NEAR_HPO,
            ),
        ),
    )


def rarity_vs_function():
    return _case(
        "conflict-rarity-vs-function-001",
        "Rarity and functional consequence favor opposite candidates.",
        "Ranking should expose rather than hide competing evidence.",
        (
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="CON2D",
                position=78011,
                consequence="intron_variant",
                impact="MODIFIER",
                af=1e-6,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="CON2C",
                position=78012,
                consequence="frameshift_variant",
                impact="HIGH",
                af=0.005,
            ),
        ),
    )


def clinical_vs_population():
    return _case(
        "conflict-clinical-vs-population-001",
        "Clinical source assertion and population frequency disagree.",
        "Clinical and population components should remain independently visible.",
        (
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="CON3D",
                position=78021,
                significance="Pathogenic",
                af=0.02,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="CON3C",
                position=78022,
                significance="Uncertain_significance",
                af=1e-6,
            ),
        ),
    )


def gene_context_vs_variant():
    context = make_gene_context(
        gene="CON4C",
        gencc_classification="Definitive",
        panelapp_confidence="Green",
    )

    return _case(
        "conflict-gene-vs-variant-001",
        "Strong gene context competes with stronger variant-level evidence.",
        "Different aggregation strategies may resolve the conflict differently.",
        (
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="CON4D",
                position=78031,
                consequence="frameshift_variant",
                impact="HIGH",
                significance="Pathogenic",
                af=1e-6,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="CON4C",
                position=78032,
                gene_context=context,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
            ),
        ),
    )


BUILDERS = (
    phenotype_vs_variant,
    rarity_vs_function,
    clinical_vs_population,
    gene_context_vs_variant,
)
