"""Controlled negative-phenotype stress cases.

Negative phenotype observations are contextual evidence.
They are not hard exclusion rules and do not establish
variant pathogenicity or diagnosis.
"""

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
            family="negative_phenotype",
            description=description,
            causal_candidate_ids=("causal",),
            expected_behavior=expected_behavior,
            tags=(
                "synthetic",
                "controlled",
                "stress",
                "negative_phenotype",
                "no_hard_exclusion",
            ),
        ),
        candidates=tuple(candidates),
    )


def absent_term_does_not_exclude():
    context = make_gene_context(
        gene="NEGC1",
        gencc_classification="Definitive",
        panelapp_confidence="Green",
    )

    return _case(
        "negative-phenotype-no-exclusion-001",
        "Causal candidate carries an explicitly absent phenotype term.",
        "Negative phenotype context must not hard-exclude the candidate.",
        (
            make_synthetic_candidate(
                candidate_id="causal",
                gene="NEGC1",
                position=80011,
                gene_context=context,
                absent_hpo_terms=(DEFAULT_NEAR_HPO,),
            ),
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="NEGD1",
                position=80012,
            ),
        ),
    )


def positive_and_negative_context():
    context = make_gene_context(
        gene="NEGC2",
        gencc_classification="Definitive",
        panelapp_confidence="Green",
    )

    return _case(
        "negative-phenotype-mixed-context-001",
        "Positive phenotype relevance coexists with an absent term.",
        "Positive and negative phenotype context must remain distinguishable.",
        (
            make_synthetic_candidate(
                candidate_id="causal",
                gene="NEGC2",
                position=80021,
                gene_context=context,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_NEAR_HPO,
                absent_hpo_terms=(DEFAULT_DISTANT_HPO,),
            ),
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="NEGD2",
                position=80022,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_DISTANT_HPO,
            ),
        ),
    )


def negative_with_variant_pressure():
    context = make_gene_context(
        gene="NEGC3",
        gencc_classification="Definitive",
        panelapp_confidence="Green",
    )

    return _case(
        "negative-phenotype-variant-pressure-001",
        "Negative phenotype context coexists with competing variant evidence.",
        "Negative phenotype must not become an absolute exclusion mechanism.",
        (
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="NEGD3",
                position=80031,
                consequence="frameshift_variant",
                impact="HIGH",
                significance="Pathogenic",
                af=1e-6,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="NEGC3",
                position=80032,
                gene_context=context,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
                absent_hpo_terms=(DEFAULT_NEAR_HPO,),
            ),
        ),
    )


def multiple_absent_terms():
    context = make_gene_context(
        gene="NEGC4",
        gencc_classification="Definitive",
        panelapp_confidence="Green",
    )

    return _case(
        "negative-phenotype-multiple-001",
        "Multiple explicitly absent HPO terms accompany the causal candidate.",
        "Multiple negative observations must remain contextual rather than binary exclusion.",
        (
            make_synthetic_candidate(
                candidate_id="causal",
                gene="NEGC4",
                position=80041,
                gene_context=context,
                absent_hpo_terms=(
                    DEFAULT_NEAR_HPO,
                    DEFAULT_DISTANT_HPO,
                ),
            ),
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="NEGD4",
                position=80042,
            ),
        ),
    )


BUILDERS = (
    absent_term_does_not_exclude,
    positive_and_negative_context,
    negative_with_variant_pressure,
    multiple_absent_terms,
)
