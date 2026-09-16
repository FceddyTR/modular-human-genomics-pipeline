"""Controlled semantic-phenotype ranking challenge cases.

These cases isolate ontology-aware phenotype behavior from ordinary
exact HPO matching. They are synthetic ranking diagnostics, not
clinical validation cases.
"""

from genomics_platform.ranking.benchmark.controlled_suite import (
    ControlledCase,
    ControlledCaseSpec,
)
from genomics_platform.ranking.benchmark.synthetic_factory import (
    DEFAULT_DISTANT_HPO,
    DEFAULT_NEAR_HPO,
    DEFAULT_PATIENT_HPO,
    make_synthetic_candidate,
)


def _case(
    case_id,
    description,
    expected_behavior,
    candidates,
):
    return ControlledCase(
        spec=ControlledCaseSpec(
            case_id=case_id,
            family="semantic_challenge",
            description=description,
            causal_candidate_ids=("causal",),
            expected_behavior=expected_behavior,
            tags=(
                "synthetic",
                "controlled",
                "semantic",
                "challenge",
            ),
        ),
        candidates=tuple(candidates),
    )


def semantic_near_vs_distant():
    """Isolate semantic discrimination without exact HPO support."""
    return _case(
        "semantic-near-vs-distant-001",
        (
            "Causal and decoy candidates have matched "
            "variant-level evidence, while only the causal "
            "phenotype association is ontology-near."
        ),
        (
            "Exact matching should not recover the ontology "
            "relationship; semantic mode should recover "
            "phenotype information for the causal candidate."
        ),
        (
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="SEM1D",
                position=80001,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_DISTANT_HPO,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="SEM1C",
                position=80002,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_NEAR_HPO,
            ),
        ),
    )


def semantic_rank_recovery():
    """Allow semantic phenotype evidence to overcome a small variant gap."""
    return _case(
        "semantic-rank-recovery-001",
        (
            "A functional-evidence decoy has a modest "
            "variant-level advantage over an ontology-near "
            "causal candidate."
        ),
        (
            "Exact phenotype ranking should retain the "
            "variant-level disadvantage, while semantic "
            "phenotype evidence should be capable of "
            "recovering the causal rank."
        ),
        (
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="SEM2D",
                position=80011,
                consequence="frameshift_variant",
                impact="HIGH",
                significance="Uncertain_significance",
                af=1e-5,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_DISTANT_HPO,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="SEM2C",
                position=80012,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_NEAR_HPO,
            ),
        ),
    )


def semantic_information_recovery():
    """Recover phenotype information without requiring a rank reversal."""
    return _case(
        "semantic-information-recovery-001",
        (
            "A strongly supported variant-level decoy competes "
            "with an ontology-near causal candidate."
        ),
        (
            "Semantic mode should recover phenotype information "
            "and reduce the evidence disadvantage without "
            "requiring the causal candidate to become rank 1."
        ),
        (
            make_synthetic_candidate(
                candidate_id="strong-decoy",
                gene="SEM3D",
                position=80021,
                consequence="frameshift_variant",
                impact="HIGH",
                significance="Pathogenic",
                af=1e-6,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_DISTANT_HPO,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="SEM3C",
                position=80022,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_NEAR_HPO,
            ),
        ),
    )


def semantic_exact_control():
    """Control for an already exact phenotype association."""
    return _case(
        "semantic-exact-control-001",
        (
            "The causal candidate already has an exact patient "
            "HPO association while the decoy is phenotype-distant."
        ),
        (
            "Ontology-aware phenotype processing should preserve "
            "the direct phenotype signal without creating an "
            "artificial reversal."
        ),
        (
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="SEM4D",
                position=80031,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_DISTANT_HPO,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="SEM4C",
                position=80032,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_PATIENT_HPO,
                exact_phenotype_match=True,
            ),
        ),
    )


BUILDERS = (
    semantic_near_vs_distant,
    semantic_rank_recovery,
    semantic_information_recovery,
    semantic_exact_control,
)
