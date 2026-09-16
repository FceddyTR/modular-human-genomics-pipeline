"""Controlled missing-evidence ranking stress cases."""

from genomics_platform.ranking.benchmark.controlled_suite import (
    ControlledCase,
    ControlledCaseSpec,
)
from genomics_platform.ranking.benchmark.synthetic_factory import (
    make_synthetic_candidate,
)


def _case(case_id, description, expected_behavior, candidates):
    return ControlledCase(
        spec=ControlledCaseSpec(
            case_id=case_id,
            family="missingness",
            description=description,
            causal_candidate_ids=("causal",),
            expected_behavior=expected_behavior,
            tags=("synthetic", "controlled", "stress", "missingness"),
        ),
        candidates=tuple(candidates),
    )


def population_not_found_vs_rare():
    return _case(
        "missingness-not-found-vs-rare-001",
        "NOT_FOUND population evidence competes with measured rarity.",
        "NOT_FOUND must not be interpreted as AF=0.",
        (
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="MIS1D",
                position=77001,
                af=None,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="MIS1C",
                position=77002,
                af=1e-6,
            ),
        ),
    )


def population_error_vs_rare():
    return _case(
        "missingness-error-vs-rare-001",
        "Population ERROR competes with measured rarity.",
        "ERROR must abstain rather than provide rarity support.",
        (
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="MIS2D",
                position=77011,
                af=None,
                population_status="ERROR",
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="MIS2C",
                position=77012,
                af=1e-6,
            ),
        ),
    )


def sparse_candidate():
    return _case(
        "missingness-sparse-candidate-001",
        "Sparse candidate competes with a fully measured candidate.",
        "Missing evidence must not create artificial positive support.",
        (
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="MIS3D",
                position=77021,
                af=None,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="MIS3C",
                position=77022,
                af=1e-5,
            ),
        ),
    )


def error_with_strong_variant():
    return _case(
        "missingness-error-strong-variant-001",
        "Strong independent variant evidence coexists with population ERROR.",
        "Population ERROR should remain separate from independent evidence.",
        (
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="MIS4D",
                position=77031,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="MIS4C",
                position=77032,
                consequence="frameshift_variant",
                impact="HIGH",
                significance="Pathogenic",
                af=None,
                population_status="ERROR",
            ),
        ),
    )


BUILDERS = (
    population_not_found_vs_rare,
    population_error_vs_rare,
    sparse_candidate,
    error_with_strong_variant,
)
