"""Initial real controlled benchmark cases."""

from __future__ import annotations

from genomics_platform.ranking.benchmark.controlled_suite import (
    ControlledCase,
    ControlledCaseSpec,
    ControlledSuiteRegistry,
)
from genomics_platform.ranking.benchmark.synthetic_factory import (
    DEFAULT_DISTANT_HPO,
    DEFAULT_NEAR_HPO,
    DEFAULT_PATIENT_HPO,
    make_synthetic_candidate,
)


def _case(
    *,
    case_id: str,
    family: str,
    description: str,
    expected_behavior: str,
    candidates,
    causal_id: str = "causal",
) -> ControlledCase:
    return ControlledCase(
        spec=ControlledCaseSpec(
            case_id=case_id,
            family=family,
            description=description,
            causal_candidate_ids=(
                causal_id,
            ),
            expected_behavior=(
                expected_behavior
            ),
            tags=(
                "synthetic",
                "controlled",
                family,
            ),
        ),
        candidates=tuple(candidates),
    )


def phenotype_semantic_near():
    return _case(
        case_id="phenotype-semantic-near-001",
        family="phenotype",
        description=(
            "Causal candidate has an ontology-near "
            "phenotype without an exact HPO match."
        ),
        expected_behavior=(
            "Semantic phenotype mode should recover "
            "information unavailable to exact matching."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="PDEC1",
                position=71001,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_DISTANT_HPO,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="PCAU1",
                position=71002,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_NEAR_HPO,
            ),
        ),
    )


def phenotype_exact_causal():
    return _case(
        case_id="phenotype-exact-001",
        family="phenotype",
        description=(
            "Causal candidate has the exact patient "
            "HPO association."
        ),
        expected_behavior=(
            "Exact phenotype mode should provide "
            "direct phenotype support."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="PDEC2",
                position=71011,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_DISTANT_HPO,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="PCAU2",
                position=71012,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_PATIENT_HPO,
                exact_phenotype_match=True,
            ),
        ),
    )


def phenotype_variant_pressure():
    return _case(
        case_id="phenotype-variant-pressure-001",
        family="phenotype",
        description=(
            "A strong variant-level decoy competes "
            "with an ontology-near causal candidate."
        ),
        expected_behavior=(
            "Semantic phenotype evidence should reduce "
            "the causal rank disadvantage."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="strong-decoy",
                gene="PDEC3",
                position=71021,
                consequence="frameshift_variant",
                impact="HIGH",
                significance="Pathogenic",
                af=1e-6,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_DISTANT_HPO,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="PCAU3",
                position=71022,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_NEAR_HPO,
            ),
        ),
    )


def phenotype_weak_decoys():
    return _case(
        case_id="phenotype-weak-decoys-001",
        family="phenotype",
        description=(
            "Ontology-near causal candidate competes "
            "with multiple weak phenotype decoys."
        ),
        expected_behavior=(
            "Semantic phenotype mode should prioritize "
            "the ontology-near candidate."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy-1",
                gene="PDEC4A",
                position=71031,
                consequence="intron_variant",
                impact="MODIFIER",
                af=0.005,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_DISTANT_HPO,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="PCAU4",
                position=71032,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_NEAR_HPO,
            ),
            make_synthetic_candidate(
                candidate_id="decoy-2",
                gene="PDEC4B",
                position=71033,
                consequence="intron_variant",
                impact="MODIFIER",
                af=0.01,
                patient_hpo=DEFAULT_PATIENT_HPO,
                associated_hpo=DEFAULT_DISTANT_HPO,
            ),
        ),
    )


def population_rare_causal():
    return _case(
        case_id="population-rare-001",
        family="population",
        description=(
            "Causal and decoy candidates differ "
            "primarily by population frequency."
        ),
        expected_behavior=(
            "Rare causal candidate should receive "
            "stronger population support."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="POP1D",
                position=72001,
                af=0.02,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="POP1C",
                position=72002,
                af=1e-6,
            ),
        ),
    )


def population_threshold():
    return _case(
        case_id="population-threshold-001",
        family="population",
        description=(
            "Candidates lie on different configured "
            "population-frequency tiers."
        ),
        expected_behavior=(
            "Lower-frequency causal candidate should "
            "receive greater population contribution."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="POP2D",
                position=72011,
                af=0.005,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="POP2C",
                position=72012,
                af=0.00005,
            ),
        ),
    )


def population_not_found():
    return _case(
        case_id="population-not-found-001",
        family="population",
        description=(
            "A gnomAD NOT_FOUND candidate competes "
            "with a numerically rare causal candidate."
        ),
        expected_behavior=(
            "NOT_FOUND must not be converted into AF=0 "
            "or stronger rarity evidence."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="not-found-decoy",
                gene="POP3D",
                position=72021,
                af=None,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="POP3C",
                position=72022,
                af=1e-6,
            ),
        ),
    )


def population_common_decoys():
    return _case(
        case_id="population-common-decoys-001",
        family="population",
        description=(
            "Rare causal candidate competes with "
            "multiple common decoys."
        ),
        expected_behavior=(
            "Population evidence should consistently "
            "favor the rare causal candidate."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy-1",
                gene="POP4A",
                position=72031,
                af=0.02,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="POP4C",
                position=72032,
                af=1e-6,
            ),
            make_synthetic_candidate(
                candidate_id="decoy-2",
                gene="POP4B",
                position=72033,
                af=0.05,
            ),
        ),
    )


def functional_frameshift_causal():
    return _case(
        case_id="functional-frameshift-001",
        family="functional",
        description=(
            "Causal candidate has a high-priority "
            "frameshift consequence."
        ),
        expected_behavior=(
            "Functional contribution should favor "
            "the frameshift candidate."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="FUN1D",
                position=73001,
                consequence="intron_variant",
                impact="MODIFIER",
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="FUN1C",
                position=73002,
                consequence="frameshift_variant",
                impact="HIGH",
            ),
        ),
    )


def functional_missense_causal():
    return _case(
        case_id="functional-missense-001",
        family="functional",
        description=(
            "Moderate missense causal candidate "
            "competes with a low-priority consequence."
        ),
        expected_behavior=(
            "Functional evidence should favor the "
            "missense candidate."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="FUN2D",
                position=73011,
                consequence="intron_variant",
                impact="MODIFIER",
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="FUN2C",
                position=73012,
                consequence="missense_variant",
                impact="MODERATE",
            ),
        ),
    )


def functional_balanced():
    return _case(
        case_id="functional-balanced-001",
        family="functional",
        description=(
            "Candidates have identical functional "
            "evidence."
        ),
        expected_behavior=(
            "Functional evidence alone should not "
            "create an artificial distinction."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="FUN3D",
                position=73021,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="FUN3C",
                position=73022,
            ),
        ),
    )


def functional_variant_pressure():
    return _case(
        case_id="functional-pressure-001",
        family="functional",
        description=(
            "Functional support favors the causal "
            "candidate while population evidence "
            "favors a decoy."
        ),
        expected_behavior=(
            "Case should expose competing ranking "
            "signals rather than encode an expected "
            "clinical conclusion."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="FUN4D",
                position=73031,
                consequence="intron_variant",
                impact="MODIFIER",
                af=1e-6,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="FUN4C",
                position=73032,
                consequence="frameshift_variant",
                impact="HIGH",
                af=0.005,
            ),
        ),
    )


def build_initial_registry():
    registry = ControlledSuiteRegistry()

    builders = (
        phenotype_semantic_near,
        phenotype_exact_causal,
        phenotype_variant_pressure,
        phenotype_weak_decoys,
        population_rare_causal,
        population_threshold,
        population_not_found,
        population_common_decoys,
        functional_frameshift_causal,
        functional_missense_causal,
        functional_balanced,
        functional_variant_pressure,
    )

    for builder in builders:
        case = builder()

        registry.register(
            case.case_id,
            builder,
        )

    return registry
