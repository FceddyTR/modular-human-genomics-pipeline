"""Controlled inheritance ranking stress cases."""

from __future__ import annotations

from genomics_platform.ranking.benchmark.controlled_suite import (
    ControlledCase,
    ControlledCaseSpec,
)
from genomics_platform.ranking.benchmark.synthetic_factory import (
    make_gene_context,
    make_synthetic_candidate,
)


def _case(
    *,
    case_id: str,
    description: str,
    expected_behavior: str,
    candidates,
) -> ControlledCase:
    return ControlledCase(
        spec=ControlledCaseSpec(
            case_id=case_id,
            family="inheritance",
            description=description,
            causal_candidate_ids=("causal",),
            expected_behavior=expected_behavior,
            tags=(
                "synthetic",
                "controlled",
                "stress",
                "inheritance",
            ),
        ),
        candidates=tuple(candidates),
    )


def ad_heterozygous_compatible():
    causal_context = make_gene_context(
        gene="INHC1",
        inheritance=("Autosomal dominant",),
    )

    return _case(
        case_id="inheritance-ad-compatible-001",
        description=(
            "Heterozygous causal candidate has "
            "autosomal-dominant source evidence."
        ),
        expected_behavior=(
            "Compatible inheritance should contribute "
            "relative ranking support without implying "
            "pathogenicity."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="INHD1",
                position=75001,
                genotype_state="heterozygous",
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="INHC1",
                position=75002,
                gene_context=causal_context,
                genotype_state="heterozygous",
            ),
        ),
    )


def ad_variant_pressure():
    causal_context = make_gene_context(
        gene="INHC2",
        inheritance=("Autosomal dominant",),
    )

    return _case(
        case_id="inheritance-ad-pressure-001",
        description=(
            "Inheritance-compatible causal candidate "
            "competes with stronger variant evidence."
        ),
        expected_behavior=(
            "The case measures whether inheritance "
            "support can reduce a variant-level "
            "ranking disadvantage."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="INHD2",
                position=75011,
                consequence="frameshift_variant",
                impact="HIGH",
                significance="Pathogenic",
                af=1e-6,
                genotype_state="heterozygous",
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="INHC2",
                position=75012,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
                gene_context=causal_context,
                genotype_state="heterozygous",
            ),
        ),
    )


def competing_inheritance_context():
    decoy_context = make_gene_context(
        gene="INHD3",
        inheritance=("Autosomal recessive",),
    )

    causal_context = make_gene_context(
        gene="INHC3",
        inheritance=("Autosomal dominant",),
    )

    return _case(
        case_id="inheritance-competing-context-001",
        description=(
            "Candidates carry different inheritance "
            "models under the same heterozygous case "
            "context."
        ),
        expected_behavior=(
            "Inheritance compatibility should remain "
            "a distinct ranking dimension."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="INHD3",
                position=75021,
                gene_context=decoy_context,
                genotype_state="heterozygous",
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="INHC3",
                position=75022,
                gene_context=causal_context,
                genotype_state="heterozygous",
            ),
        ),
    )


def inheritance_with_three_candidates():
    causal_context = make_gene_context(
        gene="INHC4",
        inheritance=("Autosomal dominant",),
    )

    return _case(
        case_id="inheritance-three-candidate-001",
        description=(
            "Inheritance-supported causal candidate "
            "competes with two variant-level decoys."
        ),
        expected_behavior=(
            "The case tests inheritance contribution "
            "with a larger candidate set."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy-1",
                gene="INHD4A",
                position=75031,
                consequence="intron_variant",
                impact="MODIFIER",
                genotype_state="heterozygous",
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="INHC4",
                position=75032,
                gene_context=causal_context,
                genotype_state="heterozygous",
            ),
            make_synthetic_candidate(
                candidate_id="decoy-2",
                gene="INHD4B",
                position=75033,
                consequence="missense_variant",
                impact="MODERATE",
                genotype_state="heterozygous",
            ),
        ),
    )


BUILDERS = (
    ad_heterozygous_compatible,
    ad_variant_pressure,
    competing_inheritance_context,
    inheritance_with_three_candidates,
)
