"""Controlled gene-disease ranking stress cases."""

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
            family="gene_disease",
            description=description,
            causal_candidate_ids=("causal",),
            expected_behavior=expected_behavior,
            tags=(
                "synthetic",
                "controlled",
                "stress",
                "gene_disease",
            ),
        ),
        candidates=tuple(candidates),
    )


def strong_causal():
    """Strong gene validity favors otherwise balanced causal."""

    causal_context = make_gene_context(
        gene="GDC1",
        gencc_classification="Definitive",
        panelapp_confidence="Green",
    )

    return _case(
        case_id="gene-disease-strong-causal-001",
        description=(
            "Balanced variants differ primarily in "
            "gene-disease evidence."
        ),
        expected_behavior=(
            "Strong gene-disease evidence should add "
            "relative prioritization support without "
            "becoming variant pathogenicity evidence."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="GDD1",
                position=74001,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="GDC1",
                position=74002,
                gene_context=causal_context,
            ),
        ),
    )


def moderate_vs_strong():
    """Compare moderate and stronger source evidence."""

    decoy_context = make_gene_context(
        gene="GDD2",
        gencc_classification="Moderate",
        panelapp_confidence=None,
    )

    causal_context = make_gene_context(
        gene="GDC2",
        gencc_classification="Definitive",
        panelapp_confidence="Green",
    )

    return _case(
        case_id="gene-disease-strength-001",
        description=(
            "Both candidates have gene-disease "
            "evidence of different source strength."
        ),
        expected_behavior=(
            "The ranking component should preserve "
            "relative gene-level evidence strength."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="GDD2",
                position=74011,
                gene_context=decoy_context,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="GDC2",
                position=74012,
                gene_context=causal_context,
            ),
        ),
    )


def gene_evidence_vs_variant_pressure():
    """Strong variant decoy competes with gene-supported causal."""

    causal_context = make_gene_context(
        gene="GDC3",
        gencc_classification="Definitive",
        panelapp_confidence="Green",
    )

    return _case(
        case_id="gene-disease-variant-pressure-001",
        description=(
            "Strong variant-level decoy competes with "
            "a causal candidate carrying stronger "
            "gene-disease context."
        ),
        expected_behavior=(
            "The case measures conflict between "
            "variant-level and gene-level signals."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="GDD3",
                position=74021,
                consequence="frameshift_variant",
                impact="HIGH",
                significance="Pathogenic",
                af=1e-6,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="GDC3",
                position=74022,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
                gene_context=causal_context,
            ),
        ),
    )


def source_disagreement():
    """Mixed gene sources test conservative aggregation."""

    decoy_context = make_gene_context(
        gene="GDD4",
        gencc_classification="Moderate",
        panelapp_confidence="Green",
    )

    causal_context = make_gene_context(
        gene="GDC4",
        gencc_classification="Definitive",
        panelapp_confidence=None,
    )

    return _case(
        case_id="gene-disease-source-mix-001",
        description=(
            "Candidates receive differently composed "
            "GenCC and PanelApp evidence."
        ),
        expected_behavior=(
            "Source evidence should remain descriptive "
            "and produce deterministic ranking input."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="GDD4",
                position=74031,
                gene_context=decoy_context,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="GDC4",
                position=74032,
                gene_context=causal_context,
            ),
        ),
    )


BUILDERS = (
    strong_causal,
    moderate_vs_strong,
    gene_evidence_vs_variant_pressure,
    source_disagreement,
)
