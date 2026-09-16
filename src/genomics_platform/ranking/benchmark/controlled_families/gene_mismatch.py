"""Controlled gene-consistency mismatch stress cases."""

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
            family="gene_mismatch",
            description=description,
            causal_candidate_ids=("causal",),
            expected_behavior=expected_behavior,
            tags=(
                "synthetic",
                "controlled",
                "stress",
                "gene_mismatch",
            ),
        ),
        candidates=tuple(candidates),
    )


def gene_disease_blocked():
    mismatch_context = make_gene_context(
        gene="GMCON1",
        gencc_classification="Definitive",
        panelapp_confidence="Green",
    )

    return _case(
        case_id="gene-mismatch-disease-block-001",
        description=(
            "A candidate carries strong gene-disease "
            "context but its CandidateCase gene does "
            "not match the evidence-contract gene."
        ),
        expected_behavior=(
            "Gene-dependent ranking support must be "
            "blocked by gene consistency mismatch."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="causal",
                gene="GMCON1",
                case_gene="GMCASE1",
                position=76001,
                gene_context=mismatch_context,
                allow_gene_mismatch=True,
            ),
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="GMDEC1",
                position=76002,
            ),
        ),
    )


def inheritance_blocked():
    mismatch_context = make_gene_context(
        gene="GMCON2",
        inheritance=("Autosomal dominant",),
    )

    return _case(
        case_id="gene-mismatch-inheritance-block-001",
        description=(
            "Inheritance-compatible source evidence "
            "is attached to a gene-mismatched case."
        ),
        expected_behavior=(
            "Inheritance ranking support must be "
            "blocked when gene consistency fails."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="causal",
                gene="GMCON2",
                case_gene="GMCASE2",
                position=76011,
                gene_context=mismatch_context,
                genotype_state="heterozygous",
                allow_gene_mismatch=True,
            ),
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="GMDEC2",
                position=76012,
                genotype_state="heterozygous",
            ),
        ),
    )


def variant_evidence_survives():
    mismatch_context = make_gene_context(
        gene="GMCON3",
        gencc_classification="Definitive",
        panelapp_confidence="Green",
    )

    return _case(
        case_id="gene-mismatch-variant-survives-001",
        description=(
            "Gene mismatch blocks gene-dependent "
            "signals while variant-level evidence "
            "remains available."
        ),
        expected_behavior=(
            "Gene mismatch must not erase independent "
            "clinical, population, or functional "
            "variant evidence."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="causal",
                gene="GMCON3",
                case_gene="GMCASE3",
                position=76021,
                consequence="frameshift_variant",
                impact="HIGH",
                significance="Pathogenic",
                af=1e-6,
                gene_context=mismatch_context,
                allow_gene_mismatch=True,
            ),
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="GMDEC3",
                position=76022,
                consequence="intron_variant",
                impact="MODIFIER",
                significance="Uncertain_significance",
                af=0.005,
            ),
        ),
    )


def matched_control():
    causal_context = make_gene_context(
        gene="GMC4",
        gencc_classification="Definitive",
        panelapp_confidence="Green",
        inheritance=("Autosomal dominant",),
    )

    return _case(
        case_id="gene-mismatch-matched-control-001",
        description=(
            "Matched control demonstrates that the "
            "same gene-dependent evidence remains "
            "eligible when consistency is preserved."
        ),
        expected_behavior=(
            "Matched gene context should remain "
            "available to gene-disease and inheritance "
            "ranking components."
        ),
        candidates=(
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="GMD4",
                position=76031,
                genotype_state="heterozygous",
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="GMC4",
                position=76032,
                gene_context=causal_context,
                genotype_state="heterozygous",
            ),
        ),
    )


BUILDERS = (
    gene_disease_blocked,
    inheritance_blocked,
    variant_evidence_survives,
    matched_control,
)
