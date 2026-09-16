"""Reusable synthetic candidate factory for controlled ranking benchmarks.

Synthetic evidence is constructed from the same production evidence
contracts used by the interpretation and ranking layers.

No benchmark truth is stored in the evidence contract or candidate.
"""

from __future__ import annotations

from dataclasses import replace

from genomics_platform.evidence.variant_evidence_contract import (
    ClinicalConditionEvidence,
    ClinicalEvidence,
    ClinicalRecordEvidence,
    FunctionalEvidence,
    PopulationDatasetEvidence,
    PopulationEvidence,
    SourceAvailability,
    TranscriptEvidence,
    VariantEvidenceContract,
    VariantIdentity,
)
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


DEFAULT_PATIENT_HPO = "HP:0040154"
DEFAULT_NEAR_HPO = "HP:0001061"
DEFAULT_DISTANT_HPO = "HP:0000987"


def make_synthetic_contract(
    *,
    gene: str,
    position: int,
    consequence: str = "missense_variant",
    impact: str = "MODERATE",
    significance: str = "Uncertain_significance",
    af: float | None = 0.00001,
    chrom: str = "1",
    ref: str = "AG",
    alt: str = "A",
    assembly: str = "GRCh38",
    population_status: str | None = None,
) -> VariantEvidenceContract:
    """Construct one synthetic production evidence contract."""

    gene = gene.upper()

    transcript = TranscriptEvidence(
        gene_symbol=gene,
        gene_id=f"SYNTHETIC:{gene}",
        transcript_id=f"SYNTHETIC:{gene}:TRANSCRIPT",
        protein_id=f"SYNTHETIC:{gene}:PROTEIN",
        feature_type="Transcript",
        biotype="protein_coding",
        consequences=(consequence,),
        impact=impact,
        canonical=True,
        variant_class="synthetic",
    )

    functional = FunctionalEvidence(
        genes=(gene,),
        gene_ids=(f"SYNTHETIC:{gene}",),
        consequences=(consequence,),
        impacts=(impact,),
        transcripts=(transcript,),
    )

    clinical = ClinicalEvidence(
        records=(
            ClinicalRecordEvidence(
                variation_id=(
                    f"SYNTHETIC-{gene}-{position}"
                ),
                allele_id=(
                    f"SYNTHETIC-ALLELE-{position}"
                ),
                vcv_accession=None,
                clinical_significance=(
                    significance,
                ),
                conflicting_significance=(),
                review_status=(
                    "synthetic_controlled_fixture"
                ),
                review_stars=0,
                conditions=(
                    ClinicalConditionEvidence(
                        name=(
                            "Synthetic controlled "
                            "condition"
                        ),
                        identifiers=(
                            "SYNTHETIC:DISEASE",
                        ),
                    ),
                ),
                hgvs=(),
                scv_accessions=(),
            ),
        ),
    )

    if af is None:
        exome = None
        gnomad_status = (
            population_status
            if population_status is not None
            else "NOT_FOUND"
        )
    else:
        exome = PopulationDatasetEvidence(
            ac=1,
            an=100000,
            af=af,
            homozygote_count=0,
            populations=(),
        )
        gnomad_status = (
            population_status
            if population_status is not None
            else "FOUND"
        )

    population = PopulationEvidence(
        exome=exome,
        genome=None,
    )

    return VariantEvidenceContract(
        identity=VariantIdentity(
            assembly=assembly,
            chrom=chrom,
            pos=position,
            ref=ref,
            alt=alt,
        ),
        functional=functional,
        clinical=clinical,
        population=population,
        availability=(
            SourceAvailability(
                source="VEP",
                status="FOUND",
            ),
            SourceAvailability(
                source="ClinVar",
                status="FOUND",
            ),
            SourceAvailability(
                source="gnomAD",
                status=gnomad_status,
            ),
        ),
    )


def make_phenotype_dimension(
    *,
    patient_hpo: str,
    associated_hpo: str,
    exact_match: bool = False,
    gene: str | None = None,
    absent_hpo_terms: tuple[str, ...] = (),
) -> EvidenceDimension:
    """Construct controlled phenotype observations."""

    matched_present = (
        [patient_hpo]
        if exact_match
        else []
    )

    unmatched_present = (
        []
        if exact_match
        else [patient_hpo]
    )

    return EvidenceDimension(
        name="phenotype",
        status=(
            "AVAILABLE"
            if exact_match
            else "PARTIAL"
        ),
        observations=(
            InterpretationObservation(
                code=(
                    "PHENOTYPE_EXACT_MATCH_SUMMARY"
                ),
                label=(
                    "Synthetic exact phenotype "
                    "summary"
                ),
                value={
                    "matched_present": (
                        matched_present
                    ),
                    "unmatched_present": (
                        unmatched_present
                    ),
                    "matched_absent": [],
                    "unmatched_absent": list(
                        absent_hpo_terms
                    ),
                },
                source="controlled_fixture",
            ),
            InterpretationObservation(
                code=(
                    "PHENOTYPE_SOURCE_ASSOCIATION"
                ),
                label=(
                    "Synthetic HPO source "
                    "association"
                ),
                value={
                    "source": "HPO",
                    "gene_symbol": gene,
                    "disease_id": (
                        "SYNTHETIC:DISEASE"
                    ),
                    "phenotype_ids": [
                        associated_hpo
                    ],
                    "present_exact_matches": (
                        matched_present
                    ),
                    "absent_exact_matches": [
                        term
                        for term in absent_hpo_terms
                        if term == associated_hpo
                    ],
                },
                source="HPO",
            ),
        ),
        limitations=(
            "Synthetic controlled benchmark "
            "phenotype evidence.",
        ),
    )


def make_gene_context(
    *,
    gene: str,
    gencc_classification: str | None = "Definitive",
    panelapp_confidence: str | None = "Green",
    inheritance: tuple[str, ...] = (),
) -> dict:
    """Construct controlled gene-level evidence context."""

    records = []

    if gencc_classification is not None:
        evidence = {
            "source": "GenCC",
            "gene_symbol": gene,
            "disease_id": "MONDO:0000001",
            "disease_label": "Synthetic disease",
            "classification": gencc_classification,
        }

        if inheritance:
            evidence["inheritance"] = list(inheritance)

        records.append(
            {
                "evidence": evidence,
                "disease_normalization": {
                    "status": "FOUND",
                    "mondo_id": "MONDO:0000001",
                    "mondo_label": "Synthetic disease",
                },
            }
        )

    if panelapp_confidence is not None:
        evidence = {
            "source": "PanelApp",
            "gene_symbol": gene,
            "disease_id": "MONDO:0000001",
            "disease_label": "Synthetic disease",
            "confidence_level": panelapp_confidence,
        }

        if inheritance:
            evidence["inheritance"] = list(inheritance)

        records.append(
            {
                "evidence": evidence,
                "disease_normalization": {
                    "status": "FOUND",
                    "mondo_id": "MONDO:0000001",
                    "mondo_label": "Synthetic disease",
                },
            }
        )

    return {
        "query": {
            "gene_symbol": gene,
        },
        "records": records,
    }


def make_synthetic_candidate(
    *,
    candidate_id: str,
    gene: str,
    position: int,
    consequence: str = "missense_variant",
    impact: str = "MODERATE",
    significance: str = "Uncertain_significance",
    af: float | None = 0.00001,
    population_status: str | None = None,
    patient_hpo: str | None = None,
    associated_hpo: str | None = None,
    exact_phenotype_match: bool = False,
    gene_context: dict | None = None,
    genotype_state: str | None = None,
    proband_sex: str | None = None,
    de_novo_status: str | None = None,
    absent_hpo_terms: tuple[str, ...] = (),
    case_gene: str | None = None,
    allow_gene_mismatch: bool = False,
):
    """Build a synthetic CandidateCase without benchmark truth."""

    contract = make_synthetic_contract(
        gene=gene,
        position=position,
        consequence=consequence,
        impact=impact,
        significance=significance,
        af=af,
        population_status=population_status,
    )

    profile = interpret_variant(
        contract,
        gene_evidence_context=gene_context,
        genotype_state=genotype_state,
        proband_sex=proband_sex,
        de_novo_status=de_novo_status,
        absent_hpo_terms=absent_hpo_terms,
    )

    present_hpo_terms: tuple[str, ...] = ()

    if (
        patient_hpo is not None
        and associated_hpo is not None
    ):
        profile = replace(
            profile,
            phenotype=make_phenotype_dimension(
                patient_hpo=patient_hpo,
                associated_hpo=associated_hpo,
                exact_match=exact_phenotype_match,
                gene=gene,
                absent_hpo_terms=absent_hpo_terms,
            ),
        )

        present_hpo_terms = (patient_hpo,)

    candidate = build_candidate_case(
        candidate_id=candidate_id,
        evidence_contract=contract,
        interpretation_profile=profile,
        gene_symbol=(
            case_gene
            if case_gene is not None
            else gene
        ),
        present_hpo_terms=present_hpo_terms,
    )

    if (
        not allow_gene_mismatch
        and candidate.gene_consistency != "MATCH"
    ):
        raise ValueError(
            "Synthetic candidate unexpectedly "
            "failed gene consistency."
        )

    return candidate
