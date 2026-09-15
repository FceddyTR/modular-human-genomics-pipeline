from __future__ import annotations

from typing import Optional

from genomics_platform.evidence.variant_evidence_contract import (
    SourceAvailability,
    VariantEvidenceContract,
)
from genomics_platform.interpretation.interpretation_profile import (
    EvidenceDimension,
    InterpretationObservation,
)


def _get_availability(
    contract: VariantEvidenceContract,
    source: str,
) -> Optional[SourceAvailability]:
    for availability in contract.availability:
        if availability.source.lower() == source.lower():
            return availability

    return None


def interpret_functional(
    contract: VariantEvidenceContract,
) -> EvidenceDimension:
    """
    Build a descriptive functional annotation dimension from VEP.

    This interpreter intentionally does NOT:
      - classify pathogenicity,
      - treat VEP IMPACT as clinical evidence strength,
      - assume loss-of-function disease mechanism,
      - assign ACMG/AMP criteria,
      - assign ranking scores or candidate tiers,
      - infer damaging effect from consequence alone.
    """

    availability = _get_availability(
        contract,
        "VEP",
    )

    if availability is None:
        return EvidenceDimension(
            name="functional",
            status="UNAVAILABLE",
            limitations=(
                "No VEP availability record is present in "
                "the evidence contract.",
            ),
        )

    if availability.status == "ERROR":
        limitations = [
            "VEP annotation failed or was unavailable.",
            (
                "Functional annotation cannot be interpreted "
                "from this lookup."
            ),
        ]

        if availability.error:
            limitations.append(
                f"Source error: {availability.error}"
            )

        return EvidenceDimension(
            name="functional",
            status="ERROR",
            limitations=tuple(limitations),
        )

    if availability.status == "NOT_FOUND":
        return EvidenceDimension(
            name="functional",
            status="UNAVAILABLE",
            observations=(
                InterpretationObservation(
                    code="VEP_NOT_FOUND",
                    label="No matching VEP annotation returned",
                    value=None,
                    source="VEP",
                    details={
                        "meaning": (
                            "No normalized functional annotation "
                            "was returned for this lookup."
                        ),
                        "does_not_mean": [
                            "benign",
                            "functionally_neutral",
                            "clinically_irrelevant",
                        ],
                    },
                ),
            ),
            limitations=(
                (
                    "Missing VEP annotation does not establish "
                    "functional neutrality or benignity."
                ),
            ),
        )

    if availability.status != "FOUND":
        return EvidenceDimension(
            name="functional",
            status="ERROR",
            limitations=(
                (
                    "Unexpected VEP availability state: "
                    f"{availability.status}"
                ),
            ),
        )

    functional = contract.functional

    has_summary = any(
        (
            functional.genes,
            functional.gene_ids,
            functional.consequences,
            functional.impacts,
        )
    )

    has_transcripts = bool(
        functional.transcripts
    )

    if not has_summary and not has_transcripts:
        return EvidenceDimension(
            name="functional",
            status="PARTIAL",
            limitations=(
                (
                    "VEP lookup status is FOUND, but no normalized "
                    "functional annotation payload is available."
                ),
                (
                    "No functional or clinical effect is inferred "
                    "from the missing payload."
                ),
            ),
        )

    observations = []

    if functional.genes or functional.gene_ids:
        observations.append(
            InterpretationObservation(
                code="VEP_GENE_CONTEXT",
                label="VEP gene context",
                value={
                    "genes": list(functional.genes),
                    "gene_ids": list(functional.gene_ids),
                },
                source="VEP",
                details={
                    "meaning": (
                        "Genes associated with normalized VEP "
                        "transcript annotations."
                    ),
                },
            )
        )

    if functional.consequences:
        observations.append(
            InterpretationObservation(
                code="VEP_CONSEQUENCES",
                label="VEP sequence consequences",
                value=list(functional.consequences),
                source="VEP",
                details={
                    "meaning": (
                        "Sequence consequence annotations reported "
                        "by VEP."
                    ),
                    "does_not_mean": (
                        "Pathogenicity is not established by "
                        "consequence category alone."
                    ),
                },
            )
        )

    if functional.impacts:
        observations.append(
            InterpretationObservation(
                code="VEP_IMPACTS",
                label="VEP impact categories",
                value=list(functional.impacts),
                source="VEP",
                details={
                    "meaning": (
                        "VEP IMPACT annotation categories."
                    ),
                    "does_not_mean": (
                        "VEP IMPACT is not this platform's clinical "
                        "evidence strength or pathogenicity class."
                    ),
                },
            )
        )

    for transcript in functional.transcripts:
        observations.append(
            InterpretationObservation(
                code="VEP_TRANSCRIPT",
                label="VEP transcript annotation",
                value={
                    "transcript_id": transcript.transcript_id,
                    "gene_symbol": transcript.gene_symbol,
                    "gene_id": transcript.gene_id,
                    "consequences": list(
                        transcript.consequences
                    ),
                    "impact": transcript.impact,
                    "feature_type": transcript.feature_type,
                    "biotype": transcript.biotype,
                    "protein_id": transcript.protein_id,
                    "canonical": transcript.canonical,
                    "variant_class": transcript.variant_class,
                },
                source="VEP",
                details={
                    "meaning": (
                        "Transcript-specific functional annotation."
                    ),
                    "clinical_interpretation": False,
                },
            )
        )

    limitations = (
        (
            "VEP consequence annotations describe predicted "
            "sequence consequences and do not independently "
            "establish pathogenicity."
        ),
        (
            "VEP IMPACT categories are annotation metadata and "
            "are not interpreted as clinical evidence strength."
        ),
        (
            "Transcript relevance, disease mechanism, gene-disease "
            "validity and phenotype compatibility are not evaluated "
            "at this layer."
        ),
    )

    return EvidenceDimension(
        name="functional",
        status="AVAILABLE",
        observations=tuple(observations),
        limitations=limitations,
    )
