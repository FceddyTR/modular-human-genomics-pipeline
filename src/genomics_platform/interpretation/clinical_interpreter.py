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


def interpret_clinical(
    contract: VariantEvidenceContract,
) -> EvidenceDimension:
    """
    Build a descriptive clinical-database evidence dimension.

    This interpreter reports what ClinVar contains.

    It intentionally does NOT:
      - adopt ClinVar significance as this platform's classification,
      - implement ACMG/AMP criteria,
      - convert review stars into pathogenicity strength,
      - infer benignity from NOT_FOUND,
      - rank variants,
      - assign candidate tiers.
    """

    availability = _get_availability(
        contract,
        "ClinVar",
    )

    if availability is None:
        return EvidenceDimension(
            name="clinical",
            status="UNAVAILABLE",
            limitations=(
                "No ClinVar availability record is present in "
                "the evidence contract.",
            ),
        )

    if availability.status == "ERROR":
        limitations = [
            "ClinVar evidence lookup failed or was unavailable.",
            (
                "No clinical-database interpretation can be made "
                "from this lookup."
            ),
        ]

        if availability.error:
            limitations.append(
                f"Source error: {availability.error}"
            )

        return EvidenceDimension(
            name="clinical",
            status="ERROR",
            limitations=tuple(limitations),
        )

    if availability.status == "NOT_FOUND":
        return EvidenceDimension(
            name="clinical",
            status="UNAVAILABLE",
            observations=(
                InterpretationObservation(
                    code="CLINVAR_NOT_FOUND",
                    label="No matching ClinVar record returned",
                    value=None,
                    source="ClinVar",
                    details={
                        "meaning": (
                            "Lookup completed without a matching "
                            "ClinVar record."
                        ),
                        "does_not_mean": [
                            "benign",
                            "likely_benign",
                            "clinically_irrelevant",
                            "no_disease_association",
                        ],
                    },
                ),
            ),
            limitations=(
                (
                    "ClinVar NOT_FOUND does not establish benignity "
                    "or clinical irrelevance."
                ),
                (
                    "Absence from ClinVar may reflect incomplete "
                    "submission or curation coverage."
                ),
            ),
        )

    if availability.status != "FOUND":
        return EvidenceDimension(
            name="clinical",
            status="ERROR",
            limitations=(
                (
                    "Unexpected ClinVar availability state: "
                    f"{availability.status}"
                ),
            ),
        )

    if not contract.clinical.records:
        return EvidenceDimension(
            name="clinical",
            status="PARTIAL",
            limitations=(
                (
                    "ClinVar lookup status is FOUND, but no clinical "
                    "records are present in the normalized contract."
                ),
                (
                    "No clinical classification is inferred from "
                    "the missing payload."
                ),
            ),
        )

    observations = []

    submitted_significances = set()
    review_statuses = set()

    for record in contract.clinical.records:
        if record.clinical_significance:
            submitted_significances.update(
                record.clinical_significance
            )

            observations.append(
                InterpretationObservation(
                    code="CLINVAR_SUBMITTED_SIGNIFICANCE",
                    label="ClinVar submitted clinical significance",
                    value=list(record.clinical_significance),
                    source="ClinVar",
                    details={
                        "variation_id": record.variation_id,
                        "vcv_accession": record.vcv_accession,
                        "interpretation": (
                            "source_reported_significance"
                        ),
                        "platform_classification": False,
                    },
                )
            )

        if record.conflicting_significance:
            observations.append(
                InterpretationObservation(
                    code="CLINVAR_CONFLICTING_SIGNIFICANCE",
                    label="ClinVar conflicting significance",
                    value=list(
                        record.conflicting_significance
                    ),
                    source="ClinVar",
                    details={
                        "variation_id": record.variation_id,
                        "vcv_accession": record.vcv_accession,
                    },
                )
            )

        if record.review_status is not None:
            review_statuses.add(record.review_status)

            observations.append(
                InterpretationObservation(
                    code="CLINVAR_REVIEW_STATUS",
                    label="ClinVar review status",
                    value={
                        "review_status": record.review_status,
                        "review_stars": record.review_stars,
                    },
                    source="ClinVar",
                    details={
                        "variation_id": record.variation_id,
                        "vcv_accession": record.vcv_accession,
                        "meaning": (
                            "ClinVar review level metadata; not "
                            "pathogenicity strength."
                        ),
                    },
                )
            )

        if record.conditions:
            observations.append(
                InterpretationObservation(
                    code="CLINVAR_CONDITIONS",
                    label="ClinVar associated conditions",
                    value=[
                        condition.to_dict()
                        for condition in record.conditions
                    ],
                    source="ClinVar",
                    details={
                        "variation_id": record.variation_id,
                        "vcv_accession": record.vcv_accession,
                        "meaning": (
                            "Conditions attached to the ClinVar "
                            "record; this does not independently "
                            "establish causality for the patient."
                        ),
                    },
                )
            )

        if record.hgvs:
            observations.append(
                InterpretationObservation(
                    code="CLINVAR_HGVS",
                    label="ClinVar HGVS expressions",
                    value=list(record.hgvs),
                    source="ClinVar",
                    details={
                        "variation_id": record.variation_id,
                    },
                )
            )

    if len(submitted_significances) > 1:
        observations.append(
            InterpretationObservation(
                code="CLINVAR_MULTIPLE_SUBMITTED_SIGNIFICANCES",
                label="Multiple ClinVar significance values observed",
                value=sorted(submitted_significances),
                source="ClinVar",
                details={
                    "meaning": (
                        "Multiple source-reported significance "
                        "values are present. No platform consensus "
                        "classification is inferred."
                    ),
                },
            )
        )

    limitations = (
        (
            "ClinVar significance values are source-reported "
            "evidence and are not adopted as this platform's "
            "independent classification."
        ),
        (
            "ClinVar review stars describe review status and are "
            "not converted into pathogenicity strength."
        ),
        (
            "Condition associations do not independently establish "
            "causality for the analyzed individual."
        ),
    )

    return EvidenceDimension(
        name="clinical",
        status="AVAILABLE",
        observations=tuple(observations),
        limitations=limitations,
    )
