from __future__ import annotations

from typing import Optional

from genomics_platform.evidence.variant_evidence_contract import (
    PopulationDatasetEvidence,
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


def _dataset_observations(
    dataset_name: str,
    dataset: PopulationDatasetEvidence,
) -> list[InterpretationObservation]:
    prefix = dataset_name.upper()

    observations = [
        InterpretationObservation(
            code=f"GNOMAD_{prefix}_SUMMARY",
            label=f"gnomAD {dataset_name} population observation",
            value={
                "ac": dataset.ac,
                "an": dataset.an,
                "af": dataset.af,
                "homozygote_count": dataset.homozygote_count,
            },
            source="gnomAD",
            details={
                "dataset": dataset_name,
                "interpretation": "descriptive_population_evidence",
            },
        )
    ]

    for population in dataset.populations:
        observations.append(
            InterpretationObservation(
                code=f"GNOMAD_{prefix}_POPULATION",
                label=(
                    f"gnomAD {dataset_name} "
                    f"population observation"
                ),
                value={
                    "population_id": population.population_id,
                    "ac": population.ac,
                    "an": population.an,
                    "af": population.af,
                    "homozygote_count": (
                        population.homozygote_count
                    ),
                },
                source="gnomAD",
                details={
                    "dataset": dataset_name,
                    "population_id": population.population_id,
                },
            )
        )

    return observations


def interpret_population(
    contract: VariantEvidenceContract,
) -> EvidenceDimension:
    """
    Build a descriptive population evidence dimension.

    This interpreter intentionally does not:
      - classify a variant as rare or common,
      - infer AF=0 from NOT_FOUND,
      - infer pathogenicity from low frequency,
      - merge exome and genome observations,
      - infer ancestry grpmax from arbitrary population IDs.

    Disease-specific frequency thresholds belong to a later
    interpretation policy layer.
    """

    availability = _get_availability(
        contract,
        "gnomAD",
    )

    if availability is None:
        return EvidenceDimension(
            name="population",
            status="UNAVAILABLE",
            limitations=(
                "No gnomAD availability record is present in "
                "the evidence contract.",
            ),
        )

    if availability.status == "ERROR":
        limitations = [
            "gnomAD evidence lookup failed or was unavailable.",
            (
                "Population evidence cannot be interpreted from "
                "this lookup."
            ),
        ]

        if availability.error:
            limitations.append(
                f"Source error: {availability.error}"
            )

        return EvidenceDimension(
            name="population",
            status="ERROR",
            limitations=tuple(limitations),
        )

    if availability.status == "NOT_FOUND":
        return EvidenceDimension(
            name="population",
            status="UNAVAILABLE",
            observations=(
                InterpretationObservation(
                    code="GNOMAD_NOT_FOUND",
                    label="No matching gnomAD record returned",
                    value=None,
                    source="gnomAD",
                    details={
                        "meaning": (
                            "Lookup completed without a matching "
                            "gnomAD record."
                        ),
                        "does_not_mean": [
                            "AF=0",
                            "rare",
                            "pathogenic",
                            "absent_from_population",
                        ],
                    },
                ),
            ),
            limitations=(
                (
                    "gnomAD NOT_FOUND does not establish allele "
                    "frequency zero or rarity."
                ),
                (
                    "Coverage and callability are not established "
                    "by this result."
                ),
            ),
        )

    if availability.status != "FOUND":
        return EvidenceDimension(
            name="population",
            status="ERROR",
            limitations=(
                (
                    "Unexpected gnomAD availability state: "
                    f"{availability.status}"
                ),
            ),
        )

    observations = []

    if contract.population.exome is not None:
        observations.extend(
            _dataset_observations(
                "exome",
                contract.population.exome,
            )
        )

    if contract.population.genome is not None:
        observations.extend(
            _dataset_observations(
                "genome",
                contract.population.genome,
            )
        )

    if not observations:
        return EvidenceDimension(
            name="population",
            status="PARTIAL",
            limitations=(
                (
                    "gnomAD lookup status is FOUND, but no exome "
                    "or genome frequency payload is available."
                ),
                (
                    "No rarity or pathogenicity inference is made "
                    "from the missing payload."
                ),
            ),
        )

    missing_datasets = []

    if contract.population.exome is None:
        missing_datasets.append("exome")

    if contract.population.genome is None:
        missing_datasets.append("genome")

    limitations = []

    if missing_datasets:
        limitations.append(
            "No frequency payload for: "
            + ", ".join(missing_datasets)
            + "."
        )

    limitations.extend(
        [
            (
                "Exome and genome observations are preserved "
                "separately."
            ),
            (
                "No disease-specific rarity threshold is applied "
                "at this layer."
            ),
            (
                "Population frequency alone does not establish "
                "pathogenicity or benignity."
            ),
        ]
    )

    return EvidenceDimension(
        name="population",
        status=(
            "PARTIAL"
            if missing_datasets
            else "AVAILABLE"
        ),
        observations=tuple(observations),
        limitations=tuple(limitations),
    )
