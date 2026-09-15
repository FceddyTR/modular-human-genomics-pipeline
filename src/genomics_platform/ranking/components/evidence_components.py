"""Evidence components for candidate prioritization.

All outputs are heuristic prioritization signals only.

They are not:
- pathogenicity probabilities,
- ACMG/AMP criteria,
- variant classifications,
- diagnoses.
"""

from __future__ import annotations

from genomics_platform.interpretation.case_model import (
    CandidateCase,
)
from genomics_platform.ranking.ranking_contract import (
    ComponentContribution,
)


def _observations(
    dimension,
    code: str,
):
    if dimension is None:
        return ()

    return tuple(
        observation
        for observation
        in dimension.observations
        if observation.code == code
    )


def unavailable(
    name: str,
    maximum: float,
    reason: str,
) -> ComponentContribution:
    return ComponentContribution(
        name=name,
        contribution=0.0,
        max_contribution=maximum,
        status="UNAVAILABLE",
        limitations=(reason,),
    )


def blocked(
    name: str,
    maximum: float,
    reason: str,
) -> ComponentContribution:
    return ComponentContribution(
        name=name,
        contribution=0.0,
        max_contribution=maximum,
        status="BLOCKED",
        limitations=(reason,),
    )


def clinical_component(
    candidate: CandidateCase,
    maximum: float,
) -> ComponentContribution:
    dimension = (
        candidate.interpretation_profile
        .clinical
    )

    if dimension.status == "ERROR":
        return ComponentContribution(
            name="clinical",
            contribution=0.0,
            max_contribution=maximum,
            status="ERROR",
            limitations=(
                "Clinical evidence interpretation "
                "contains an error.",
            ),
        )

    observations = _observations(
        dimension,
        "CLINVAR_SUBMITTED_SIGNIFICANCE",
    )

    if not observations:
        return unavailable(
            "clinical",
            maximum,
            "No ClinVar submitted significance "
            "is available.",
        )

    values = set()

    for observation in observations:
        raw = observation.value

        if isinstance(raw, (list, tuple)):
            values.update(
                str(value)
                .strip()
                .lower()
                .replace(" ", "_")
                for value in raw
            )

        elif raw is not None:
            values.add(
                str(raw)
                .strip()
                .lower()
                .replace(" ", "_")
            )

    if values & {
        "pathogenic",
        "likely_pathogenic",
    }:
        fraction = 1.0
        reason = (
            "ClinVar contains submitted "
            "Pathogenic/Likely pathogenic "
            "significance."
        )

    elif values & {
        "uncertain_significance",
        "vus",
    }:
        fraction = 0.25
        reason = (
            "ClinVar contains submitted VUS "
            "significance."
        )

    elif values & {
        "benign",
        "likely_benign",
    }:
        fraction = 0.0
        reason = (
            "ClinVar contains submitted "
            "Benign/Likely benign significance."
        )

    else:
        fraction = 0.0
        reason = (
            "ClinVar significance is available "
            "but has no v0.2 prioritization rule."
        )

    contribution = maximum * fraction

    return ComponentContribution(
        name="clinical",
        contribution=round(
            contribution,
            6,
        ),
        max_contribution=maximum,
        status=(
            "SUPPORTING"
            if contribution > 0
            else "NEUTRAL"
        ),
        reasons=(reason,),
        limitations=(
            "ClinVar assertions are source values, "
            "not platform classifications.",
            "No ACMG/AMP classification is performed.",
        ),
    )


def population_component(
    candidate: CandidateCase,
    maximum: float,
) -> ComponentContribution:
    dimension = (
        candidate.interpretation_profile
        .population
    )

    if dimension.status == "ERROR":
        return ComponentContribution(
            name="population",
            contribution=0.0,
            max_contribution=maximum,
            status="ERROR",
            limitations=(
                "Population evidence interpretation "
                "contains an error.",
            ),
        )

    observations = (
        _observations(
            dimension,
            "GNOMAD_EXOME_SUMMARY",
        )
        + _observations(
            dimension,
            "GNOMAD_GENOME_SUMMARY",
        )
    )

    frequencies = []

    for observation in observations:
        value = observation.value or {}
        af = value.get("af")

        if isinstance(af, (int, float)):
            frequencies.append(
                float(af)
            )

    if not frequencies:
        return unavailable(
            "population",
            maximum,
            "No numeric gnomAD AF is available. "
            "NOT_FOUND is not interpreted as AF=0.",
        )

    max_af = max(frequencies)

    if max_af <= 0.0001:
        fraction = 1.0
        label = "very low observed AF"

    elif max_af <= 0.001:
        fraction = 0.75
        label = "low observed AF"

    elif max_af <= 0.01:
        fraction = 0.25
        label = "intermediate observed AF"

    else:
        fraction = 0.0
        label = "higher observed AF"

    contribution = maximum * fraction

    return ComponentContribution(
        name="population",
        contribution=round(
            contribution,
            6,
        ),
        max_contribution=maximum,
        status=(
            "SUPPORTING"
            if contribution > 0
            else "NEUTRAL"
        ),
        reasons=(
            f"Maximum observed AF={max_af:g}; "
            f"{label}.",
        ),
        limitations=(
            "Population rarity is not pathogenicity.",
            "Coverage, ancestry and callability are "
            "not fully modeled in v0.2.",
            "gnomAD NOT_FOUND is never converted "
            "to AF=0.",
        ),
    )


def functional_component(
    candidate: CandidateCase,
    maximum: float,
) -> ComponentContribution:
    dimension = (
        candidate.interpretation_profile
        .functional
    )

    if dimension.status == "ERROR":
        return ComponentContribution(
            name="functional",
            contribution=0.0,
            max_contribution=maximum,
            status="ERROR",
            limitations=(
                "Functional annotation contains "
                "an error.",
            ),
        )

    observations = _observations(
        dimension,
        "VEP_CONSEQUENCES",
    )

    consequences = set()

    for observation in observations:
        raw = observation.value

        if isinstance(raw, (list, tuple)):
            consequences.update(
                str(value).lower()
                for value in raw
            )

    if not consequences:
        return unavailable(
            "functional",
            maximum,
            "No VEP consequence is available.",
        )

    high = {
        "stop_gained",
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "start_lost",
        "stop_lost",
    }

    moderate = {
        "missense_variant",
        "inframe_deletion",
        "inframe_insertion",
        "protein_altering_variant",
    }

    if consequences & high:
        fraction = 1.0
        reason = (
            "Variant has a v0.2 high-priority "
            "VEP consequence."
        )

    elif consequences & moderate:
        fraction = 0.6
        reason = (
            "Variant has a v0.2 moderate-priority "
            "VEP consequence."
        )

    else:
        fraction = 0.1
        reason = (
            "Variant is annotated but has no "
            "v0.2 high/moderate-priority consequence."
        )

    contribution = maximum * fraction

    return ComponentContribution(
        name="functional",
        contribution=round(
            contribution,
            6,
        ),
        max_contribution=maximum,
        status=(
            "SUPPORTING"
            if contribution > 0
            else "NEUTRAL"
        ),
        reasons=(reason,),
        limitations=(
            "VEP consequence is not a pathogenicity "
            "classification.",
            "Mechanism-of-disease compatibility is "
            "not modeled in v0.2.",
        ),
    )


def gene_disease_component(
    candidate: CandidateCase,
    maximum: float,
) -> ComponentContribution:
    if candidate.gene_consistency != "MATCH":
        return blocked(
            "gene_disease",
            maximum,
            "Gene-dependent evidence blocked because "
            f"gene consistency is "
            f"{candidate.gene_consistency}.",
        )

    dimension = (
        candidate.interpretation_profile
        .gene_disease
    )

    if dimension is None:
        return unavailable(
            "gene_disease",
            maximum,
            "Gene-disease context was not evaluated.",
        )

    records = _observations(
        dimension,
        "GENE_DISEASE_SOURCE_RECORD",
    )

    if not records:
        return unavailable(
            "gene_disease",
            maximum,
            "No gene-disease source records "
            "are available.",
        )

    gencc = []
    panelapp = []
    orphanet_assessed = False

    for observation in records:
        value = observation.value or {}
        evidence = value.get(
            "evidence"
        ) or {}

        source = (
            value.get("source")
            or observation.source
        )

        if source == "GenCC":
            classification = (
                evidence.get(
                    "classification"
                )
            )

            if classification:
                gencc.append(
                    str(classification)
                    .strip()
                    .lower()
                )

        elif source == "PanelApp":
            confidence = (
                evidence.get(
                    "confidence"
                )
            )

            if confidence:
                panelapp.append(
                    str(confidence)
                    .strip()
                    .lower()
                )

        elif source == "Orphanet":
            confidence = (
                evidence.get(
                    "confidence"
                )
            )

            if (
                isinstance(
                    confidence,
                    str,
                )
                and confidence.lower()
                == "assessed"
            ):
                orphanet_assessed = True

    fraction = 0.0
    reasons = []

    if any(
        value in {
            "definitive",
            "strong",
        }
        for value in gencc
    ):
        fraction = max(
            fraction,
            1.0,
        )

        reasons.append(
            "GenCC contains Strong/Definitive "
            "gene-disease validity evidence."
        )

    elif "moderate" in gencc:
        fraction = max(
            fraction,
            0.7,
        )

        reasons.append(
            "GenCC contains Moderate "
            "gene-disease validity evidence."
        )

    if "green" in panelapp:
        fraction = max(
            fraction,
            0.7,
        )

        reasons.append(
            "PanelApp contains Green "
            "gene-level evidence."
        )

    elif "amber" in panelapp:
        fraction = max(
            fraction,
            0.4,
        )

        reasons.append(
            "PanelApp contains Amber "
            "gene-level evidence."
        )

    if orphanet_assessed:
        fraction = max(
            fraction,
            0.4,
        )

        reasons.append(
            "Orphanet contains an assessed "
            "gene-disease association."
        )

    contribution = maximum * fraction

    return ComponentContribution(
        name="gene_disease",
        contribution=round(
            contribution,
            6,
        ),
        max_contribution=maximum,
        status=(
            "SUPPORTING"
            if contribution > 0
            else "NEUTRAL"
        ),
        reasons=tuple(reasons),
        limitations=(
            "Gene-disease validity is not variant "
            "pathogenicity evidence.",
            "PanelApp confidence is gene/panel "
            "context, not ACMG classification.",
        ),
    )


def inheritance_component(
    candidate: CandidateCase,
    maximum: float,
) -> ComponentContribution:
    if candidate.gene_consistency != "MATCH":
        return blocked(
            "inheritance",
            maximum,
            "Gene-dependent evidence blocked because "
            f"gene consistency is "
            f"{candidate.gene_consistency}.",
        )

    dimension = (
        candidate.interpretation_profile
        .inheritance
    )

    if dimension is None:
        return unavailable(
            "inheritance",
            maximum,
            "Inheritance context was not evaluated.",
        )

    records = _observations(
        dimension,
        "INHERITANCE_SOURCE_EVIDENCE",
    )

    if not records:
        return unavailable(
            "inheritance",
            maximum,
            "No inheritance source evidence "
            "is available.",
        )

    states = [
        (record.value or {}).get(
            "genotype_compatibility"
        )
        for record in records
    ]

    compatible = sum(
        state == "COMPATIBLE"
        for state in states
    )

    partial = sum(
        state == "PARTIAL"
        for state in states
    )

    evaluable = (
        compatible
        + partial
    )

    if evaluable == 0:
        return unavailable(
            "inheritance",
            maximum,
            "Inheritance evidence exists but "
            "compatibility is not evaluable.",
        )

    fraction = (
        compatible
        + (0.5 * partial)
    ) / len(records)

    contribution = maximum * fraction

    return ComponentContribution(
        name="inheritance",
        contribution=round(
            contribution,
            6,
        ),
        max_contribution=maximum,
        status=(
            "SUPPORTING"
            if contribution > 0
            else "NEUTRAL"
        ),
        reasons=(
            f"{compatible} compatible and "
            f"{partial} partial inheritance "
            f"observation(s) among "
            f"{len(records)} record(s).",
        ),
        limitations=(
            "Inheritance compatibility does not "
            "establish pathogenicity.",
            "Phase, second allele, CNV/SV and "
            "penetrance may remain unresolved.",
        ),
    )
