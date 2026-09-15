from __future__ import annotations

from typing import Any, Iterable, Optional

from genomics_platform.evidence.variant_evidence_bundle import (
    EvidenceComponent,
    VariantEvidenceBundle,
)
from genomics_platform.evidence.variant_evidence_contract import (
    ClinicalConditionEvidence,
    ClinicalEvidence,
    ClinicalRecordEvidence,
    EvidenceProvenance,
    FunctionalEvidence,
    PopulationDatasetEvidence,
    PopulationEvidence,
    PopulationObservation,
    SourceAvailability,
    TranscriptEvidence,
    VariantEvidenceContract,
    VariantIdentity,
)


def _clean_optional(value: Any) -> Optional[str]:
    if value is None:
        return None

    text = str(value).strip()

    return text if text else None


def _as_tuple(value: Any) -> tuple:
    if value is None:
        return ()

    if isinstance(value, tuple):
        return value

    if isinstance(value, list):
        return tuple(value)

    return (value,)


def _unique_sorted(values: Iterable[Any]) -> tuple[str, ...]:
    cleaned = {
        str(value).strip()
        for value in values
        if value is not None and str(value).strip()
    }

    return tuple(sorted(cleaned))


def _extract_vep_annotations(component: EvidenceComponent) -> list[dict]:
    if component.status != "FOUND":
        return []

    data = component.data

    if data is None:
        return []

    if isinstance(data, list):
        return [
            item
            for item in data
            if isinstance(item, dict)
        ]

    if isinstance(data, tuple):
        return [
            item
            for item in data
            if isinstance(item, dict)
        ]

    if isinstance(data, dict):
        if isinstance(data.get("transcript_annotations"), list):
            return [
                item
                for item in data["transcript_annotations"]
                if isinstance(item, dict)
            ]

        if isinstance(data.get("variants"), list):
            annotations = []

            for variant in data["variants"]:
                if not isinstance(variant, dict):
                    continue

                for item in variant.get("transcript_annotations", []):
                    if isinstance(item, dict):
                        annotations.append(item)

            return annotations

        return [data]

    return []


def _build_functional_evidence(
    component: EvidenceComponent,
) -> FunctionalEvidence:
    annotations = _extract_vep_annotations(component)

    transcripts = []

    genes = []
    gene_ids = []
    consequences = []
    impacts = []

    for annotation in annotations:
        normalized = annotation.get("_normalized", annotation)

        if not isinstance(normalized, dict):
            continue

        gene_symbol = _clean_optional(normalized.get("symbol"))
        gene_id = _clean_optional(normalized.get("gene_id"))
        transcript_id = _clean_optional(normalized.get("transcript_id"))
        protein_id = _clean_optional(normalized.get("protein_id"))
        feature_type = _clean_optional(normalized.get("feature_type"))
        biotype = _clean_optional(normalized.get("biotype"))
        impact = _clean_optional(normalized.get("impact"))
        variant_class = _clean_optional(normalized.get("variant_class"))

        raw_consequence = _clean_optional(
            normalized.get("consequence")
        )

        transcript_consequences = ()

        if raw_consequence:
            transcript_consequences = tuple(
                term
                for term in raw_consequence.split("&")
                if term
            )

        canonical = normalized.get("canonical", False)

        if isinstance(canonical, str):
            canonical = canonical.upper() == "YES"
        else:
            canonical = bool(canonical)

        transcript = TranscriptEvidence(
            gene_symbol=gene_symbol,
            gene_id=gene_id,
            transcript_id=transcript_id,
            protein_id=protein_id,
            feature_type=feature_type,
            biotype=biotype,
            consequences=transcript_consequences,
            impact=impact,
            canonical=canonical,
            variant_class=variant_class,
        )

        transcripts.append(transcript)

        if gene_symbol:
            genes.append(gene_symbol)

        if gene_id:
            gene_ids.append(gene_id)

        consequences.extend(transcript_consequences)

        if impact:
            impacts.append(impact)

    return FunctionalEvidence(
        genes=_unique_sorted(genes),
        gene_ids=_unique_sorted(gene_ids),
        consequences=_unique_sorted(consequences),
        impacts=_unique_sorted(impacts),
        transcripts=tuple(transcripts),
    )


def _build_clinical_evidence(
    component: EvidenceComponent,
) -> ClinicalEvidence:
    if component.status != "FOUND":
        return ClinicalEvidence()

    data = component.data

    if data is None:
        return ClinicalEvidence()

    if isinstance(data, dict):
        if isinstance(data.get("records"), list):
            records = data["records"]
        else:
            records = [data]
    elif isinstance(data, (list, tuple)):
        records = list(data)
    else:
        records = []

    output_records = []

    for record in records:
        if hasattr(record, "to_dict"):
            record = record.to_dict()

        if not isinstance(record, dict):
            continue

        conditions = []

        for condition in record.get("conditions", []) or []:
            if hasattr(condition, "to_dict"):
                condition = condition.to_dict()

            if not isinstance(condition, dict):
                continue

            conditions.append(
                ClinicalConditionEvidence(
                    name=_clean_optional(condition.get("name")),
                    identifiers=_unique_sorted(
                        condition.get("identifiers", []) or []
                    ),
                )
            )

        output_records.append(
            ClinicalRecordEvidence(
                variation_id=_clean_optional(
                    record.get("variation_id")
                ),
                allele_id=_clean_optional(
                    record.get("allele_id")
                ),
                vcv_accession=_clean_optional(
                    record.get("vcv_accession")
                ),
                clinical_significance=_unique_sorted(
                    record.get("clinical_significance", []) or []
                ),
                conflicting_significance=_unique_sorted(
                    record.get("conflicting_significance", []) or []
                ),
                review_status=_clean_optional(
                    record.get("review_status")
                ),
                review_stars=record.get("review_stars"),
                conditions=tuple(conditions),
                hgvs=_unique_sorted(
                    record.get("hgvs", []) or []
                ),
                scv_accessions=_unique_sorted(
                    record.get("scv_accessions", []) or []
                ),
            )
        )

    return ClinicalEvidence(
        records=tuple(output_records)
    )


def _build_population_dataset(
    data: Any,
) -> Optional[PopulationDatasetEvidence]:
    if data is None:
        return None

    if hasattr(data, "to_dict"):
        data = data.to_dict()

    if not isinstance(data, dict):
        return None

    populations = []

    for population in data.get("populations", []) or []:
        if hasattr(population, "to_dict"):
            population = population.to_dict()

        if not isinstance(population, dict):
            continue

        population_id = (
            population.get("id")
            or population.get("population_id")
        )

        if not population_id:
            continue

        populations.append(
            PopulationObservation(
                population_id=str(population_id),
                ac=population.get("ac"),
                an=population.get("an"),
                af=population.get("af"),
                homozygote_count=population.get(
                    "homozygote_count"
                ),
            )
        )

    return PopulationDatasetEvidence(
        ac=data.get("ac"),
        an=data.get("an"),
        af=data.get("af"),
        homozygote_count=data.get("homozygote_count"),
        populations=tuple(populations),
    )


def _build_population_evidence(
    component: EvidenceComponent,
) -> PopulationEvidence:
    if component.status != "FOUND":
        return PopulationEvidence()

    data = component.data

    if hasattr(data, "to_dict"):
        data = data.to_dict()

    if not isinstance(data, dict):
        return PopulationEvidence()

    return PopulationEvidence(
        exome=_build_population_dataset(
            data.get("exome")
        ),
        genome=_build_population_dataset(
            data.get("genome")
        ),
    )


def _availability(
    component: EvidenceComponent,
) -> SourceAvailability:
    return SourceAvailability(
        source=component.source,
        status=component.status,
        error=component.error,
    )


def _provenance_from_component(
    component: EvidenceComponent,
) -> EvidenceProvenance:
    provenance = dict(component.provenance or {})

    source_release = provenance.pop(
        "source_release",
        provenance.pop("release", None),
    )

    source_dataset = provenance.pop(
        "source_dataset",
        provenance.pop("dataset", None),
    )

    access_method = provenance.pop(
        "access_method",
        None,
    )

    return EvidenceProvenance(
        source=component.source,
        source_release=_clean_optional(source_release),
        source_dataset=_clean_optional(source_dataset),
        access_method=_clean_optional(access_method),
        metadata=provenance,
    )


def build_variant_evidence_contract(
    bundle: VariantEvidenceBundle,
) -> VariantEvidenceContract:
    """
    Convert a source-oriented VariantEvidenceBundle into the stable,
    source-normalized Variant Evidence Contract.

    This transformation performs normalization only.

    It does NOT:
      - classify pathogenicity,
      - apply ACMG/AMP criteria,
      - infer rarity from NOT_FOUND,
      - rank candidates,
      - assign tiers,
      - calculate pathogenicity probabilities.
    """

    identity = VariantIdentity(
        assembly=bundle.assembly,
        chrom=bundle.chrom,
        pos=bundle.pos,
        ref=bundle.ref,
        alt=bundle.alt,
    )

    functional = _build_functional_evidence(
        bundle.vep
    )

    clinical = _build_clinical_evidence(
        bundle.clinvar
    )

    population = _build_population_evidence(
        bundle.gnomad
    )

    availability = (
        _availability(bundle.vep),
        _availability(bundle.clinvar),
        _availability(bundle.gnomad),
    )

    provenance = (
        _provenance_from_component(bundle.vep),
        _provenance_from_component(bundle.clinvar),
        _provenance_from_component(bundle.gnomad),
    )

    return VariantEvidenceContract(
        identity=identity,
        functional=functional,
        clinical=clinical,
        population=population,
        availability=availability,
        provenance=provenance,
    )
