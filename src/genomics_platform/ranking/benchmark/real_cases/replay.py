"""Replay frozen production evidence snapshots.

Only canonical constructor fields are reconstructed. Derived serialization
fields are intentionally ignored and recomputed by production contracts.
"""

from __future__ import annotations

from typing import Any

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


def _mapping(
    value: Any,
    *,
    field_name: str,
) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(
            f"{field_name} must be an object."
        )
    return value


def _sequence(
    value: Any,
    *,
    field_name: str,
) -> list[Any] | tuple[Any, ...]:
    if not isinstance(
        value,
        (list, tuple),
    ):
        raise ValueError(
            f"{field_name} must be an array or tuple."
        )
    return value


def _tuple_strings(
    value: Any,
    *,
    field_name: str,
) -> tuple[str, ...]:
    values = _sequence(
        value,
        field_name=field_name,
    )

    if not all(
        isinstance(item, str)
        for item in values
    ):
        raise ValueError(
            f"{field_name} must contain strings."
        )

    return tuple(values)


def _transcript(
    payload: dict[str, Any],
) -> TranscriptEvidence:
    return TranscriptEvidence(
        gene_symbol=payload.get("gene_symbol"),
        gene_id=payload.get("gene_id"),
        transcript_id=payload.get("transcript_id"),
        protein_id=payload.get("protein_id"),
        feature_type=payload.get("feature_type"),
        biotype=payload.get("biotype"),
        consequences=_tuple_strings(
            payload.get("consequences", []),
            field_name="transcript.consequences",
        ),
        impact=payload.get("impact"),
        canonical=payload.get(
            "canonical",
            False,
        ),
        variant_class=payload.get(
            "variant_class"
        ),
    )


def _functional(
    payload: dict[str, Any],
) -> FunctionalEvidence:
    transcripts = _sequence(
        payload.get("transcripts", []),
        field_name="functional.transcripts",
    )

    return FunctionalEvidence(
        genes=_tuple_strings(
            payload.get("genes", []),
            field_name="functional.genes",
        ),
        gene_ids=_tuple_strings(
            payload.get("gene_ids", []),
            field_name="functional.gene_ids",
        ),
        consequences=_tuple_strings(
            payload.get("consequences", []),
            field_name="functional.consequences",
        ),
        impacts=_tuple_strings(
            payload.get("impacts", []),
            field_name="functional.impacts",
        ),
        transcripts=tuple(
            _transcript(
                _mapping(
                    item,
                    field_name="functional.transcript",
                )
            )
            for item in transcripts
        ),
    )


def _condition(
    payload: dict[str, Any],
) -> ClinicalConditionEvidence:
    return ClinicalConditionEvidence(
        name=payload.get("name"),
        identifiers=_tuple_strings(
            payload.get("identifiers", []),
            field_name="condition.identifiers",
        ),
    )


def _clinical_record(
    payload: dict[str, Any],
) -> ClinicalRecordEvidence:
    conditions = _sequence(
        payload.get("conditions", []),
        field_name="clinical_record.conditions",
    )

    return ClinicalRecordEvidence(
        variation_id=payload.get(
            "variation_id"
        ),
        allele_id=payload.get(
            "allele_id"
        ),
        vcv_accession=payload.get(
            "vcv_accession"
        ),
        clinical_significance=_tuple_strings(
            payload.get(
                "clinical_significance",
                [],
            ),
            field_name=(
                "clinical_record."
                "clinical_significance"
            ),
        ),
        conflicting_significance=_tuple_strings(
            payload.get(
                "conflicting_significance",
                [],
            ),
            field_name=(
                "clinical_record."
                "conflicting_significance"
            ),
        ),
        review_status=payload.get(
            "review_status"
        ),
        review_stars=payload.get(
            "review_stars"
        ),
        conditions=tuple(
            _condition(
                _mapping(
                    item,
                    field_name=(
                        "clinical_record.condition"
                    ),
                )
            )
            for item in conditions
        ),
        hgvs=_tuple_strings(
            payload.get("hgvs", []),
            field_name="clinical_record.hgvs",
        ),
        scv_accessions=_tuple_strings(
            payload.get(
                "scv_accessions",
                [],
            ),
            field_name=(
                "clinical_record.scv_accessions"
            ),
        ),
    )


def _clinical(
    payload: dict[str, Any],
) -> ClinicalEvidence:
    records = _sequence(
        payload.get("records", []),
        field_name="clinical.records",
    )

    return ClinicalEvidence(
        records=tuple(
            _clinical_record(
                _mapping(
                    item,
                    field_name="clinical.record",
                )
            )
            for item in records
        )
    )


def _population_observation(
    payload: dict[str, Any],
) -> PopulationObservation:
    population_id = payload.get(
        "population_id"
    )

    if not isinstance(
        population_id,
        str,
    ) or not population_id.strip():
        raise ValueError(
            "population_id is required."
        )

    return PopulationObservation(
        population_id=population_id,
        ac=payload.get("ac"),
        an=payload.get("an"),
        af=payload.get("af"),
        homozygote_count=payload.get(
            "homozygote_count"
        ),
    )


def _population_dataset(
    payload: dict[str, Any],
) -> PopulationDatasetEvidence:
    populations = _sequence(
        payload.get("populations", []),
        field_name="population.populations",
    )

    return PopulationDatasetEvidence(
        ac=payload.get("ac"),
        an=payload.get("an"),
        af=payload.get("af"),
        homozygote_count=payload.get(
            "homozygote_count"
        ),
        populations=tuple(
            _population_observation(
                _mapping(
                    item,
                    field_name=(
                        "population.observation"
                    ),
                )
            )
            for item in populations
        ),
    )


def _population(
    payload: dict[str, Any],
) -> PopulationEvidence:
    exome = payload.get("exome")
    genome = payload.get("genome")

    return PopulationEvidence(
        exome=(
            _population_dataset(
                _mapping(
                    exome,
                    field_name="population.exome",
                )
            )
            if exome is not None
            else None
        ),
        genome=(
            _population_dataset(
                _mapping(
                    genome,
                    field_name="population.genome",
                )
            )
            if genome is not None
            else None
        ),
    )


def _availability(
    payload: Any,
) -> tuple[SourceAvailability, ...]:
    if isinstance(payload, dict):
        values = list(
            payload.values()
        )
    elif isinstance(payload, list):
        values = payload
    else:
        raise ValueError(
            "availability must be an object or array."
        )

    return tuple(
        SourceAvailability(
            source=_mapping(
                item,
                field_name="availability entry",
            ).get("source", ""),
            status=item.get("status", ""),
            error=item.get("error"),
        )
        for item in values
    )


def _provenance(
    payload: Any,
) -> tuple[EvidenceProvenance, ...]:
    values = _sequence(
        payload,
        field_name="provenance",
    )

    output = []

    for item in values:
        item = _mapping(
            item,
            field_name="provenance entry",
        )

        metadata = item.get(
            "metadata",
            {},
        )

        if not isinstance(
            metadata,
            dict,
        ):
            raise ValueError(
                "provenance metadata must "
                "be an object."
            )

        output.append(
            EvidenceProvenance(
                source=item.get(
                    "source",
                    "",
                ),
                source_release=item.get(
                    "source_release"
                ),
                source_dataset=item.get(
                    "source_dataset"
                ),
                access_method=item.get(
                    "access_method"
                ),
                metadata=metadata,
            )
        )

    return tuple(output)


def replay_variant_evidence_contract(
    payload: dict[str, Any],
) -> VariantEvidenceContract:
    """Reconstruct a production VariantEvidenceContract."""

    payload = _mapping(
        payload,
        field_name="evidence_contract",
    )

    identity = _mapping(
        payload.get("identity"),
        field_name="identity",
    )

    contract = VariantEvidenceContract(
        identity=VariantIdentity(
            assembly=identity.get(
                "assembly",
                "",
            ),
            chrom=identity.get(
                "chrom",
                "",
            ),
            pos=identity.get("pos"),
            ref=identity.get(
                "ref",
                "",
            ),
            alt=identity.get(
                "alt",
                "",
            ),
        ),
        functional=_functional(
            _mapping(
                payload.get("functional"),
                field_name="functional",
            )
        ),
        clinical=_clinical(
            _mapping(
                payload.get("clinical"),
                field_name="clinical",
            )
        ),
        population=_population(
            _mapping(
                payload.get("population"),
                field_name="population",
            )
        ),
        availability=_availability(
            payload.get(
                "availability",
                {},
            )
        ),
        provenance=_provenance(
            payload.get(
                "provenance",
                [],
            )
        ),
        schema_version=str(
            payload.get(
                "schema_version",
                "1.0",
            )
        ),
    )

    serialized_key = payload.get(
        "variant_key"
    )

    if (
        serialized_key is not None
        and serialized_key
        != contract.identity.variant_key
    ):
        raise ValueError(
            "Serialized variant_key does not "
            "match reconstructed identity."
        )

    return contract


# ---------------------------------------------------------------------------
# Interpretation and CandidateCase replay
# ---------------------------------------------------------------------------

from genomics_platform.interpretation.case_model import (
    CandidateCase,
    CaseContext,
)
from genomics_platform.interpretation.interpretation_profile import (
    EvidenceDimension,
    InterpretationObservation,
    InterpretationProfile,
)


def _interpretation_observation(
    payload: dict[str, Any],
) -> InterpretationObservation:
    payload = _mapping(
        payload,
        field_name="interpretation observation",
    )

    details = payload.get(
        "details",
        {},
    )

    if not isinstance(
        details,
        dict,
    ):
        raise ValueError(
            "interpretation observation details "
            "must be an object."
        )

    code = payload.get("code")
    label = payload.get("label")

    if not isinstance(
        code,
        str,
    ) or not code.strip():
        raise ValueError(
            "interpretation observation code "
            "is required."
        )

    if not isinstance(
        label,
        str,
    ) or not label.strip():
        raise ValueError(
            "interpretation observation label "
            "is required."
        )

    return InterpretationObservation(
        code=code,
        label=label,
        value=payload.get("value"),
        source=payload.get("source"),
        details=details,
    )


def _evidence_dimension(
    payload: dict[str, Any],
) -> EvidenceDimension:
    payload = _mapping(
        payload,
        field_name="interpretation dimension",
    )

    observations = _sequence(
        payload.get(
            "observations",
            [],
        ),
        field_name=(
            "interpretation dimension observations"
        ),
    )

    return EvidenceDimension(
        name=payload.get(
            "name",
            "",
        ),
        status=payload.get(
            "status",
            "",
        ),
        observations=tuple(
            _interpretation_observation(
                _mapping(
                    item,
                    field_name=(
                        "interpretation observation"
                    ),
                )
            )
            for item in observations
        ),
        limitations=_tuple_strings(
            payload.get(
                "limitations",
                [],
            ),
            field_name=(
                "interpretation dimension limitations"
            ),
        ),
    )


def replay_interpretation_profile(
    payload: dict[str, Any],
) -> InterpretationProfile:
    """Reconstruct a production InterpretationProfile."""

    payload = _mapping(
        payload,
        field_name="interpretation_profile",
    )

    dimensions = _mapping(
        payload.get("dimensions"),
        field_name=(
            "interpretation_profile.dimensions"
        ),
    )

    def optional_dimension(
        name: str,
    ) -> EvidenceDimension | None:
        value = dimensions.get(name)

        if value is None:
            return None

        return _evidence_dimension(
            _mapping(
                value,
                field_name=(
                    f"interpretation dimension {name}"
                ),
            )
        )

    for required in (
        "population",
        "clinical",
        "functional",
    ):
        if dimensions.get(required) is None:
            raise ValueError(
                "Required interpretation dimension "
                f"is missing: {required}"
            )

    return InterpretationProfile(
        variant_key=str(
            payload.get(
                "variant_key",
                "",
            )
        ),
        population=_evidence_dimension(
            _mapping(
                dimensions["population"],
                field_name=(
                    "interpretation dimension population"
                ),
            )
        ),
        clinical=_evidence_dimension(
            _mapping(
                dimensions["clinical"],
                field_name=(
                    "interpretation dimension clinical"
                ),
            )
        ),
        functional=_evidence_dimension(
            _mapping(
                dimensions["functional"],
                field_name=(
                    "interpretation dimension functional"
                ),
            )
        ),
        phenotype=optional_dimension(
            "phenotype"
        ),
        inheritance=optional_dimension(
            "inheritance"
        ),
        gene_disease=optional_dimension(
            "gene_disease"
        ),
        schema_version=str(
            payload.get(
                "schema_version",
                "1.0",
            )
        ),
    )


def replay_case_context(
    payload: dict[str, Any],
) -> CaseContext:
    """Reconstruct production patient/case context."""

    payload = _mapping(
        payload,
        field_name="case_context",
    )

    return CaseContext(
        present_hpo_terms=_tuple_strings(
            payload.get(
                "present_hpo_terms",
                [],
            ),
            field_name=(
                "case_context.present_hpo_terms"
            ),
        ),
        absent_hpo_terms=_tuple_strings(
            payload.get(
                "absent_hpo_terms",
                [],
            ),
            field_name=(
                "case_context.absent_hpo_terms"
            ),
        ),
        proband_sex=payload.get(
            "proband_sex"
        ),
        genotype_state=payload.get(
            "genotype_state"
        ),
        de_novo_status=payload.get(
            "de_novo_status"
        ),
        family_context=payload.get(
            "family_context"
        ),
    )


def replay_candidate_case(
    payload: dict[str, Any],
) -> CandidateCase:
    """Reconstruct one frozen production CandidateCase."""

    payload = _mapping(
        payload,
        field_name="candidate",
    )

    evidence_payload = _mapping(
        payload.get("evidence_contract"),
        field_name="candidate.evidence_contract",
    )

    profile_payload = _mapping(
        payload.get("interpretation_profile"),
        field_name=(
            "candidate.interpretation_profile"
        ),
    )

    context_payload = _mapping(
        payload.get(
            "case_context",
            {},
        ),
        field_name="candidate.case_context",
    )

    candidate = CandidateCase(
        candidate_id=str(
            payload.get(
                "candidate_id",
                "",
            )
        ),
        evidence_contract=(
            replay_variant_evidence_contract(
                evidence_payload
            )
        ),
        interpretation_profile=(
            replay_interpretation_profile(
                profile_payload
            )
        ),
        case_context=replay_case_context(
            context_payload
        ),
        gene_symbol=payload.get(
            "gene_symbol"
        ),
        schema_version=str(
            payload.get(
                "schema_version",
                "1.0",
            )
        ),
    )

    serialized_key = payload.get(
        "variant_key"
    )

    if (
        serialized_key is not None
        and serialized_key
        != candidate.variant_key
    ):
        raise ValueError(
            "Serialized candidate variant_key "
            "does not match reconstructed candidate."
        )

    return candidate
