from __future__ import annotations

from typing import Any, Iterable, Optional

from genomics_platform.evidence.variant_evidence_bundle import (
    EvidenceComponent,
    VariantEvidenceBundle,
)


def canonical_chrom(chrom: str) -> str:
    value = str(chrom).strip()

    if value.lower().startswith("chr"):
        value = value[3:]

    if value.upper() == "M":
        return "MT"

    return value


def _serialize(value: Any) -> Any:
    if value is None:
        return None

    if hasattr(value, "to_dict"):
        return value.to_dict()

    return value


def _clinvar_component(result: Any) -> EvidenceComponent:
    """
    Accept either:

    1. Native ClinVar store result:
       {
           "status": "FOUND" | "NOT_FOUND",
           "interpretation": ...,
           "query": {...},
           "records": [...]
       }

    2. A plain iterable of records, retained for synthetic/unit tests.

    Missing input is ERROR, never NOT_FOUND.
    """

    if result is None:
        return EvidenceComponent(
            source="ClinVar",
            status="ERROR",
            error="ClinVar lookup result was not supplied.",
        )

    # Native clinvar_store.query_variant() contract.
    if isinstance(result, dict) and "status" in result and "records" in result:
        status = str(result.get("status", "ERROR")).upper()
        records = result.get("records") or []

        provenance = {
            "interpretation": result.get("interpretation"),
            "query": result.get("query"),
            "record_count": len(records),
        }

        releases = sorted(
            {
                str(record.get("source_release"))
                for record in records
                if isinstance(record, dict)
                and record.get("source_release")
            }
        )

        if releases:
            provenance["releases"] = releases

        if status == "FOUND":
            return EvidenceComponent(
                source="ClinVar",
                status="FOUND",
                data=[_serialize(record) for record in records],
                provenance=provenance,
            )

        if status == "NOT_FOUND":
            return EvidenceComponent(
                source="ClinVar",
                status="NOT_FOUND",
                data=[],
                provenance=provenance,
            )

        return EvidenceComponent(
            source="ClinVar",
            status="ERROR",
            data=[_serialize(record) for record in records],
            error=(
                f"Unexpected ClinVar lookup status: {status}"
            ),
            provenance=provenance,
        )

    # Synthetic/unit-test compatibility.
    try:
        records = list(result)
    except TypeError:
        return EvidenceComponent(
            source="ClinVar",
            status="ERROR",
            error=(
                "Unsupported ClinVar result type: "
                f"{type(result).__name__}"
            ),
        )

    if not records:
        return EvidenceComponent(
            source="ClinVar",
            status="NOT_FOUND",
            data=[],
            provenance={
                "record_count": 0,
                "input_mode": "record_iterable",
            },
        )

    serialized = [_serialize(record) for record in records]

    releases = sorted(
        {
            str(record.get("source_release"))
            for record in serialized
            if isinstance(record, dict)
            and record.get("source_release")
        }
    )

    return EvidenceComponent(
        source="ClinVar",
        status="FOUND",
        data=serialized,
        provenance={
            "record_count": len(serialized),
            "releases": releases,
            "input_mode": "record_iterable",
        },
    )


def _gnomad_component(record: Any) -> EvidenceComponent:
    if record is None:
        return EvidenceComponent(
            source="gnomAD",
            status="ERROR",
            error="gnomAD lookup result was not supplied.",
        )

    status = str(
        getattr(record, "lookup_status", "ERROR")
    ).upper()

    error = getattr(record, "error", None)

    if status not in {"FOUND", "NOT_FOUND", "ERROR"}:
        error = f"Unexpected gnomAD lookup status: {status}"
        status = "ERROR"

    provenance = {
        "dataset": getattr(record, "source_dataset", None),
        "release": getattr(record, "source_release", None),
        "access_method": getattr(record, "access_method", None),
    }

    return EvidenceComponent(
        source="gnomAD",
        status=status,
        data=_serialize(record),
        error=error,
        provenance=provenance,
    )


def _vep_component(
    annotations: Optional[Iterable[Any]],
    *,
    source_release: Optional[str] = None,
    provenance: Optional[dict[str, Any]] = None,
) -> EvidenceComponent:

    base_provenance = {
        "release": source_release,
    }

    if provenance:
        base_provenance.update(provenance)

    if annotations is None:
        return EvidenceComponent(
            source="VEP",
            status="ERROR",
            error="VEP annotation result was not supplied.",
            provenance=base_provenance,
        )

    annotations = list(annotations)

    if not annotations:
        return EvidenceComponent(
            source="VEP",
            status="NOT_FOUND",
            data=[],
            provenance=base_provenance,
        )

    return EvidenceComponent(
        source="VEP",
        status="FOUND",
        data=[_serialize(annotation) for annotation in annotations],
        provenance=base_provenance,
    )


def build_variant_evidence_bundle(
    *,
    assembly: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    vep_annotations: Optional[Iterable[Any]],
    clinvar_records: Any,
    gnomad_record: Any,
    vep_release: Optional[str] = None,
    vep_provenance: Optional[dict[str, Any]] = None,
) -> VariantEvidenceBundle:

    chrom = canonical_chrom(chrom)
    pos = int(pos)
    ref = str(ref).upper()
    alt = str(alt).upper()

    return VariantEvidenceBundle(
        assembly=assembly,
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        vep=_vep_component(
            vep_annotations,
            source_release=vep_release,
            provenance=vep_provenance,
        ),
        clinvar=_clinvar_component(clinvar_records),
        gnomad=_gnomad_component(gnomad_record),
    )
