from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional


SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class EvidenceComponent:
    source: str
    status: str
    data: Optional[Any] = None
    error: Optional[str] = None
    provenance: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "source": self.source,
            "status": self.status,
            "data": self.data,
            "error": self.error,
            "provenance": self.provenance,
        }


@dataclass(frozen=True)
class VariantEvidenceBundle:
    assembly: str
    chrom: str
    pos: int
    ref: str
    alt: str

    vep: EvidenceComponent
    clinvar: EvidenceComponent
    gnomad: EvidenceComponent

    schema_version: str = SCHEMA_VERSION

    @property
    def variant_key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @property
    def evidence_complete(self) -> bool:
        """
        Technical completeness only.

        This does NOT mean that the variant has sufficient clinical evidence.
        NOT_FOUND is a valid completed lookup; ERROR is not.
        """
        return all(
            component.status != "ERROR"
            for component in (self.vep, self.clinvar, self.gnomad)
        )

    @property
    def error_sources(self) -> list[str]:
        return [
            component.source
            for component in (self.vep, self.clinvar, self.gnomad)
            if component.status == "ERROR"
        ]

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "assembly": self.assembly,
            "variant_key": self.variant_key,
            "variant": {
                "chrom": self.chrom,
                "pos": self.pos,
                "ref": self.ref,
                "alt": self.alt,
            },
            "evidence_complete": self.evidence_complete,
            "error_sources": self.error_sources,
            "evidence": {
                "functional": self.vep.to_dict(),
                "clinical": self.clinvar.to_dict(),
                "population": self.gnomad.to_dict(),
            },
            "semantics": {
                "bundle": (
                    "This object aggregates evidence and does not perform "
                    "clinical classification or candidate ranking."
                ),
                "not_found": (
                    "NOT_FOUND means a source lookup completed without a matching "
                    "record; it is not equivalent to benign, rare, or AF=0."
                ),
                "error": (
                    "ERROR represents unavailable or failed evidence and must not "
                    "be interpreted as NOT_FOUND."
                ),
                "completeness": (
                    "evidence_complete describes technical lookup completion only, "
                    "not clinical evidence sufficiency."
                ),
            },
        }
