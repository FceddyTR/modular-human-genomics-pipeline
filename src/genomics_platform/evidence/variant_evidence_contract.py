from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Optional


SCHEMA_VERSION = "1.0"

VALID_AVAILABILITY_STATUSES = {"FOUND", "NOT_FOUND", "ERROR"}


@dataclass(frozen=True)
class VariantIdentity:
    assembly: str
    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def variant_key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    def to_dict(self) -> dict[str, Any]:
        return {
            "assembly": self.assembly,
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "variant_key": self.variant_key,
        }


@dataclass(frozen=True)
class TranscriptEvidence:
    gene_symbol: Optional[str] = None
    gene_id: Optional[str] = None
    transcript_id: Optional[str] = None
    protein_id: Optional[str] = None
    feature_type: Optional[str] = None
    biotype: Optional[str] = None
    consequences: tuple[str, ...] = field(default_factory=tuple)
    impact: Optional[str] = None
    canonical: bool = False
    variant_class: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class FunctionalEvidence:
    genes: tuple[str, ...] = field(default_factory=tuple)
    gene_ids: tuple[str, ...] = field(default_factory=tuple)
    consequences: tuple[str, ...] = field(default_factory=tuple)
    impacts: tuple[str, ...] = field(default_factory=tuple)
    transcripts: tuple[TranscriptEvidence, ...] = field(default_factory=tuple)

    def to_dict(self) -> dict[str, Any]:
        return {
            "genes": list(self.genes),
            "gene_ids": list(self.gene_ids),
            "consequences": list(self.consequences),
            "impacts": list(self.impacts),
            "canonical_transcript_count": sum(
                1 for transcript in self.transcripts if transcript.canonical
            ),
            "transcripts": [
                transcript.to_dict() for transcript in self.transcripts
            ],
        }


@dataclass(frozen=True)
class ClinicalConditionEvidence:
    name: Optional[str] = None
    identifiers: tuple[str, ...] = field(default_factory=tuple)

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "identifiers": list(self.identifiers),
        }


@dataclass(frozen=True)
class ClinicalRecordEvidence:
    variation_id: Optional[str] = None
    allele_id: Optional[str] = None
    vcv_accession: Optional[str] = None
    clinical_significance: tuple[str, ...] = field(default_factory=tuple)
    conflicting_significance: tuple[str, ...] = field(default_factory=tuple)
    review_status: Optional[str] = None
    review_stars: Optional[int] = None
    conditions: tuple[ClinicalConditionEvidence, ...] = field(
        default_factory=tuple
    )
    hgvs: tuple[str, ...] = field(default_factory=tuple)
    scv_accessions: tuple[str, ...] = field(default_factory=tuple)

    def to_dict(self) -> dict[str, Any]:
        return {
            "variation_id": self.variation_id,
            "allele_id": self.allele_id,
            "vcv_accession": self.vcv_accession,
            "clinical_significance": list(self.clinical_significance),
            "conflicting_significance": list(self.conflicting_significance),
            "review_status": self.review_status,
            "review_stars": self.review_stars,
            "conditions": [
                condition.to_dict() for condition in self.conditions
            ],
            "hgvs": list(self.hgvs),
            "scv_accessions": list(self.scv_accessions),
        }


@dataclass(frozen=True)
class ClinicalEvidence:
    records: tuple[ClinicalRecordEvidence, ...] = field(default_factory=tuple)

    def to_dict(self) -> dict[str, Any]:
        return {
            "record_count": len(self.records),
            "records": [record.to_dict() for record in self.records],
        }


@dataclass(frozen=True)
class PopulationObservation:
    population_id: str
    ac: Optional[int] = None
    an: Optional[int] = None
    af: Optional[float] = None
    homozygote_count: Optional[int] = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class PopulationDatasetEvidence:
    ac: Optional[int] = None
    an: Optional[int] = None
    af: Optional[float] = None
    homozygote_count: Optional[int] = None
    populations: tuple[PopulationObservation, ...] = field(default_factory=tuple)

    def to_dict(self) -> dict[str, Any]:
        return {
            "ac": self.ac,
            "an": self.an,
            "af": self.af,
            "homozygote_count": self.homozygote_count,
            "populations": [
                population.to_dict() for population in self.populations
            ],
        }


@dataclass(frozen=True)
class PopulationEvidence:
    exome: Optional[PopulationDatasetEvidence] = None
    genome: Optional[PopulationDatasetEvidence] = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "exome": self.exome.to_dict() if self.exome else None,
            "genome": self.genome.to_dict() if self.genome else None,
        }


@dataclass(frozen=True)
class SourceAvailability:
    source: str
    status: str
    error: Optional[str] = None

    def __post_init__(self) -> None:
        normalized_status = self.status.upper()

        if normalized_status not in VALID_AVAILABILITY_STATUSES:
            raise ValueError(
                f"Invalid availability status for {self.source}: "
                f"{self.status}"
            )

        object.__setattr__(self, "status", normalized_status)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class EvidenceProvenance:
    source: str
    source_release: Optional[str] = None
    source_dataset: Optional[str] = None
    access_method: Optional[str] = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class VariantEvidenceContract:
    identity: VariantIdentity
    functional: FunctionalEvidence
    clinical: ClinicalEvidence
    population: PopulationEvidence
    availability: tuple[SourceAvailability, ...]
    provenance: tuple[EvidenceProvenance, ...] = field(default_factory=tuple)
    schema_version: str = SCHEMA_VERSION

    @property
    def evidence_complete(self) -> bool:
        return all(
            availability.status != "ERROR"
            for availability in self.availability
        )

    @property
    def error_sources(self) -> list[str]:
        return [
            availability.source
            for availability in self.availability
            if availability.status == "ERROR"
        ]

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "variant_key": self.identity.variant_key,
            "identity": self.identity.to_dict(),
            "functional": self.functional.to_dict(),
            "clinical": self.clinical.to_dict(),
            "population": self.population.to_dict(),
            "availability": {
                availability.source: availability.to_dict()
                for availability in self.availability
            },
            "provenance": [
                provenance.to_dict() for provenance in self.provenance
            ],
            "evidence_complete": self.evidence_complete,
            "error_sources": self.error_sources,
            "semantics": {
                "contract": (
                    "Normalized evidence contract for downstream "
                    "interpretation. This contract does not perform "
                    "clinical classification or candidate ranking."
                ),
                "not_found": (
                    "NOT_FOUND means that a source lookup completed "
                    "without a matching record. It does not imply benignity, "
                    "AF=0, rarity, or absence of biological relevance."
                ),
                "error": (
                    "ERROR means that evidence was unavailable or the "
                    "lookup failed. ERROR must not be interpreted as "
                    "NOT_FOUND."
                ),
                "population": (
                    "Population frequency is evidence only. Exome and "
                    "genome observations remain separate and no "
                    "pathogenicity conclusion is derived here."
                ),
                "clinical": (
                    "Clinical database records are preserved as evidence. "
                    "Review level and submitted significance are not "
                    "converted into an ACMG/AMP classification here."
                ),
                "functional": (
                    "Functional consequences describe predicted molecular "
                    "effects and do not establish pathogenicity."
                ),
            },
        }
