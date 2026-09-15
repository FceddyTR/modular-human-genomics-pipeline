#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import List, Optional


SCHEMA_VERSION = "1.1"


@dataclass
class ClinVarCondition:
    name: Optional[str] = None
    identifiers: List[str] = field(default_factory=list)


@dataclass
class MolecularConsequence:
    so_id: Optional[str] = None
    consequence: Optional[str] = None


@dataclass
class ClinVarRecord:
    assembly: str
    chrom: str
    pos: int
    ref: str
    alt: str

    variation_id: Optional[str] = None
    allele_id: Optional[str] = None
    vcv_accession: Optional[str] = None

    gene_symbols: List[str] = field(default_factory=list)
    gene_ids: List[str] = field(default_factory=list)

    clinical_significance: List[str] = field(default_factory=list)
    conflicting_significance: List[str] = field(default_factory=list)

    review_status: Optional[str] = None
    review_stars: Optional[int] = None

    conditions: List[ClinVarCondition] = field(default_factory=list)

    hgvs: List[str] = field(default_factory=list)
    molecular_consequences: List[MolecularConsequence] = field(
        default_factory=list
    )
    scv_accessions: List[str] = field(default_factory=list)

    variant_type: Optional[str] = None
    variant_type_so: Optional[str] = None
    origin: Optional[str] = None

    source: str = "ClinVar"
    source_release: Optional[str] = None
    source_reference: str = "GRCh38"
    schema_version: str = SCHEMA_VERSION

    @property
    def variant_key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    def to_dict(self):
        payload = asdict(self)
        payload["variant_key"] = self.variant_key
        return payload
