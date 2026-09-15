from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import List, Optional


SCHEMA_VERSION = "1.1"


@dataclass
class EvidenceRecord:
    """
    Normalized evidence record.

    This structure stores evidence and provenance only.
    It does NOT perform ACMG/AMP classification and does NOT establish
    variant pathogenicity.
    """

    source: str
    evidence_type: str

    gene_symbol: Optional[str] = None
    gene_id: Optional[str] = None

    disease_id: Optional[str] = None
    disease_name: Optional[str] = None

    phenotype_ids: List[str] = field(default_factory=list)

    inheritance: List[str] = field(default_factory=list)

    classification: Optional[str] = None
    confidence: Optional[str] = None

    # Frequency is distinct from evidence confidence.
    # Examples:
    #   3/3
    #   HP:0040281
    #   HP:0040283
    frequency: Optional[str] = None

    source_record_id: Optional[str] = None
    source_version: Optional[str] = None
    source_release: Optional[str] = None

    citations: List[str] = field(default_factory=list)

    evidence_description: Optional[str] = None

    license: Optional[str] = None

    def to_dict(self):
        return {
            "schema_version": SCHEMA_VERSION,
            **asdict(self),
        }
