from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field
from typing import List

from .evidence_record import EvidenceRecord


BUNDLE_SCHEMA_VERSION = "1.0"


@dataclass
class EvidenceBundle:
    gene_symbol: str
    records: List[EvidenceRecord] = field(default_factory=list)

    def add(self, record: EvidenceRecord):
        self.records.append(record)

    def summary(self):
        sources = Counter(record.source for record in self.records)
        evidence_types = Counter(record.evidence_type for record in self.records)

        diseases = sorted(
            {
                record.disease_id
                for record in self.records
                if record.disease_id
            }
        )

        phenotypes = sorted(
            {
                phenotype
                for record in self.records
                for phenotype in record.phenotype_ids
            }
        )

        return {
            "gene_symbol": self.gene_symbol,
            "record_count": len(self.records),
            "sources": dict(sorted(sources.items())),
            "evidence_types": dict(sorted(evidence_types.items())),
            "disease_ids": diseases,
            "phenotype_ids": phenotypes,
        }

    def to_dict(self):
        return {
            "schema_version": BUNDLE_SCHEMA_VERSION,
            "scope": {
                "clinical_classification": False,
                "candidate_ranking": False,
                "statement": (
                    "Evidence aggregation only. "
                    "This structure does not establish pathogenicity."
                ),
            },
            "gene_symbol": self.gene_symbol,
            "summary": self.summary(),
            "records": [record.to_dict() for record in self.records],
        }
