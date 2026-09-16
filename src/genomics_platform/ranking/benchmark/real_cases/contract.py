"""Truth-blind contracts for public or de-identified benchmark cases."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


REAL_CASE_SCHEMA_VERSION = "0.1"


@dataclass(frozen=True)
class CaseProvenance:
    """Provenance for a public or de-identified benchmark case."""

    source_name: str
    source_reference: str
    source_type: str
    accessed_date: str | None = None
    notes: str | None = None

    def __post_init__(self) -> None:
        source_name = self.source_name.strip()
        source_reference = self.source_reference.strip()
        source_type = self.source_type.strip()

        if not source_name:
            raise ValueError(
                "source_name is required."
            )

        if not source_reference:
            raise ValueError(
                "source_reference is required."
            )

        if not source_type:
            raise ValueError(
                "source_type is required."
            )

        object.__setattr__(
            self,
            "source_name",
            source_name,
        )
        object.__setattr__(
            self,
            "source_reference",
            source_reference,
        )
        object.__setattr__(
            self,
            "source_type",
            source_type,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "source_name": self.source_name,
            "source_reference": self.source_reference,
            "source_type": self.source_type,
            "accessed_date": self.accessed_date,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class RealCaseInput:
    """Truth-blind ranking input for one benchmark case.

    This object may contain only information available to the ranking
    workflow. Benchmark truth is intentionally stored elsewhere.
    """

    case_id: str
    assembly: str
    present_hpo_terms: tuple[str, ...]
    absent_hpo_terms: tuple[str, ...]
    candidates: tuple[Any, ...]
    provenance: CaseProvenance

    def __post_init__(self) -> None:
        case_id = self.case_id.strip()
        assembly = self.assembly.strip()

        if not case_id:
            raise ValueError(
                "case_id is required."
            )

        if not assembly:
            raise ValueError(
                "assembly is required."
            )

        if not self.candidates:
            raise ValueError(
                "Real case must contain candidates."
            )

        present = tuple(
            dict.fromkeys(
                term.strip()
                for term in self.present_hpo_terms
                if term.strip()
            )
        )

        absent = tuple(
            dict.fromkeys(
                term.strip()
                for term in self.absent_hpo_terms
                if term.strip()
            )
        )

        overlap = (
            set(present)
            & set(absent)
        )

        if overlap:
            raise ValueError(
                "HPO terms cannot be both present "
                "and absent: "
                f"{sorted(overlap)!r}"
            )

        candidate_ids = tuple(
            candidate.candidate_id
            for candidate in self.candidates
        )

        if len(candidate_ids) != len(
            set(candidate_ids)
        ):
            raise ValueError(
                "Candidate IDs must be unique "
                "within a real case."
            )

        object.__setattr__(
            self,
            "case_id",
            case_id,
        )
        object.__setattr__(
            self,
            "assembly",
            assembly,
        )
        object.__setattr__(
            self,
            "present_hpo_terms",
            present,
        )
        object.__setattr__(
            self,
            "absent_hpo_terms",
            absent,
        )

    def ranking_candidates(self) -> tuple[Any, ...]:
        """Return candidates without benchmark truth."""
        return self.candidates

    def to_dict(self) -> dict[str, Any]:
        """Serialize metadata without candidate internals or truth."""
        return {
            "schema_version": REAL_CASE_SCHEMA_VERSION,
            "case_id": self.case_id,
            "assembly": self.assembly,
            "present_hpo_terms": list(
                self.present_hpo_terms
            ),
            "absent_hpo_terms": list(
                self.absent_hpo_terms
            ),
            "candidate_ids": [
                candidate.candidate_id
                for candidate in self.candidates
            ],
            "provenance": (
                self.provenance.to_dict()
            ),
            "semantic_boundaries": {
                "synthetic_controlled_case": False,
                "public_or_deidentified_case": True,
                "clinical_validation": False,
                "truth_present_in_ranking_input": False,
                "truth_used_for_ranking": False,
                "truth_used_for_benchmark": True,
                "pathogenicity_probability": False,
                "diagnosis": False,
            },
        }
