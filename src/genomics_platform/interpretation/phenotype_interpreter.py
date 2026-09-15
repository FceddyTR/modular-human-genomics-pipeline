"""Descriptive phenotype interpretation.

This module compares patient-supplied HPO terms with HPO associations
already present in a Gene Evidence Bundle.

v0.1 intentionally performs exact HPO identifier matching only.

It does NOT:
- calculate phenotype similarity scores,
- perform ontology ancestor/descendant matching,
- rank genes or variants,
- classify variants,
- diagnose disease,
- treat absence of an exact match as disease exclusion.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Set

from genomics_platform.interpretation.interpretation_profile import (
    EvidenceDimension,
    InterpretationObservation,
)


ENGINE_VERSION = "0.1.0"


def _normalize_hpo_terms(
    terms: Optional[Iterable[str]],
) -> List[str]:
    """Normalize, validate and de-duplicate HPO identifiers."""

    if terms is None:
        return []

    if isinstance(terms, str):
        terms = [terms]

    normalized: List[str] = []
    seen: Set[str] = set()

    for term in terms:
        if not isinstance(term, str):
            raise TypeError(
                "HPO terms must be strings."
            )

        value = term.strip().upper()

        if not value:
            continue

        if not value.startswith("HP:"):
            raise ValueError(
                f"Invalid HPO identifier: {term!r}"
            )

        suffix = value[3:]

        if len(suffix) != 7 or not suffix.isdigit():
            raise ValueError(
                f"Invalid HPO identifier: {term!r}"
            )

        if value not in seen:
            seen.add(value)
            normalized.append(value)

    return normalized


def _extract_hpo_records(
    gene_evidence_bundle: Dict[str, Any],
) -> List[Dict[str, Any]]:
    """Return HPO source records from a Gene Evidence Bundle."""

    records = gene_evidence_bundle.get("records", [])

    if not isinstance(records, list):
        return []

    hpo_records: List[Dict[str, Any]] = []

    for record in records:
        if not isinstance(record, dict):
            continue

        evidence = record.get("evidence")

        if not isinstance(evidence, dict):
            continue

        if evidence.get("source") != "HPO":
            continue

        hpo_records.append(record)

    return hpo_records


def _associated_terms(
    hpo_records: Iterable[Dict[str, Any]],
) -> Set[str]:
    terms: Set[str] = set()

    for record in hpo_records:
        evidence = record.get("evidence", {})
        phenotype_ids = evidence.get("phenotype_ids", [])

        if not isinstance(phenotype_ids, list):
            continue

        for term in phenotype_ids:
            if isinstance(term, str):
                terms.add(term.strip().upper())

    return terms


def interpret_phenotype(
    gene_evidence_bundle: Optional[Dict[str, Any]],
    present_hpo_terms: Optional[Iterable[str]] = None,
    absent_hpo_terms: Optional[Iterable[str]] = None,
) -> EvidenceDimension:
    """Interpret phenotype context descriptively.

    Parameters
    ----------
    gene_evidence_bundle:
        Existing Gene Evidence Bundle.

    present_hpo_terms:
        HPO terms observed in the patient.

    absent_hpo_terms:
        HPO terms explicitly assessed and absent in the patient.

    Returns
    -------
    EvidenceDimension
        Descriptive phenotype evidence only.
    """

    present = _normalize_hpo_terms(
        present_hpo_terms
    )
    absent = _normalize_hpo_terms(
        absent_hpo_terms
    )

    overlap = set(present) & set(absent)

    if overlap:
        raise ValueError(
            "The same HPO term cannot be both present "
            f"and absent: {sorted(overlap)}"
        )

    if gene_evidence_bundle is None:
        return EvidenceDimension(
            name="phenotype",
            status="UNAVAILABLE",
            observations=[],
            limitations=[
                "No gene evidence context was supplied.",
                "Phenotype evidence was not evaluated.",
                "No phenotype ranking or diagnostic inference was performed.",
            ],
        )

    if not isinstance(gene_evidence_bundle, dict):
        raise TypeError(
            "gene_evidence_bundle must be a dictionary or None."
        )

    if not present and not absent:
        return EvidenceDimension(
            name="phenotype",
            status="UNAVAILABLE",
            observations=[],
            limitations=[
                "No patient HPO terms were supplied.",
                "Missing phenotype context does not imply phenotype mismatch.",
                "No phenotype ranking or diagnostic inference was performed.",
            ],
        )

    hpo_records = _extract_hpo_records(
        gene_evidence_bundle
    )

    associated = _associated_terms(
        hpo_records
    )

    present_set = set(present)
    absent_set = set(absent)

    matched_present = [
        term for term in present
        if term in associated
    ]

    unmatched_present = [
        term for term in present
        if term not in associated
    ]

    matched_absent = [
        term for term in absent
        if term in associated
    ]

    unmatched_absent = [
        term for term in absent
        if term not in associated
    ]

    observations: List[
        InterpretationObservation
    ] = []

    observations.append(
        InterpretationObservation(
            code="PHENOTYPE_PATIENT_CONTEXT",
            label="Patient phenotype context",
            value={
                "present_hpo_terms": present,
                "absent_hpo_terms": absent,
            },
            source="patient_context",
            details={
                "phenotype_similarity_score": False,
                "candidate_ranking": False,
                "diagnosis": False,
            },
        )
    )

    observations.append(
        InterpretationObservation(
            code="PHENOTYPE_EXACT_MATCH_SUMMARY",
            label="Exact HPO match summary",
            value={
                "matched_present": matched_present,
                "unmatched_present": unmatched_present,
                "matched_absent": matched_absent,
                "unmatched_absent": unmatched_absent,
            },
            source="phenotype_interpreter",
            details={
                "matching_method": "exact_hpo_identifier",
                "ontology_similarity": False,
                "ancestor_descendant_matching": False,
                "candidate_ranking": False,
            },
        )
    )

    for record in hpo_records:
        evidence = record.get("evidence", {})
        phenotype_ids = evidence.get(
            "phenotype_ids",
            [],
        )

        record_terms = {
            term.strip().upper()
            for term in phenotype_ids
            if isinstance(term, str)
        }

        present_matches = sorted(
            record_terms & present_set
        )

        absent_matches = sorted(
            record_terms & absent_set
        )

        observations.append(
            InterpretationObservation(
                code="PHENOTYPE_SOURCE_ASSOCIATION",
                label="HPO source association",
                value={
                    "source": "HPO",
                    "gene_symbol": evidence.get(
                        "gene_symbol"
                    ),
                    "disease_id": evidence.get(
                        "disease_id"
                    ),
                    "phenotype_ids": phenotype_ids,
                    "evidence_description": (
                        evidence.get(
                            "evidence_description"
                        )
                    ),
                    "present_exact_matches": (
                        present_matches
                    ),
                    "absent_exact_matches": (
                        absent_matches
                    ),
                },
                source="HPO",
                details={
                    "source_record_id": (
                        evidence.get(
                            "source_record_id"
                        )
                    ),
                    "disease_normalization": (
                        record.get(
                            "disease_normalization"
                        )
                    ),
                    "variant_pathogenicity": False,
                    "candidate_ranking": False,
                    "diagnosis": False,
                },
            )
        )

    if not hpo_records:
        status = "UNAVAILABLE"
    elif unmatched_present or matched_absent:
        status = "PARTIAL"
    else:
        status = "AVAILABLE"

    limitations = [
        "Phenotype interpretation v0.1 uses exact HPO identifier matching only.",
        "No ontology ancestor, descendant, semantic similarity, or information-content matching is performed.",
        "Absence of an exact HPO match does not exclude a gene or disease.",
        "An HPO term explicitly absent in the patient is preserved as negative phenotype evidence but does not automatically exclude a gene or disease.",
        "HPO source associations are descriptive evidence and do not establish variant pathogenicity.",
        "Phenotype evidence does not constitute ACMG/AMP variant classification.",
        "Phenotype evidence does not constitute a diagnosis.",
        "No candidate ranking or phenotype similarity score is produced.",
        "Some HPO associations may encode inheritance concepts; these are preserved here and may be interpreted separately by the inheritance layer.",
    ]

    return EvidenceDimension(
        name="phenotype",
        status=status,
        observations=observations,
        limitations=limitations,
    )
