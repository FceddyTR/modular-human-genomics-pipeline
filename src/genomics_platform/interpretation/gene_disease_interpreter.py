from __future__ import annotations

from typing import Any, Dict

from genomics_platform.interpretation.interpretation_profile import (
    EvidenceDimension,
    InterpretationObservation,
)


def interpret_gene_disease(
    gene_evidence_bundle: Dict[str, Any] | None,
) -> EvidenceDimension:
    """
    Interpret an already-built gene evidence bundle descriptively.

    This layer reports gene-disease evidence context only.

    It intentionally does NOT:
      - classify variants,
      - implement ACMG/AMP,
      - assign pathogenicity,
      - assign candidate tiers,
      - rank genes or variants,
      - convert GenCC classifications into variant classifications,
      - convert PanelApp confidence into variant pathogenicity.
    """

    if gene_evidence_bundle is None:
        return EvidenceDimension(
            name="gene_disease",
            status="UNAVAILABLE",
            observations=(),
            limitations=(
                "No gene evidence bundle was supplied.",
                "Absence of gene-disease evidence does not imply that "
                "the gene or variant is clinically irrelevant.",
            ),
        )

    if not isinstance(gene_evidence_bundle, dict):
        raise TypeError(
            "interpret_gene_disease expects a gene evidence bundle dict"
        )

    query = gene_evidence_bundle.get("query") or {}
    summary = gene_evidence_bundle.get("summary") or {}
    records = gene_evidence_bundle.get("records") or []
    normalized_diseases = (
        gene_evidence_bundle.get("normalized_diseases") or []
    )

    gene_symbol = query.get("gene_symbol")

    observations = []

    if gene_symbol:
        observations.append(
            InterpretationObservation(
                code="GENE_DISEASE_GENE",
                label="Gene evaluated for disease associations",
                value=gene_symbol,
                source="GeneEvidenceBundle",
                details={},
            )
        )

    source_counts = summary.get("source_counts") or {}

    if source_counts:
        observations.append(
            InterpretationObservation(
                code="GENE_DISEASE_SOURCE_COUNTS",
                label="Gene-disease evidence records by source",
                value=dict(source_counts),
                source="GeneEvidenceBundle",
                details={
                    "interpretation": (
                        "Record counts describe source availability "
                        "and are not evidence-strength scores."
                    )
                },
            )
        )

    if normalized_diseases:
        observations.append(
            InterpretationObservation(
                code="GENE_DISEASE_NORMALIZED_DISEASES",
                label="MONDO-normalized disease associations",
                value=[
                    {
                        "mondo_id": disease.get("mondo_id"),
                        "mondo_label": disease.get("mondo_label"),
                    }
                    for disease in normalized_diseases
                ],
                source="MONDO",
                details={
                    "normalization_only": True,
                    "variant_classification": False,
                },
            )
        )

    for item in records:
        if not isinstance(item, dict):
            continue

        evidence = item.get("evidence") or {}
        normalization = (
            item.get("disease_normalization") or {}
        )

        source = evidence.get("source", "UNKNOWN")

        observations.append(
            InterpretationObservation(
                code="GENE_DISEASE_SOURCE_RECORD",
                label="Gene-disease source record",
                value={
                    "source": source,
                    "disease_id": evidence.get("disease_id"),
                    "disease_label": (
                        evidence.get("disease_label")
                        or evidence.get("disease_name")
                    ),
                    "normalized_mondo_id": (
                        normalization.get("mondo_id")
                    ),
                    "normalized_mondo_label": (
                        normalization.get("mondo_label")
                    ),
                    "normalization_status": (
                        normalization.get("status")
                    ),
                    "evidence": evidence,
                },
                source=source,
                details={
                    "gene_disease_context_only": True,
                    "variant_pathogenicity": False,
                    "candidate_ranking": False,
                },
            )
        )

    record_count = summary.get(
        "record_count",
        len(records),
    )

    if record_count == 0:
        return EvidenceDimension(
            name="gene_disease",
            status="UNAVAILABLE",
            observations=tuple(observations),
            limitations=(
                "No gene-disease records were present in the supplied "
                "evidence bundle.",
                "No records does not mean the gene has no biological "
                "or clinical relevance.",
                "Source coverage and database scope may be incomplete.",
            ),
        )

    normalization_counts = (
        summary.get("disease_normalization") or {}
    )

    normalization_failures = sum(
        count
        for status, count in normalization_counts.items()
        if status not in {"FOUND"}
    )

    status = (
        "PARTIAL"
        if normalization_failures > 0
        else "AVAILABLE"
    )

    limitations = [
        (
            "Gene-disease associations describe evidence about a gene "
            "and disease relationship; they do not classify the "
            "patient's variant."
        ),
        (
            "GenCC classifications must not be interpreted as ACMG/AMP "
            "variant classifications."
        ),
        (
            "PanelApp confidence levels describe panel/gene evidence "
            "and must not be interpreted as variant pathogenicity."
        ),
        (
            "Orphanet and HPO associations describe disease context "
            "and do not establish pathogenicity of a specific variant."
        ),
        (
            "MONDO mapping normalizes disease identifiers and does not "
            "establish causality."
        ),
    ]

    if normalization_failures:
        limitations.append(
            "One or more disease identifiers could not be normalized "
            "to MONDO; raw source evidence remains preserved."
        )

    return EvidenceDimension(
        name="gene_disease",
        status=status,
        observations=tuple(observations),
        limitations=tuple(limitations),
    )
