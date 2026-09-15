"""Ontology-aware phenotype contribution for candidate prioritization.

This component uses normalized Lin Best Match Average (Lin-BMA)
between patient HPO terms and candidate-associated HPO terms.

It provides phenotype relevance for candidate prioritization only.

It does NOT provide:
- diagnosis
- variant pathogenicity
- ACMG/AMP classification
- disease probability
"""

from __future__ import annotations

from genomics_platform.interpretation.case_model import (
    CandidateCase,
)
from genomics_platform.phenotype.hpo_semantic_similarity import (
    HPOSemanticSimilarity,
)
from genomics_platform.ranking.ranking_contract import (
    ComponentContribution,
)


def _candidate_hpo_terms(
    candidate: CandidateCase,
) -> tuple[str, ...]:
    """Extract HPO terms associated with the candidate gene.

    Terms are taken only from phenotype source-association
    observations already present in the interpretation profile.
    """

    dimension = (
        candidate.interpretation_profile
        .phenotype
    )

    if dimension is None:
        return ()

    terms: list[str] = []

    for observation in dimension.observations:
        if (
            observation.code
            != "PHENOTYPE_SOURCE_ASSOCIATION"
        ):
            continue

        value = observation.value or {}

        evidence = (
            value.get("evidence")
            if isinstance(value, dict)
            else None
        )

        if not isinstance(evidence, dict):
            evidence = value

        hpo_id = evidence.get("hpo_id")

        if (
            isinstance(hpo_id, str)
            and hpo_id.startswith("HP:")
        ):
            terms.append(
                hpo_id.strip().upper()
            )

        phenotype_ids = evidence.get(
            "phenotype_ids"
        )

        if isinstance(
            phenotype_ids,
            (list, tuple),
        ):
            for term in phenotype_ids:
                if (
                    isinstance(term, str)
                    and term.startswith("HP:")
                ):
                    terms.append(
                        term.strip().upper()
                    )

    return tuple(
        dict.fromkeys(terms)
    )


def semantic_phenotype_component(
    candidate: CandidateCase,
    maximum: float,
    semantic_engine: HPOSemanticSimilarity,
) -> ComponentContribution:
    """Calculate ontology-aware phenotype contribution."""

    if maximum == 0:
        return ComponentContribution(
            name="phenotype",
            contribution=0.0,
            max_contribution=0.0,
            status="UNAVAILABLE",
            limitations=(
                "Phenotype evidence is intentionally "
                "disabled for this ranking mode.",
            ),
        )

    if candidate.gene_consistency != "MATCH":
        return ComponentContribution(
            name="phenotype",
            contribution=0.0,
            max_contribution=maximum,
            status="BLOCKED",
            limitations=(
                "Semantic phenotype evidence blocked "
                "because gene consistency is "
                f"{candidate.gene_consistency}.",
            ),
        )

    patient_terms = tuple(
        candidate.case_context
        .present_hpo_terms
    )

    if not patient_terms:
        return ComponentContribution(
            name="phenotype",
            contribution=0.0,
            max_contribution=maximum,
            status="UNAVAILABLE",
            limitations=(
                "No present patient HPO terms "
                "are available.",
            ),
        )

    candidate_terms = (
        _candidate_hpo_terms(
            candidate
        )
    )

    if not candidate_terms:
        return ComponentContribution(
            name="phenotype",
            contribution=0.0,
            max_contribution=maximum,
            status="UNAVAILABLE",
            limitations=(
                "No candidate-associated phenotype "
                "terms are available.",
            ),
        )

    similarity = (
        semantic_engine.lin_bma(
            patient_terms,
            candidate_terms,
        )
    )

    fraction = min(
        1.0,
        max(
            0.0,
            float(
                similarity.bma_similarity
            ),
        ),
    )

    contribution = (
        maximum
        * fraction
    )

    return ComponentContribution(
        name="phenotype",
        contribution=round(
            contribution,
            6,
        ),
        max_contribution=maximum,
        status=(
            "SUPPORTING"
            if contribution > 0
            else "NEUTRAL"
        ),
        reasons=(
            "Ontology-aware phenotype similarity "
            f"Lin-BMA={fraction:.6f} across "
            f"{len(patient_terms)} patient and "
            f"{len(candidate_terms)} candidate "
            "HPO term(s).",
        ),
        limitations=(
            "Semantic phenotype similarity is "
            "candidate prioritization context only.",
            "Lin-BMA is not disease probability.",
            "Phenotype similarity does not establish "
            "diagnosis or variant pathogenicity.",
            "Explicitly absent phenotypes are not "
            "used as a hard penalty in v0.1.",
        ),
    )
