"""Mode-aware candidate prioritization engine."""

from __future__ import annotations

from typing import Iterable

from genomics_platform.interpretation.case_model import (
    CandidateCase,
)
from genomics_platform.ranking.ranking_contract import (
    ComponentContribution,
    RankedCandidate,
    RankingResult,
)
from genomics_platform.ranking.ranking_modes import (
    RankingMode,
    normalize_ranking_mode,
)
from genomics_platform.ranking.components.evidence_components import (
    clinical_component,
    population_component,
    functional_component,
    gene_disease_component,
    inheritance_component,
)


ENGINE_VERSION = "0.2.0"

MODE_WEIGHTS = {
    RankingMode.VARIANT_FIRST: {
        "clinical": 25.0,
        "population": 20.0,
        "functional": 15.0,
        "gene_disease": 20.0,
        "inheritance": 20.0,
        "phenotype": 0.0,
    },
    RankingMode.PHENOTYPE: {
        "clinical": 20.0,
        "population": 15.0,
        "functional": 10.0,
        "gene_disease": 20.0,
        "inheritance": 15.0,
        "phenotype": 20.0,
    },
}


def _phenotype_component(
    candidate: CandidateCase,
    maximum: float,
) -> ComponentContribution:
    if maximum == 0:
        return ComponentContribution(
            name="phenotype",
            contribution=0.0,
            max_contribution=0.0,
            status="UNAVAILABLE",
            limitations=(
                "Phenotype evidence is intentionally "
                "not used in variant-first mode.",
            ),
        )

    if candidate.gene_consistency != "MATCH":
        return ComponentContribution(
            name="phenotype",
            contribution=0.0,
            max_contribution=maximum,
            status="BLOCKED",
            limitations=(
                "Phenotype evidence blocked because "
                f"gene consistency is "
                f"{candidate.gene_consistency}.",
            ),
        )

    dimension = (
        candidate.interpretation_profile
        .phenotype
    )

    if dimension is None:
        return ComponentContribution(
            name="phenotype",
            contribution=0.0,
            max_contribution=maximum,
            status="UNAVAILABLE",
            limitations=(
                "Phenotype context was not evaluated.",
            ),
        )

    summaries = tuple(
        observation
        for observation
        in dimension.observations
        if observation.code
        == "PHENOTYPE_EXACT_MATCH_SUMMARY"
    )

    if not summaries:
        return ComponentContribution(
            name="phenotype",
            contribution=0.0,
            max_contribution=maximum,
            status="UNAVAILABLE",
            limitations=(
                "No phenotype exact-match summary "
                "is available.",
            ),
        )

    value = summaries[0].value or {}

    matched = tuple(
        value.get(
            "matched_present"
        ) or ()
    )

    unmatched = tuple(
        value.get(
            "unmatched_present"
        ) or ()
    )

    matched_absent = tuple(
        value.get(
            "matched_absent"
        ) or ()
    )

    total = (
        len(matched)
        + len(unmatched)
    )

    if total == 0:
        return ComponentContribution(
            name="phenotype",
            contribution=0.0,
            max_contribution=maximum,
            status="UNAVAILABLE",
            limitations=(
                "No present HPO terms are available.",
            ),
        )

    fraction = (
        len(matched)
        / total
    )

    contribution = (
        maximum
        * fraction
    )

    limitations = [
        "Phenotype v0.2 uses exact HPO matching only.",
        "Phenotype match does not establish diagnosis "
        "or variant pathogenicity.",
    ]

    if matched_absent:
        limitations.append(
            f"{len(matched_absent)} explicitly absent "
            "phenotype(s) matched; recorded without "
            "hard exclusion."
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
            f"{len(matched)}/{total} present "
            "HPO term(s) matched exactly.",
        ),
        limitations=tuple(limitations),
    )


def _evaluate(
    candidate: CandidateCase,
    mode: RankingMode,
):
    weights = MODE_WEIGHTS[mode]

    warnings = []

    if candidate.gene_consistency == "MISMATCH":
        warnings.append(
            "DATA_INTEGRITY_WARNING: candidate gene "
            "does not match functional annotation "
            "gene(s)."
        )

    elif candidate.gene_consistency == "UNKNOWN":
        warnings.append(
            "GENE_CONTEXT_WARNING: candidate gene "
            "consistency could not be evaluated."
        )

    components = (
        clinical_component(
            candidate,
            weights["clinical"],
        ),
        population_component(
            candidate,
            weights["population"],
        ),
        functional_component(
            candidate,
            weights["functional"],
        ),
        gene_disease_component(
            candidate,
            weights["gene_disease"],
        ),
        inheritance_component(
            candidate,
            weights["inheritance"],
        ),
        _phenotype_component(
            candidate,
            weights["phenotype"],
        ),
    )

    score = round(
        sum(
            component.contribution
            for component in components
        ),
        6,
    )

    return (
        score,
        components,
        tuple(warnings),
    )


def rank_candidates(
    candidates: Iterable[CandidateCase],
    mode: str | RankingMode = (
        RankingMode.VARIANT_FIRST
    ),
) -> RankingResult:
    mode = normalize_ranking_mode(
        mode
    )

    candidates = tuple(candidates)

    if any(
        not isinstance(
            candidate,
            CandidateCase,
        )
        for candidate in candidates
    ):
        raise TypeError(
            "All ranking inputs must be "
            "CandidateCase objects."
        )

    evaluated = []

    for candidate in candidates:
        score, components, warnings = (
            _evaluate(
                candidate,
                mode,
            )
        )

        evaluated.append(
            (
                candidate,
                score,
                components,
                warnings,
            )
        )

    evaluated.sort(
        key=lambda item: (
            -item[1],
            item[0].candidate_id,
        )
    )

    ranked = tuple(
        RankedCandidate(
            candidate_id=(
                candidate.candidate_id
            ),
            variant_key=(
                candidate.variant_key
            ),
            gene_symbol=(
                candidate.gene_symbol
            ),
            rank=rank,
            prioritization_score=score,
            components=components,
            warnings=warnings,
        )
        for rank, (
            candidate,
            score,
            components,
            warnings,
        )
        in enumerate(
            evaluated,
            start=1,
        )
    )

    return RankingResult(
        mode=mode,
        candidates=ranked,
        engine_version=ENGINE_VERSION,
    )
