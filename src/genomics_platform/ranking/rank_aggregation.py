"""Experimental rank-based candidate aggregation.

This module provides alternative aggregation strategies for comparing
candidate evidence profiles.

The geometric-rank strategy is inspired by the general rank-aggregation
principle used by MINI-EX: candidates are ranked independently within
each metric/component and the resulting ranks are combined using their
geometric mean.

Important semantic differences for this genomics implementation:

- UNAVAILABLE is not negative biological evidence.
- BLOCKED is not negative biological evidence.
- ERROR is a technical/evidence interpretation failure.
- Missing component values therefore abstain rather than being assigned
  an artificial worst biological rank.
- Benchmark truth is never consumed by this module.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from genomics_platform.ranking.ranking_contract import (
    RankingResult,
)


class RankAggregationMethod(str, Enum):
    GEOMETRIC_MEAN = "geometric_mean"
    ARITHMETIC_MEAN = "arithmetic_mean"


ELIGIBLE_STATUSES = frozenset(
    {
        "SUPPORTING",
        "NEUTRAL",
    }
)

ABSTAIN_STATUSES = frozenset(
    {
        "UNAVAILABLE",
        "BLOCKED",
        "ERROR",
    }
)


@dataclass(frozen=True)
class RankComponentEntry:
    candidate_id: str
    component_name: str
    normalized_evidence: float | None
    component_rank: float | None
    status: str
    abstained: bool


@dataclass(frozen=True)
class RankAggregatedCandidate:
    candidate_id: str
    variant_key: str
    gene_symbol: str | None
    rank: int
    aggregate_rank_score: float
    components_participated: int
    components_available: int
    evidence_coverage: float
    component_entries: tuple[RankComponentEntry, ...]
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class RankAggregationResult:
    source_mode: str
    method: RankAggregationMethod
    candidates: tuple[RankAggregatedCandidate, ...]
    component_names: tuple[str, ...]
    engine_version: str = "0.1.0"
    schema_version: str = "1.0"

    def to_dict(self) -> dict:
        return {
            "engine": "experimental_rank_aggregation",
            "engine_version": self.engine_version,
            "schema_version": self.schema_version,
            "source_mode": self.source_mode,
            "method": self.method.value,
            "component_names": list(
                self.component_names
            ),
            "semantic_boundaries": {
                "relative_prioritization_only": True,
                "pathogenicity_probability": False,
                "acmg_classification": False,
                "diagnosis": False,
                "clinical_validation": False,
                "benchmark_truth_used": False,
            },
            "candidates": [
                {
                    "candidate_id": c.candidate_id,
                    "variant_key": c.variant_key,
                    "gene_symbol": c.gene_symbol,
                    "rank": c.rank,
                    "aggregate_rank_score": (
                        c.aggregate_rank_score
                    ),
                    "components_participated": (
                        c.components_participated
                    ),
                    "components_available": (
                        c.components_available
                    ),
                    "evidence_coverage": (
                        c.evidence_coverage
                    ),
                    "warnings": list(c.warnings),
                    "component_entries": [
                        {
                            "candidate_id": (
                                e.candidate_id
                            ),
                            "component_name": (
                                e.component_name
                            ),
                            "normalized_evidence": (
                                e.normalized_evidence
                            ),
                            "component_rank": (
                                e.component_rank
                            ),
                            "status": e.status,
                            "abstained": e.abstained,
                        }
                        for e in c.component_entries
                    ],
                }
                for c in self.candidates
            ],
        }


def _component_names(
    result: RankingResult,
) -> tuple[str, ...]:
    names = []

    for candidate in result.candidates:
        for component in candidate.components:
            if component.name not in names:
                names.append(component.name)

    return tuple(names)


def _normalize_component(component) -> float | None:
    if component.status not in ELIGIBLE_STATUSES:
        return None

    if component.max_contribution <= 0:
        return None

    value = (
        component.contribution
        / component.max_contribution
    )

    return max(
        0.0,
        min(1.0, float(value)),
    )


def _average_competition_ranks(
    values: list[tuple[str, float]],
) -> dict[str, float]:
    """Assign average ranks, with larger evidence being better."""

    if not values:
        return {}

    ordered = sorted(
        values,
        key=lambda item: (
            -item[1],
            item[0],
        ),
    )

    ranks: dict[str, float] = {}

    index = 0
    while index < len(ordered):
        value = ordered[index][1]
        end = index + 1

        while (
            end < len(ordered)
            and math.isclose(
                ordered[end][1],
                value,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
        ):
            end += 1

        first_rank = index + 1
        last_rank = end

        average_rank = (
            first_rank + last_rank
        ) / 2.0

        for position in range(index, end):
            candidate_id = ordered[position][0]
            ranks[candidate_id] = average_rank

        index = end

    return ranks


def _aggregate_ranks(
    ranks: list[float],
    method: RankAggregationMethod,
) -> float:
    if not ranks:
        return math.inf

    if method == RankAggregationMethod.GEOMETRIC_MEAN:
        return math.exp(
            sum(math.log(rank) for rank in ranks)
            / len(ranks)
        )

    if method == RankAggregationMethod.ARITHMETIC_MEAN:
        return sum(ranks) / len(ranks)

    raise ValueError(
        f"Unsupported rank aggregation method: {method}"
    )


def aggregate_component_ranks(
    result: RankingResult,
    *,
    method: RankAggregationMethod = (
        RankAggregationMethod.GEOMETRIC_MEAN
    ),
) -> RankAggregationResult:
    """Aggregate candidate component ranks.

    Smaller aggregate rank score is better.

    Missing/unavailable/error/blocked evidence abstains and is reflected
    through evidence coverage. It is not converted into negative
    biological evidence.
    """

    component_names = _component_names(result)

    candidate_by_id = {
        candidate.candidate_id: candidate
        for candidate in result.candidates
    }

    normalized: dict[
        tuple[str, str],
        float | None,
    ] = {}

    statuses: dict[
        tuple[str, str],
        str,
    ] = {}

    for candidate in result.candidates:
        components = {
            component.name: component
            for component in candidate.components
        }

        for name in component_names:
            component = components.get(name)

            if component is None:
                normalized[
                    (candidate.candidate_id, name)
                ] = None
                statuses[
                    (candidate.candidate_id, name)
                ] = "UNAVAILABLE"
                continue

            normalized[
                (candidate.candidate_id, name)
            ] = _normalize_component(component)

            statuses[
                (candidate.candidate_id, name)
            ] = component.status

    component_ranks: dict[
        tuple[str, str],
        float | None,
    ] = {}

    for name in component_names:
        values = []

        for candidate in result.candidates:
            candidate_id = candidate.candidate_id

            value = normalized[
                (candidate_id, name)
            ]

            if value is not None:
                values.append(
                    (candidate_id, value)
                )

        ranks = _average_competition_ranks(
            values
        )

        for candidate in result.candidates:
            candidate_id = candidate.candidate_id

            component_ranks[
                (candidate_id, name)
            ] = ranks.get(candidate_id)

    provisional = []

    for candidate in result.candidates:
        entries = []
        ranks = []

        for name in component_names:
            key = (
                candidate.candidate_id,
                name,
            )

            value = normalized[key]
            component_rank = (
                component_ranks[key]
            )
            status = statuses[key]

            abstained = (
                component_rank is None
            )

            if component_rank is not None:
                ranks.append(component_rank)

            entries.append(
                RankComponentEntry(
                    candidate_id=(
                        candidate.candidate_id
                    ),
                    component_name=name,
                    normalized_evidence=value,
                    component_rank=component_rank,
                    status=status,
                    abstained=abstained,
                )
            )

        participated = len(ranks)
        available = len(component_names)

        coverage = (
            participated / available
            if available
            else 0.0
        )

        score = _aggregate_ranks(
            ranks,
            method,
        )

        warnings = []

        if participated == 0:
            warnings.append(
                "No eligible evidence components "
                "participated in rank aggregation."
            )

        if coverage < 1.0:
            warnings.append(
                "Some evidence components abstained; "
                "aggregate rank score is based on "
                "partial evidence coverage."
            )

        provisional.append(
            {
                "candidate": candidate,
                "score": score,
                "participated": participated,
                "available": available,
                "coverage": coverage,
                "entries": tuple(entries),
                "warnings": tuple(warnings),
            }
        )

    ordered = sorted(
        provisional,
        key=lambda item: (
            item["score"],
            -item["coverage"],
            item["candidate"].candidate_id,
        ),
    )

    output = []

    previous_score = None
    previous_rank = None

    for index, item in enumerate(
        ordered,
        start=1,
    ):
        score = item["score"]

        if (
            previous_score is not None
            and (
                (
                    math.isinf(score)
                    and math.isinf(
                        previous_score
                    )
                )
                or math.isclose(
                    score,
                    previous_score,
                    rel_tol=0.0,
                    abs_tol=1e-12,
                )
            )
        ):
            rank = previous_rank
        else:
            rank = index

        candidate = item["candidate"]

        output.append(
            RankAggregatedCandidate(
                candidate_id=(
                    candidate.candidate_id
                ),
                variant_key=candidate.variant_key,
                gene_symbol=candidate.gene_symbol,
                rank=rank,
                aggregate_rank_score=score,
                components_participated=(
                    item["participated"]
                ),
                components_available=(
                    item["available"]
                ),
                evidence_coverage=(
                    item["coverage"]
                ),
                component_entries=(
                    item["entries"]
                ),
                warnings=item["warnings"],
            )
        )

        previous_score = score
        previous_rank = rank

    return RankAggregationResult(
        source_mode=result.mode.value,
        method=method,
        candidates=tuple(output),
        component_names=component_names,
    )
