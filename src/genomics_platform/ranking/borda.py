"""Borda rank aggregation for candidate prioritization.

This module aggregates existing evidence-component contributions
into a candidate ranking.

It does NOT:
- classify variants,
- estimate pathogenicity probability,
- perform ACMG/AMP classification,
- diagnose disease,
- use benchmark truth.

Unavailable, blocked, and error component observations abstain from
that component ballot rather than being interpreted as negative
biological evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

from genomics_platform.ranking.ranking_contract import (
    RankedCandidate,
    RankingResult,
)


ENGINE_VERSION = "0.1.0"
SCHEMA_VERSION = "1.0"

ABSTAIN_STATUSES = frozenset(
    {
        "UNAVAILABLE",
        "BLOCKED",
        "ERROR",
    }
)

ELIGIBLE_STATUSES = frozenset(
    {
        "SUPPORTING",
        "NEUTRAL",
    }
)


@dataclass(frozen=True)
class BordaBallotEntry:
    """One candidate's result within one evidence ballot."""

    candidate_id: str
    component_name: str
    normalized_evidence: float | None
    rank: float | None
    borda_points: float
    status: str
    abstained: bool

    def to_dict(self) -> dict[str, Any]:
        return {
            "candidate_id": self.candidate_id,
            "component_name": self.component_name,
            "normalized_evidence": (
                self.normalized_evidence
            ),
            "rank": self.rank,
            "borda_points": self.borda_points,
            "status": self.status,
            "abstained": self.abstained,
        }


@dataclass(frozen=True)
class BordaCandidate:
    """Aggregated Borda result for one candidate."""

    candidate_id: str
    variant_key: str
    gene_symbol: str | None
    rank: int
    borda_score: float
    ballots_participated: int
    ballots_available: int
    ballot_entries: tuple[
        BordaBallotEntry,
        ...
    ]
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "candidate_id": self.candidate_id,
            "variant_key": self.variant_key,
            "gene_symbol": self.gene_symbol,
            "rank": self.rank,
            "borda_score": self.borda_score,
            "ballots_participated": (
                self.ballots_participated
            ),
            "ballots_available": (
                self.ballots_available
            ),
            "ballot_entries": [
                entry.to_dict()
                for entry in self.ballot_entries
            ],
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class BordaResult:
    """Candidate ranking produced by Borda aggregation."""

    source_mode: str
    candidates: tuple[
        BordaCandidate,
        ...
    ]
    component_names: tuple[str, ...]
    engine_version: str = ENGINE_VERSION
    schema_version: str = SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "engine": {
                "name": "borda_rank_aggregation",
                "version": self.engine_version,
            },
            "strategy": "borda",
            "source_mode": self.source_mode,
            "candidate_count": len(
                self.candidates
            ),
            "component_names": list(
                self.component_names
            ),
            "candidates": [
                candidate.to_dict()
                for candidate in self.candidates
            ],
            "aggregation_semantics": {
                "ballot_unit": (
                    "evidence_component"
                ),
                "component_strength": (
                    "normalized contribution"
                ),
                "tie_handling": (
                    "average rank and average "
                    "Borda points"
                ),
                "missing_evidence": (
                    "UNAVAILABLE, BLOCKED, and "
                    "ERROR abstain"
                ),
                "candidate_score": (
                    "mean Borda points across "
                    "participating ballots, "
                    "adjusted by evidence coverage"
                ),
                "coverage_adjustment": (
                    "participating ballots divided "
                    "by available component ballots"
                ),
            },
            "semantic_boundaries": {
                "relative_prioritization_only": True,
                "pathogenicity_probability": False,
                "acmg_amp_classification": False,
                "variant_classification": False,
                "diagnosis": False,
                "candidate_tier": False,
                "benchmark_truth_used": False,
            },
        }


def _component_map(
    candidate: RankedCandidate,
) -> dict[str, Any]:
    return {
        component.name: component
        for component in candidate.components
    }


def _normalized_component_value(
    component: Any,
) -> float | None:
    """Return comparable [0, 1] evidence strength.

    None means the candidate abstains from this ballot.
    """

    if component.status in ABSTAIN_STATUSES:
        return None

    if component.status not in ELIGIBLE_STATUSES:
        raise ValueError(
            "Unsupported component status: "
            f"{component.status!r}"
        )

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


def _average_rank(
    ordered_values: list[
        tuple[str, float]
    ],
) -> dict[str, float]:
    """Assign average ranks to tied evidence values."""

    result: dict[str, float] = {}

    index = 0

    while index < len(ordered_values):
        value = ordered_values[index][1]

        end = index + 1

        while (
            end < len(ordered_values)
            and ordered_values[end][1] == value
        ):
            end += 1

        first_rank = index + 1
        last_rank = end

        average = (
            first_rank + last_rank
        ) / 2.0

        for candidate_id, _ in (
            ordered_values[index:end]
        ):
            result[candidate_id] = average

        index = end

    return result


def _ballot(
    candidates: tuple[
        RankedCandidate,
        ...
    ],
    component_name: str,
) -> dict[str, BordaBallotEntry]:
    """Create one component ballot.

    Only candidates with evaluable component evidence
    participate.

    For m participating candidates, ordinary untied
    Borda points are:

        rank 1 -> m - 1
        rank 2 -> m - 2
        ...
        rank m -> 0

    Tied candidates receive the average points of the
    positions occupied by the tie.
    """

    candidate_components = {
        candidate.candidate_id: (
            _component_map(candidate).get(
                component_name
            )
        )
        for candidate in candidates
    }

    eligible: list[
        tuple[str, float]
    ] = []

    for candidate_id, component in (
        candidate_components.items()
    ):
        if component is None:
            continue

        value = _normalized_component_value(
            component
        )

        if value is not None:
            eligible.append(
                (
                    candidate_id,
                    value,
                )
            )

    eligible.sort(
        key=lambda item: (
            -item[1],
            item[0],
        )
    )

    ranks = _average_rank(
        eligible
    )

    participant_count = len(
        eligible
    )

    entries: dict[
        str,
        BordaBallotEntry,
    ] = {}

    for candidate in candidates:
        candidate_id = (
            candidate.candidate_id
        )

        component = (
            candidate_components[
                candidate_id
            ]
        )

        if component is None:
            entries[candidate_id] = (
                BordaBallotEntry(
                    candidate_id=candidate_id,
                    component_name=component_name,
                    normalized_evidence=None,
                    rank=None,
                    borda_points=0.0,
                    status="UNAVAILABLE",
                    abstained=True,
                )
            )
            continue

        normalized = (
            _normalized_component_value(
                component
            )
        )

        if normalized is None:
            entries[candidate_id] = (
                BordaBallotEntry(
                    candidate_id=candidate_id,
                    component_name=component_name,
                    normalized_evidence=None,
                    rank=None,
                    borda_points=0.0,
                    status=component.status,
                    abstained=True,
                )
            )
            continue

        rank = ranks[candidate_id]

        points = max(
            0.0,
            participant_count - rank,
        )

        entries[candidate_id] = (
            BordaBallotEntry(
                candidate_id=candidate_id,
                component_name=component_name,
                normalized_evidence=normalized,
                rank=rank,
                borda_points=points,
                status=component.status,
                abstained=False,
            )
        )

    return entries


def aggregate_borda(
    result: RankingResult,
) -> BordaResult:
    """Aggregate a completed component-based ranking.

    Benchmark truth is neither accepted nor accessible here.
    """

    candidates = tuple(
        result.candidates
    )

    if not candidates:
        return BordaResult(
            source_mode=result.mode.value,
            candidates=(),
            component_names=(),
        )

    component_names = tuple(
        sorted(
            {
                component.name
                for candidate in candidates
                for component
                in candidate.components
                if component.max_contribution > 0
            }
        )
    )

    ballots = {
        component_name: _ballot(
            candidates,
            component_name,
        )
        for component_name
        in component_names
    }

    evaluated = []

    for candidate in candidates:
        entries = tuple(
            ballots[component_name][
                candidate.candidate_id
            ]
            for component_name
            in component_names
        )

        participated = sum(
            not entry.abstained
            for entry in entries
        )

        points = sum(
            entry.borda_points
            for entry in entries
            if not entry.abstained
        )

        # First calculate mean Borda performance
        # across ballots in which the candidate had
        # evaluable evidence.
        #
        # Missing evidence is NOT interpreted as
        # negative biological evidence.
        mean_borda = (
            points / participated
            if participated
            else 0.0
        )

        # Apply a separate evidence-completeness
        # adjustment.
        #
        # Without this adjustment, a candidate with
        # only one available component could win that
        # single ballot and receive an artificially
        # inflated mean score relative to candidates
        # evaluated across the full evidence model.
        #
        # Coverage therefore modifies confidence in
        # the aggregate ranking signal; it does not
        # convert UNAVAILABLE/BLOCKED/ERROR into
        # negative biological evidence.
        coverage = (
            participated / len(component_names)
            if component_names
            else 0.0
        )

        score = (
            mean_borda * coverage
        )

        warnings = list(
            candidate.warnings
        )

        if participated < len(
            component_names
        ):
            warnings.append(
                "BORDA_PARTIAL_EVIDENCE: "
                f"{participated}/"
                f"{len(component_names)} "
                "component ballots participated."
            )

        evaluated.append(
            (
                candidate,
                score,
                participated,
                entries,
                tuple(warnings),
            )
        )

    evaluated.sort(
        key=lambda item: (
            -item[1],
            item[0].candidate_id,
        )
    )

    ranked_candidates = []

    previous_score = None
    previous_rank = None

    for position, (
        candidate,
        score,
        participated,
        entries,
        warnings,
    ) in enumerate(
        evaluated,
        start=1,
    ):
        rounded_score = round(
            score,
            6,
        )

        # Competition ranking:
        #
        #   1, 1, 3
        #
        # Candidate ID remains only a deterministic
        # presentation-order tie-break. It must not
        # create a biological ranking distinction
        # between equal aggregate Borda scores.
        if (
            previous_score is not None
            and rounded_score
            == previous_score
        ):
            rank = previous_rank
        else:
            rank = position

        ranked_candidates.append(
            BordaCandidate(
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
                borda_score=rounded_score,
                ballots_participated=participated,
                ballots_available=len(
                    component_names
                ),
                ballot_entries=entries,
                warnings=warnings,
            )
        )

        previous_score = rounded_score
        previous_rank = rank

    ranked = tuple(
        ranked_candidates
    )

    return BordaResult(
        source_mode=result.mode.value,
        candidates=ranked,
        component_names=component_names,
    )
