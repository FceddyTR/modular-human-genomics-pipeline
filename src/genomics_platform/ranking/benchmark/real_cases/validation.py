"""Structural and truth-separation guards for real-case benchmarks."""

from __future__ import annotations

from dataclasses import fields

from genomics_platform.ranking.benchmark.truth_set import (
    RankingTruth,
)

from .contract import RealCaseInput


FORBIDDEN_RANKING_FIELDS = frozenset(
    {
        "truth",
        "causal_candidate_id",
        "causal_candidate_ids",
        "is_causal",
        "expected_rank",
        "expected_ranking",
        "benchmark_label",
    }
)


def validate_real_case_input(
    case: RealCaseInput,
) -> None:
    """Validate invariants of a truth-blind ranking input."""

    field_names = {
        field.name
        for field in fields(case)
    }

    leaked = (
        field_names
        & FORBIDDEN_RANKING_FIELDS
    )

    if leaked:
        raise ValueError(
            "Benchmark truth fields are forbidden "
            "from real-case ranking input: "
            f"{sorted(leaked)!r}"
        )

    candidate_ids = [
        candidate.candidate_id
        for candidate in case.candidates
    ]

    if len(candidate_ids) != len(
        set(candidate_ids)
    ):
        raise ValueError(
            "Candidate IDs must be unique."
        )


def validate_case_truth_pair(
    case: RealCaseInput,
    truth: RankingTruth,
) -> None:
    """Validate truth only at the post-ranking benchmark boundary."""

    validate_real_case_input(case)

    if case.case_id != truth.case_id:
        raise ValueError(
            "Real-case ID does not match "
            "benchmark truth case ID."
        )

    candidate_ids = {
        candidate.candidate_id
        for candidate in case.candidates
    }

    missing = (
        set(truth.causal_candidate_ids)
        - candidate_ids
    )

    if missing:
        raise ValueError(
            "Benchmark truth references candidates "
            "absent from the ranked cohort: "
            f"{sorted(missing)!r}"
        )
