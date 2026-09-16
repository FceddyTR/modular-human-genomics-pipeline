"""Controlled synthetic benchmark suite for candidate prioritization.

The suite is designed to exercise ranking behavior under controlled
evidence perturbations.

Important semantic boundary:
benchmark truth is stored separately from candidate construction and
is introduced only after ranking has completed.

Synthetic benchmark performance is not clinical validation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Sequence


SUITE_VERSION = "0.1.0"


@dataclass(frozen=True)
class ControlledCaseSpec:
    """Metadata describing one controlled ranking experiment."""

    case_id: str
    family: str
    description: str
    causal_candidate_ids: tuple[str, ...]
    expected_behavior: str
    tags: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.case_id.strip():
            raise ValueError("case_id is required.")

        if not self.family.strip():
            raise ValueError("family is required.")

        if not self.causal_candidate_ids:
            raise ValueError(
                "At least one causal candidate ID is required."
            )

        if len(set(self.causal_candidate_ids)) != len(
            self.causal_candidate_ids
        ):
            raise ValueError(
                "causal_candidate_ids must be unique."
            )

    def to_dict(self) -> dict[str, Any]:
        return {
            "case_id": self.case_id,
            "family": self.family,
            "description": self.description,
            "causal_candidate_ids": list(
                self.causal_candidate_ids
            ),
            "expected_behavior": self.expected_behavior,
            "tags": list(self.tags),
            "suite_version": SUITE_VERSION,
            "semantic_boundaries": {
                "synthetic_controlled_case": True,
                "clinical_validation": False,
                "truth_used_for_ranking": False,
                "truth_used_for_benchmark": True,
                "pathogenicity_probability": False,
                "diagnosis": False,
            },
        }


@dataclass(frozen=True)
class ControlledCase:
    """Candidate cohort plus separately held benchmark truth."""

    spec: ControlledCaseSpec
    candidates: tuple[Any, ...]

    def __post_init__(self) -> None:
        if not self.candidates:
            raise ValueError(
                "Controlled case must contain candidates."
            )

        candidate_ids = tuple(
            candidate.candidate_id
            for candidate in self.candidates
        )

        if len(set(candidate_ids)) != len(candidate_ids):
            raise ValueError(
                "Candidate IDs must be unique within a case."
            )

        missing_truth = (
            set(self.spec.causal_candidate_ids)
            - set(candidate_ids)
        )

        if missing_truth:
            raise ValueError(
                "Causal candidate IDs are absent from cohort: "
                f"{sorted(missing_truth)!r}"
            )

    @property
    def case_id(self) -> str:
        return self.spec.case_id

    @property
    def family(self) -> str:
        return self.spec.family

    def ranking_candidates(self) -> tuple[Any, ...]:
        """Return candidates without exposing benchmark truth."""
        return self.candidates

    def causal_candidate_ids(self) -> tuple[str, ...]:
        """Expose truth only to the benchmark layer."""
        return self.spec.causal_candidate_ids


class ControlledSuiteRegistry:
    """Registry of independently constructed controlled cases."""

    def __init__(self) -> None:
        self._builders: dict[
            str,
            Callable[[], ControlledCase],
        ] = {}

    def register(
        self,
        case_id: str,
        builder: Callable[[], ControlledCase],
    ) -> None:
        normalized = case_id.strip()

        if not normalized:
            raise ValueError("case_id is required.")

        if normalized in self._builders:
            raise ValueError(
                f"Duplicate controlled case: {normalized}"
            )

        self._builders[normalized] = builder

    def case_ids(self) -> tuple[str, ...]:
        return tuple(sorted(self._builders))

    def build(
        self,
        case_id: str,
    ) -> ControlledCase:
        try:
            builder = self._builders[case_id]
        except KeyError as exc:
            raise KeyError(
                f"Unknown controlled case: {case_id}"
            ) from exc

        case = builder()

        if case.case_id != case_id:
            raise ValueError(
                "Registered case ID does not match "
                "constructed case ID."
            )

        return case

    def build_all(self) -> tuple[ControlledCase, ...]:
        return tuple(
            self.build(case_id)
            for case_id in self.case_ids()
        )

    def family_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}

        for case in self.build_all():
            counts[case.family] = (
                counts.get(case.family, 0) + 1
            )

        return dict(sorted(counts.items()))


def validate_suite(
    cases: Sequence[ControlledCase],
) -> None:
    """Validate suite-level structural invariants."""

    if not cases:
        raise ValueError(
            "Controlled benchmark suite is empty."
        )

    case_ids = [
        case.case_id
        for case in cases
    ]

    if len(case_ids) != len(set(case_ids)):
        raise ValueError(
            "Controlled benchmark case IDs must be unique."
        )

    for case in cases:
        candidate_ids = {
            candidate.candidate_id
            for candidate in case.candidates
        }

        causal_ids = set(
            case.causal_candidate_ids()
        )

        if not causal_ids <= candidate_ids:
            raise ValueError(
                f"Truth mismatch in case {case.case_id!r}."
            )
