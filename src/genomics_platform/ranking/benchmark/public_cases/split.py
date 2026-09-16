"""Leakage-aware deterministic public benchmark splitting."""

from __future__ import annotations

from dataclasses import dataclass

from .cohort import (
    CohortSelection,
    PublicCaseRecord,
)


@dataclass(frozen=True)
class PublicBenchmarkSplit:
    development: tuple[
        PublicCaseRecord, ...
    ]
    locked_evaluation: tuple[
        PublicCaseRecord, ...
    ]
    split_version: str = "0.1.0"

    def __post_init__(self) -> None:
        all_cases = (
            self.development
            + self.locked_evaluation
        )

        case_ids = [
            case.case_id
            for case in all_cases
        ]

        if len(case_ids) != len(
            set(case_ids)
        ):
            raise ValueError(
                "Case overlap across benchmark "
                "splits."
            )

        development_publications = {
            case.publication_id
            for case in self.development
        }

        locked_publications = {
            case.publication_id
            for case in self.locked_evaluation
        }

        if (
            development_publications
            & locked_publications
        ):
            raise ValueError(
                "Publication leakage across "
                "benchmark splits."
            )

        development_variants = {
            key
            for case in self.development
            for key in (
                case.causal_variant_keys
            )
        }

        locked_variants = {
            key
            for case
            in self.locked_evaluation
            for key in (
                case.causal_variant_keys
            )
        }

        if (
            development_variants
            & locked_variants
        ):
            raise ValueError(
                "Causal variant leakage across "
                "benchmark splits."
            )

    @property
    def development_genes(
        self,
    ) -> frozenset[str]:
        return frozenset(
            gene
            for case in self.development
            for gene in case.causal_genes
        )

    @property
    def locked_genes(
        self,
    ) -> frozenset[str]:
        return frozenset(
            gene
            for case
            in self.locked_evaluation
            for gene in case.causal_genes
        )

    @property
    def locked_seen_genes(
        self,
    ) -> frozenset[str]:
        return (
            self.locked_genes
            & self.development_genes
        )

    @property
    def locked_unseen_genes(
        self,
    ) -> frozenset[str]:
        return (
            self.locked_genes
            - self.development_genes
        )


def split_public_benchmark(
    cohort: CohortSelection,
    *,
    development_size: int = 300,
    locked_size: int = 200,
) -> PublicBenchmarkSplit:
    """Split by publication groups without leakage."""

    if development_size <= 0:
        raise ValueError(
            "development_size must be positive."
        )

    if locked_size <= 0:
        raise ValueError(
            "locked_size must be positive."
        )

    if (
        development_size
        + locked_size
        != len(cohort.cases)
    ):
        raise ValueError(
            "Requested split sizes must equal "
            "cohort size."
        )

    publication_groups: dict[
        str,
        list[PublicCaseRecord],
    ] = {}

    for case in cohort.cases:
        publication_groups.setdefault(
            case.publication_id,
            [],
        ).append(case)

    groups = [
        tuple(group)
        for _, group in sorted(
            publication_groups.items()
        )
    ]

    # Large groups are placed first so exact
    # capacities are easier to satisfy.
    groups.sort(
        key=lambda group: (
            -len(group),
            group[0].publication_id,
        )
    )

    development: list[
        PublicCaseRecord
    ] = []

    locked: list[
        PublicCaseRecord
    ] = []

    for group in groups:
        dev_remaining = (
            development_size
            - len(development)
        )

        locked_remaining = (
            locked_size
            - len(locked)
        )

        size = len(group)

        can_dev = size <= dev_remaining
        can_locked = (
            size <= locked_remaining
        )

        if not can_dev and not can_locked:
            raise ValueError(
                "Publication-disjoint split "
                "cannot satisfy requested exact "
                "sizes with this cohort."
            )

        if can_dev and not can_locked:
            development.extend(group)
            continue

        if can_locked and not can_dev:
            locked.extend(group)
            continue

        # Prefer the side with the larger
        # proportional remaining capacity.
        dev_fraction = (
            dev_remaining
            / development_size
        )

        locked_fraction = (
            locked_remaining
            / locked_size
        )

        if dev_fraction >= locked_fraction:
            development.extend(group)
        else:
            locked.extend(group)

    if (
        len(development)
        != development_size
        or len(locked)
        != locked_size
    ):
        raise ValueError(
            "Publication-disjoint split did not "
            "reach requested exact sizes."
        )

    return PublicBenchmarkSplit(
        development=tuple(development),
        locked_evaluation=tuple(locked),
    )
