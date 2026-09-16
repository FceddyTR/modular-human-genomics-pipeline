"""Deterministic diversity-first public benchmark cohort construction."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from .phenopacket import PublicCaseEligibility


@dataclass(frozen=True)
class PublicCaseRecord:
    """Eligible public case plus provenance needed for cohort selection."""

    eligibility: PublicCaseEligibility
    publication_id: str
    source_reference: str

    def __post_init__(self) -> None:
        if not self.eligibility.eligible:
            raise ValueError(
                "PublicCaseRecord requires an "
                "eligible case."
            )

        publication_id = (
            self.publication_id.strip()
        )
        source_reference = (
            self.source_reference.strip()
        )

        if not publication_id:
            raise ValueError(
                "publication_id is required."
            )

        if not source_reference:
            raise ValueError(
                "source_reference is required."
            )

        object.__setattr__(
            self,
            "publication_id",
            publication_id,
        )
        object.__setattr__(
            self,
            "source_reference",
            source_reference,
        )

    @property
    def case_id(self) -> str:
        return self.eligibility.case_id

    @property
    def causal_genes(self) -> tuple[str, ...]:
        return self.eligibility.causal_genes

    @property
    def causal_variant_keys(
        self,
    ) -> tuple[str, ...]:
        return tuple(
            sorted(
                variant.variant_key
                for variant
                in self.eligibility
                .causal_variants
            )
        )

    @property
    def causal_variant_set_key(
        self,
    ) -> tuple[str, ...]:
        return self.causal_variant_keys


@dataclass(frozen=True)
class CohortSelection:
    cases: tuple[PublicCaseRecord, ...]
    target_size: int
    selector_version: str = "0.1.0"

    def __post_init__(self) -> None:
        if self.target_size <= 0:
            raise ValueError(
                "target_size must be positive."
            )

        if len(self.cases) != self.target_size:
            raise ValueError(
                "Cohort size does not match "
                "target_size."
            )

        case_ids = [
            case.case_id
            for case in self.cases
        ]

        if len(case_ids) != len(
            set(case_ids)
        ):
            raise ValueError(
                "Duplicate case_id in cohort."
            )

        variant_sets = [
            case.causal_variant_set_key
            for case in self.cases
        ]

        if len(variant_sets) != len(
            set(variant_sets)
        ):
            raise ValueError(
                "Duplicate exact causal variant "
                "set in cohort."
            )

    @property
    def unique_gene_count(self) -> int:
        return len(
            {
                gene
                for case in self.cases
                for gene in case.causal_genes
            }
        )

    @property
    def unique_publication_count(
        self,
    ) -> int:
        return len(
            {
                case.publication_id
                for case in self.cases
            }
        )


def _case_sort_key(
    case: PublicCaseRecord,
) -> tuple:
    return (
        case.publication_id,
        case.causal_genes,
        case.causal_variant_set_key,
        case.case_id,
        case.source_reference,
    )


def _deduplicate_records(
    records: Iterable[PublicCaseRecord],
) -> tuple[PublicCaseRecord, ...]:
    by_case_id: dict[
        str,
        PublicCaseRecord,
    ] = {}

    for record in records:
        existing = by_case_id.get(
            record.case_id
        )

        if (
            existing is not None
            and existing != record
        ):
            raise ValueError(
                "Conflicting records for "
                f"case_id {record.case_id!r}."
            )

        by_case_id[
            record.case_id
        ] = record

    ordered = sorted(
        by_case_id.values(),
        key=_case_sort_key,
    )

    seen_variant_sets: set[
        tuple[str, ...]
    ] = set()

    unique: list[
        PublicCaseRecord
    ] = []

    for record in ordered:
        key = (
            record.causal_variant_set_key
        )

        if key in seen_variant_sets:
            continue

        seen_variant_sets.add(key)
        unique.append(record)

    return tuple(unique)


def select_diverse_cohort(
    records: Iterable[PublicCaseRecord],
    *,
    target_size: int,
) -> CohortSelection:
    """Select a deterministic diversity-first cohort.

    Greedy priorities:
    1. introduce a previously unseen causal gene
    2. introduce a previously unseen publication
    3. prefer lower current publication representation
    4. prefer lower current gene representation
    5. deterministic lexical identity
    """

    if target_size <= 0:
        raise ValueError(
            "target_size must be positive."
        )

    remaining = list(
        _deduplicate_records(records)
    )

    if len(remaining) < target_size:
        raise ValueError(
            "Insufficient unique eligible cases "
            "for requested cohort size."
        )

    selected: list[
        PublicCaseRecord
    ] = []

    seen_genes: set[str] = set()
    seen_publications: set[str] = set()

    gene_counts: dict[str, int] = {}
    publication_counts: dict[str, int] = {}

    while len(selected) < target_size:
        def priority(
            case: PublicCaseRecord,
        ) -> tuple:
            new_gene_count = sum(
                gene not in seen_genes
                for gene in case.causal_genes
            )

            new_publication = (
                case.publication_id
                not in seen_publications
            )

            publication_count = (
                publication_counts.get(
                    case.publication_id,
                    0,
                )
            )

            gene_burden = sum(
                gene_counts.get(
                    gene,
                    0,
                )
                for gene
                in case.causal_genes
            )

            return (
                -new_gene_count,
                -int(new_publication),
                publication_count,
                gene_burden,
                _case_sort_key(case),
            )

        chosen = min(
            remaining,
            key=priority,
        )

        remaining.remove(chosen)
        selected.append(chosen)

        seen_publications.add(
            chosen.publication_id
        )

        publication_counts[
            chosen.publication_id
        ] = (
            publication_counts.get(
                chosen.publication_id,
                0,
            )
            + 1
        )

        for gene in chosen.causal_genes:
            seen_genes.add(gene)

            gene_counts[gene] = (
                gene_counts.get(
                    gene,
                    0,
                )
                + 1
            )

    return CohortSelection(
        cases=tuple(selected),
        target_size=target_size,
    )
