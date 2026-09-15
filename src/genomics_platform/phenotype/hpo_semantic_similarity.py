"""Ontology-aware HPO semantic similarity.

Implements:
- disease-corpus information content lookup
- MICA: Most Informative Common Ancestor
- Resnik term similarity
- Best Match Average (BMA) set similarity

Semantic similarity is phenotype relevance only.
It is not diagnosis, pathogenicity, ACMG classification,
or a probability of disease.
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable

from genomics_platform.phenotype.hpo_ontology import (
    HPOOntology,
)


@dataclass(frozen=True)
class TermSimilarity:
    first: str
    second: str
    mica_terms: tuple[str, ...]
    mica_information_content: float
    similarity: float

    def to_dict(self) -> dict:
        return {
            "first": self.first,
            "second": self.second,
            "mica_terms": list(
                self.mica_terms
            ),
            "mica_information_content": (
                self.mica_information_content
            ),
            "similarity": self.similarity,
        }


@dataclass(frozen=True)
class SetSimilarity:
    patient_terms: tuple[str, ...]
    candidate_terms: tuple[str, ...]
    patient_to_candidate: float
    candidate_to_patient: float
    bma_similarity: float
    best_matches_patient: tuple[
        TermSimilarity,
        ...,
    ]
    best_matches_candidate: tuple[
        TermSimilarity,
        ...,
    ]

    def to_dict(self) -> dict:
        return {
            "patient_terms": list(
                self.patient_terms
            ),
            "candidate_terms": list(
                self.candidate_terms
            ),
            "patient_to_candidate": (
                self.patient_to_candidate
            ),
            "candidate_to_patient": (
                self.candidate_to_patient
            ),
            "bma_similarity": (
                self.bma_similarity
            ),
            "best_matches_patient": [
                item.to_dict()
                for item
                in self.best_matches_patient
            ],
            "best_matches_candidate": [
                item.to_dict()
                for item
                in self.best_matches_candidate
            ],
            "semantic_boundaries": {
                "phenotype_similarity_only": True,
                "diagnosis": False,
                "variant_classification": False,
                "pathogenicity_probability": False,
            },
        }


def _normalize_terms(
    values: Iterable[str],
) -> tuple[str, ...]:
    if isinstance(values, str):
        values = (values,)

    normalized = []

    for value in values:
        term = str(value).strip().upper()

        if term and term not in normalized:
            normalized.append(term)

    return tuple(normalized)


class HPOSemanticSimilarity:
    def __init__(
        self,
        ontology_database: str | Path,
        information_content_database: str | Path,
    ) -> None:
        self.ontology_database = str(
            ontology_database
        )

        self.information_content_database = str(
            information_content_database
        )

        if not Path(
            self.information_content_database
        ).exists():
            raise FileNotFoundError(
                self.information_content_database
            )

        self.ontology = HPOOntology(
            self.ontology_database
        )

    def _ic_connect(
        self,
    ) -> sqlite3.Connection:
        connection = sqlite3.connect(
            self.information_content_database
        )

        connection.row_factory = sqlite3.Row

        return connection

    @lru_cache(maxsize=65536)
    def information_content(
        self,
        hpo_id: str,
    ) -> float | None:
        with self._ic_connect() as connection:
            row = connection.execute(
                """
                SELECT information_content
                FROM information_content
                WHERE hpo_id = ?
                """,
                (hpo_id,),
            ).fetchone()

        if row is None:
            return None

        return float(
            row["information_content"]
        )

    @lru_cache(maxsize=131072)
    def resnik(
        self,
        first: str,
        second: str,
    ) -> TermSimilarity:
        common = (
            self.ontology.common_ancestors(
                first,
                second,
                include_self=True,
            )
        )

        scored = []

        for term in common:
            ic = self.information_content(
                term
            )

            if ic is not None:
                scored.append(
                    (
                        term,
                        ic,
                    )
                )

        if not scored:
            return TermSimilarity(
                first=first,
                second=second,
                mica_terms=(),
                mica_information_content=0.0,
                similarity=0.0,
            )

        maximum = max(
            ic
            for _, ic in scored
        )

        mica_terms = tuple(
            sorted(
                term
                for term, ic
                in scored
                if abs(
                    ic - maximum
                ) < 1e-12
            )
        )

        return TermSimilarity(
            first=first,
            second=second,
            mica_terms=mica_terms,
            mica_information_content=maximum,
            similarity=maximum,
        )

    @lru_cache(maxsize=131072)
    def lin(
        self,
        first: str,
        second: str,
    ) -> TermSimilarity:
        """Normalized Lin semantic similarity.

        sim_Lin(a,b) =
            2 * IC(MICA(a,b))
            / (IC(a) + IC(b))

        Returns a bounded phenotype similarity in [0, 1].
        """

        first_ic = self.information_content(
            first
        )

        second_ic = self.information_content(
            second
        )

        resnik = self.resnik(
            first,
            second,
        )

        if (
            first_ic is None
            or second_ic is None
        ):
            return TermSimilarity(
                first=first,
                second=second,
                mica_terms=(
                    resnik.mica_terms
                ),
                mica_information_content=(
                    resnik.mica_information_content
                ),
                similarity=0.0,
            )

        denominator = (
            first_ic
            + second_ic
        )

        if denominator <= 0:
            similarity = (
                1.0
                if first == second
                else 0.0
            )

        else:
            similarity = (
                2.0
                * resnik.mica_information_content
                / denominator
            )

        similarity = min(
            1.0,
            max(
                0.0,
                similarity,
            ),
        )

        return TermSimilarity(
            first=first,
            second=second,
            mica_terms=(
                resnik.mica_terms
            ),
            mica_information_content=(
                resnik.mica_information_content
            ),
            similarity=similarity,
        )

    def _best_match_lin(
        self,
        query: str,
        targets: tuple[str, ...],
    ) -> TermSimilarity:
        if not targets:
            return TermSimilarity(
                first=query,
                second="",
                mica_terms=(),
                mica_information_content=0.0,
                similarity=0.0,
            )

        matches = tuple(
            self.lin(
                query,
                target,
            )
            for target in targets
        )

        return max(
            matches,
            key=lambda item: (
                item.similarity,
                item.second,
            ),
        )

    def lin_bma(
        self,
        patient_terms: Iterable[str],
        candidate_terms: Iterable[str],
    ) -> SetSimilarity:
        """Best Match Average using normalized Lin similarity."""

        patient = _normalize_terms(
            patient_terms
        )

        candidate = _normalize_terms(
            candidate_terms
        )

        if not patient or not candidate:
            return SetSimilarity(
                patient_terms=patient,
                candidate_terms=candidate,
                patient_to_candidate=0.0,
                candidate_to_patient=0.0,
                bma_similarity=0.0,
                best_matches_patient=(),
                best_matches_candidate=(),
            )

        patient_best = tuple(
            self._best_match_lin(
                term,
                candidate,
            )
            for term in patient
        )

        candidate_best = tuple(
            self._best_match_lin(
                term,
                patient,
            )
            for term in candidate
        )

        patient_score = (
            sum(
                item.similarity
                for item in patient_best
            )
            / len(patient_best)
        )

        candidate_score = (
            sum(
                item.similarity
                for item in candidate_best
            )
            / len(candidate_best)
        )

        bma_score = (
            patient_score
            + candidate_score
        ) / 2.0

        return SetSimilarity(
            patient_terms=patient,
            candidate_terms=candidate,
            patient_to_candidate=patient_score,
            candidate_to_patient=candidate_score,
            bma_similarity=bma_score,
            best_matches_patient=patient_best,
            best_matches_candidate=candidate_best,
        )


    def _best_match(
        self,
        query: str,
        targets: tuple[str, ...],
    ) -> TermSimilarity:
        if not targets:
            return TermSimilarity(
                first=query,
                second="",
                mica_terms=(),
                mica_information_content=0.0,
                similarity=0.0,
            )

        matches = tuple(
            self.resnik(
                query,
                target,
            )
            for target in targets
        )

        return max(
            matches,
            key=lambda item: (
                item.similarity,
                item.second,
            ),
        )

    def bma(
        self,
        patient_terms: Iterable[str],
        candidate_terms: Iterable[str],
    ) -> SetSimilarity:
        patient = _normalize_terms(
            patient_terms
        )

        candidate = _normalize_terms(
            candidate_terms
        )

        if not patient or not candidate:
            return SetSimilarity(
                patient_terms=patient,
                candidate_terms=candidate,
                patient_to_candidate=0.0,
                candidate_to_patient=0.0,
                bma_similarity=0.0,
                best_matches_patient=(),
                best_matches_candidate=(),
            )

        patient_best = tuple(
            self._best_match(
                term,
                candidate,
            )
            for term in patient
        )

        candidate_best = tuple(
            self._best_match(
                term,
                patient,
            )
            for term in candidate
        )

        patient_score = (
            sum(
                item.similarity
                for item in patient_best
            )
            / len(patient_best)
        )

        candidate_score = (
            sum(
                item.similarity
                for item in candidate_best
            )
            / len(candidate_best)
        )

        bma_score = (
            patient_score
            + candidate_score
        ) / 2.0

        return SetSimilarity(
            patient_terms=patient,
            candidate_terms=candidate,
            patient_to_candidate=(
                patient_score
            ),
            candidate_to_patient=(
                candidate_score
            ),
            bma_similarity=bma_score,
            best_matches_patient=(
                patient_best
            ),
            best_matches_candidate=(
                candidate_best
            ),
        )
