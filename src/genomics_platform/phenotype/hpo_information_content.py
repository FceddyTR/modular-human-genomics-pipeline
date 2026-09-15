"""Disease-corpus information content for HPO terms.

Information content is derived from unique disease-HPO annotations.

Rules:
- annotation unit = unique disease
- duplicate gene associations do not increase frequency
- only terms under HP:0000118 Phenotypic abnormality are used
- annotations are propagated to phenotype ancestors
- IC(t) = -ln(P(t))
- P(t) = propagated disease count / phenotype corpus disease count

This module provides phenotype informativeness only.
It does not perform diagnosis, variant classification, or ranking.
"""

from __future__ import annotations

import math
import sqlite3
from collections import defaultdict
from pathlib import Path

from genomics_platform.phenotype.hpo_ontology import (
    HPOOntology,
)


SCHEMA_VERSION = "1.0"
PHENOTYPE_ROOT = "HP:0000118"


def create_ic_schema(
    connection: sqlite3.Connection,
) -> None:
    connection.executescript(
        """
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );

        CREATE TABLE information_content (
            hpo_id TEXT PRIMARY KEY,
            direct_disease_count INTEGER NOT NULL,
            propagated_disease_count INTEGER NOT NULL,
            probability REAL NOT NULL,
            information_content REAL NOT NULL
        );

        CREATE INDEX idx_ic_value
        ON information_content(information_content);

        CREATE INDEX idx_ic_propagated_count
        ON information_content(propagated_disease_count);
        """
    )


def _load_unique_disease_annotations(
    association_database: str | Path,
) -> dict[str, set[str]]:
    connection = sqlite3.connect(
        str(association_database)
    )

    try:
        rows = connection.execute(
            """
            SELECT DISTINCT
                disease_id,
                hpo_id
            FROM gene_phenotype
            WHERE disease_id IS NOT NULL
              AND disease_id != ''
              AND hpo_id IS NOT NULL
              AND hpo_id != ''
            """
        ).fetchall()

    finally:
        connection.close()

    disease_terms: dict[
        str,
        set[str],
    ] = defaultdict(set)

    for disease_id, hpo_id in rows:
        disease_terms[
            str(disease_id)
        ].add(
            str(hpo_id)
        )

    return dict(disease_terms)


def build_information_content_database(
    association_database: str | Path,
    ontology_database: str | Path,
    output_database: str | Path,
) -> dict:
    association_database = Path(
        association_database
    )

    ontology_database = Path(
        ontology_database
    )

    output_database = Path(
        output_database
    )

    if not association_database.exists():
        raise FileNotFoundError(
            association_database
        )

    if not ontology_database.exists():
        raise FileNotFoundError(
            ontology_database
        )

    ontology = HPOOntology(
        ontology_database
    )

    # ancestors(root) does not give descendants.
    # Determine phenotype membership by asking whether
    # PHENOTYPE_ROOT is an ancestor of each annotated term.
    disease_annotations = (
        _load_unique_disease_annotations(
            association_database
        )
    )

    direct_term_diseases: dict[
        str,
        set[str],
    ] = defaultdict(set)

    propagated_term_diseases: dict[
        str,
        set[str],
    ] = defaultdict(set)

    phenotype_diseases = set()
    retained_direct_annotations = 0
    excluded_nonphenotype_annotations = 0
    unknown_ontology_terms = 0

    for (
        disease_id,
        annotated_terms,
    ) in disease_annotations.items():

        retained_for_disease = set()

        for hpo_id in annotated_terms:
            ancestors = ontology.ancestors(
                hpo_id,
                include_self=True,
            )

            if not ancestors:
                unknown_ontology_terms += 1
                continue

            if PHENOTYPE_ROOT not in ancestors:
                excluded_nonphenotype_annotations += 1
                continue

            retained_for_disease.add(
                hpo_id
            )

            direct_term_diseases[
                hpo_id
            ].add(
                disease_id
            )

        if not retained_for_disease:
            continue

        phenotype_diseases.add(
            disease_id
        )

        retained_direct_annotations += len(
            retained_for_disease
        )

        propagated_terms = set()

        for hpo_id in retained_for_disease:
            ancestors = ontology.ancestors(
                hpo_id,
                include_self=True,
            )

            for ancestor in ancestors:
                ancestor_ancestors = (
                    ontology.ancestors(
                        ancestor,
                        include_self=True,
                    )
                )

                if (
                    PHENOTYPE_ROOT
                    in ancestor_ancestors
                ):
                    propagated_terms.add(
                        ancestor
                    )

        for hpo_id in propagated_terms:
            propagated_term_diseases[
                hpo_id
            ].add(
                disease_id
            )

    corpus_size = len(
        phenotype_diseases
    )

    if corpus_size == 0:
        raise ValueError(
            "No phenotype diseases were retained."
        )

    output_database.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if output_database.exists():
        output_database.unlink()

    connection = sqlite3.connect(
        output_database
    )

    try:
        create_ic_schema(
            connection
        )

        rows = []

        for (
            hpo_id,
            diseases,
        ) in propagated_term_diseases.items():

            propagated_count = len(
                diseases
            )

            probability = (
                propagated_count
                / corpus_size
            )

            ic = -math.log(
                probability
            )

            rows.append(
                (
                    hpo_id,
                    len(
                        direct_term_diseases.get(
                            hpo_id,
                            set(),
                        )
                    ),
                    propagated_count,
                    probability,
                    ic,
                )
            )

        connection.executemany(
            """
            INSERT INTO information_content (
                hpo_id,
                direct_disease_count,
                propagated_disease_count,
                probability,
                information_content
            )
            VALUES (?, ?, ?, ?, ?)
            """,
            rows,
        )

        metadata = {
            "schema_version": (
                SCHEMA_VERSION
            ),
            "method": (
                "disease_level_ancestor_propagated"
            ),
            "phenotype_root": (
                PHENOTYPE_ROOT
            ),
            "association_database": (
                association_database.name
            ),
            "ontology_database": (
                ontology_database.name
            ),
            "phenotype_disease_count": (
                str(corpus_size)
            ),
            "retained_direct_annotation_count": (
                str(
                    retained_direct_annotations
                )
            ),
            "excluded_nonphenotype_annotation_count": (
                str(
                    excluded_nonphenotype_annotations
                )
            ),
            "unknown_ontology_annotation_count": (
                str(
                    unknown_ontology_terms
                )
            ),
            "ic_term_count": str(
                len(rows)
            ),
        }

        connection.executemany(
            """
            INSERT INTO metadata (
                key,
                value
            )
            VALUES (?, ?)
            """,
            tuple(
                metadata.items()
            ),
        )

        connection.commit()

    finally:
        connection.close()

    return {
        "schema_version": SCHEMA_VERSION,
        "phenotype_root": PHENOTYPE_ROOT,
        "phenotype_disease_count": (
            corpus_size
        ),
        "ic_term_count": len(rows),
        "retained_direct_annotation_count": (
            retained_direct_annotations
        ),
        "excluded_nonphenotype_annotation_count": (
            excluded_nonphenotype_annotations
        ),
        "unknown_ontology_annotation_count": (
            unknown_ontology_terms
        ),
        "database": str(
            output_database
        ),
    }
