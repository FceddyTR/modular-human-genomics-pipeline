"""HPO ontology SQLite store.

This database is intentionally separate from the existing
gene-phenotype association database.

It stores ontology structure only:
- HPO terms
- direct is_a relationships
- ontology metadata

It does not perform phenotype ranking or diagnosis.
"""

from __future__ import annotations

import re
import sqlite3
from pathlib import Path


SCHEMA_VERSION = "1.0"

HPO_RE = re.compile(r"^HP:\d{7}$")


def create_schema(
    connection: sqlite3.Connection,
) -> None:
    connection.executescript(
        """
        PRAGMA foreign_keys = ON;

        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS terms (
            hpo_id TEXT PRIMARY KEY,
            name TEXT,
            is_obsolete INTEGER NOT NULL DEFAULT 0
        );

        CREATE TABLE IF NOT EXISTS edges (
            child_id TEXT NOT NULL,
            parent_id TEXT NOT NULL,
            relation TEXT NOT NULL DEFAULT 'is_a',

            PRIMARY KEY (
                child_id,
                parent_id,
                relation
            )
        );

        CREATE INDEX IF NOT EXISTS idx_hpo_edges_child
        ON edges(child_id);

        CREATE INDEX IF NOT EXISTS idx_hpo_edges_parent
        ON edges(parent_id);
        """
    )


def _parse_term(
    lines: list[str],
) -> dict | None:
    hpo_id = None
    name = None
    is_obsolete = False
    parents = []

    for line in lines:
        if line.startswith("id: "):
            hpo_id = line[4:].strip()

        elif line.startswith("name: "):
            name = line[6:].strip()

        elif line == "is_obsolete: true":
            is_obsolete = True

        elif line.startswith("is_a: "):
            parent = (
                line[6:]
                .split(" ! ", 1)[0]
                .strip()
            )

            if HPO_RE.match(parent):
                parents.append(parent)

    if (
        hpo_id is None
        or not HPO_RE.match(hpo_id)
    ):
        return None

    return {
        "hpo_id": hpo_id,
        "name": name,
        "is_obsolete": is_obsolete,
        "parents": tuple(
            dict.fromkeys(parents)
        ),
    }


def parse_obo(
    obo_path: str | Path,
):
    obo_path = Path(obo_path)

    current = []
    in_term = False

    with obo_path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            if line == "[Term]":
                if current:
                    term = _parse_term(
                        current
                    )

                    if term is not None:
                        yield term

                current = []
                in_term = True
                continue

            if line.startswith("[") and line.endswith("]"):
                if in_term and current:
                    term = _parse_term(
                        current
                    )

                    if term is not None:
                        yield term

                current = []
                in_term = False
                continue

            if in_term:
                current.append(line)

    if in_term and current:
        term = _parse_term(current)

        if term is not None:
            yield term


def build_ontology_database(
    obo_path: str | Path,
    database_path: str | Path,
) -> dict:
    obo_path = Path(obo_path)
    database_path = Path(database_path)

    if not obo_path.exists():
        raise FileNotFoundError(
            obo_path
        )

    database_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if database_path.exists():
        database_path.unlink()

    connection = sqlite3.connect(
        database_path
    )

    try:
        create_schema(connection)

        term_count = 0
        edge_count = 0

        for term in parse_obo(
            obo_path
        ):
            connection.execute(
                """
                INSERT INTO terms (
                    hpo_id,
                    name,
                    is_obsolete
                )
                VALUES (?, ?, ?)
                """,
                (
                    term["hpo_id"],
                    term["name"],
                    int(
                        term["is_obsolete"]
                    ),
                ),
            )

            term_count += 1

            for parent in term[
                "parents"
            ]:
                connection.execute(
                    """
                    INSERT OR IGNORE INTO edges (
                        child_id,
                        parent_id,
                        relation
                    )
                    VALUES (?, ?, 'is_a')
                    """,
                    (
                        term["hpo_id"],
                        parent,
                    ),
                )

                edge_count += 1

        connection.executemany(
            """
            INSERT INTO metadata (
                key,
                value
            )
            VALUES (?, ?)
            """,
            (
                (
                    "schema_version",
                    SCHEMA_VERSION,
                ),
                (
                    "source_file",
                    obo_path.name,
                ),
                (
                    "term_count",
                    str(term_count),
                ),
                (
                    "edge_count",
                    str(edge_count),
                ),
            ),
        )

        connection.commit()

        return {
            "schema_version": (
                SCHEMA_VERSION
            ),
            "term_count": term_count,
            "edge_count": edge_count,
            "database": str(
                database_path
            ),
        }

    finally:
        connection.close()
