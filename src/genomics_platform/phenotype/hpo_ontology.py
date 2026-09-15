"""Read-only traversal utilities for the HPO ontology graph."""

from __future__ import annotations

import sqlite3
from functools import lru_cache
from pathlib import Path


ROOT_TERM = "HP:0000001"


class HPOOntology:
    """Read-only HPO ontology traversal layer."""

    def __init__(
        self,
        database_path: str | Path,
    ) -> None:
        self.database_path = str(
            database_path
        )

        if not Path(
            self.database_path
        ).exists():
            raise FileNotFoundError(
                self.database_path
            )

    def _connect(
        self,
    ) -> sqlite3.Connection:
        connection = sqlite3.connect(
            self.database_path
        )

        connection.row_factory = (
            sqlite3.Row
        )

        return connection

    def term_exists(
        self,
        hpo_id: str,
    ) -> bool:
        with self._connect() as connection:
            row = connection.execute(
                """
                SELECT 1
                FROM terms
                WHERE hpo_id = ?
                LIMIT 1
                """,
                (hpo_id,),
            ).fetchone()

        return row is not None

    def term_name(
        self,
        hpo_id: str,
    ) -> str | None:
        with self._connect() as connection:
            row = connection.execute(
                """
                SELECT name
                FROM terms
                WHERE hpo_id = ?
                """,
                (hpo_id,),
            ).fetchone()

        if row is None:
            return None

        return row["name"]

    @lru_cache(maxsize=65536)
    def parents(
        self,
        hpo_id: str,
    ) -> tuple[str, ...]:
        with self._connect() as connection:
            rows = connection.execute(
                """
                SELECT parent_id
                FROM edges
                WHERE child_id = ?
                ORDER BY parent_id
                """,
                (hpo_id,),
            ).fetchall()

        return tuple(
            row["parent_id"]
            for row in rows
        )

    @lru_cache(maxsize=65536)
    def ancestors(
        self,
        hpo_id: str,
        include_self: bool = True,
    ) -> frozenset[str]:
        if not self.term_exists(
            hpo_id
        ):
            return frozenset()

        visited = set()

        if include_self:
            visited.add(hpo_id)

        stack = list(
            self.parents(hpo_id)
        )

        while stack:
            current = stack.pop()

            if current in visited:
                continue

            visited.add(current)

            stack.extend(
                self.parents(current)
            )

        return frozenset(
            visited
        )

    def common_ancestors(
        self,
        first: str,
        second: str,
        include_self: bool = True,
    ) -> frozenset[str]:
        return (
            self.ancestors(
                first,
                include_self,
            )
            & self.ancestors(
                second,
                include_self,
            )
        )

    @lru_cache(maxsize=65536)
    def depth(
        self,
        hpo_id: str,
    ) -> int | None:
        """Shortest distance from any root-like term.

        This is a structural diagnostic metric only.
        It is NOT information content.
        """

        if not self.term_exists(
            hpo_id
        ):
            return None

        parents = self.parents(
            hpo_id
        )

        if not parents:
            return 0

        parent_depths = [
            self.depth(parent)
            for parent in parents
        ]

        valid = [
            value
            for value in parent_depths
            if value is not None
        ]

        if not valid:
            return 0

        return 1 + min(valid)

    def most_specific_common_ancestors(
        self,
        first: str,
        second: str,
    ) -> tuple[str, ...]:
        """Return structurally deepest common ancestors.

        This is NOT yet MICA.

        MICA must later be selected by information
        content, not ontology depth.
        """

        common = self.common_ancestors(
            first,
            second,
            include_self=True,
        )

        if not common:
            return ()

        depths = {
            term: self.depth(term)
            for term in common
        }

        valid = {
            term: depth
            for term, depth
            in depths.items()
            if depth is not None
        }

        if not valid:
            return ()

        maximum = max(
            valid.values()
        )

        return tuple(
            sorted(
                term
                for term, depth
                in valid.items()
                if depth == maximum
            )
        )
