"""Ranking modes for candidate prioritization."""

from __future__ import annotations

from enum import Enum


class RankingMode(str, Enum):
    """Supported candidate prioritization modes."""

    VARIANT_FIRST = "variant_first"
    PHENOTYPE = "phenotype"


def normalize_ranking_mode(
    value: str | RankingMode,
) -> RankingMode:
    if isinstance(value, RankingMode):
        return value

    if not isinstance(value, str):
        raise TypeError(
            "ranking mode must be str or RankingMode."
        )

    normalized = value.strip().lower()

    try:
        return RankingMode(normalized)
    except ValueError as exc:
        allowed = ", ".join(
            mode.value
            for mode in RankingMode
        )

        raise ValueError(
            f"Unsupported ranking mode: {value!r}. "
            f"Allowed: {allowed}"
        ) from exc
