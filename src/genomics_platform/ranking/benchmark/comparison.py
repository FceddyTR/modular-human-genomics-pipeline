"""Compare ranking modes without leaking truth into ranking."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from genomics_platform.ranking.benchmark.metrics import (
    RankingMetrics,
)


@dataclass(frozen=True)
class RankingComparison:
    """Comparison of baseline and phenotype-aware ranking."""

    baseline: RankingMetrics
    phenotype: RankingMetrics

    @property
    def causal_rank_improvement(
        self,
    ) -> int | None:
        """Return positive delta when phenotype improves rank.

        Example:
            baseline causal rank = 8
            phenotype causal rank = 2

            improvement = 8 - 2 = +6

        Truth is used only by the benchmark layer.
        It is never passed into the ranking engine.
        """

        baseline_rank = (
            self.baseline.best_causal_rank
        )

        phenotype_rank = (
            self.phenotype.best_causal_rank
        )

        if (
            baseline_rank is None
            or phenotype_rank is None
        ):
            return None

        return (
            baseline_rank
            - phenotype_rank
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "baseline": (
                self.baseline.to_dict()
            ),
            "phenotype": (
                self.phenotype.to_dict()
            ),
            "causal_rank_improvement": (
                self.causal_rank_improvement
            ),
            "interpretation": {
                "positive_delta_means": (
                    "Phenotype mode moved the "
                    "best known causal candidate "
                    "closer to rank 1."
                ),
                "truth_used_for_ranking": False,
            },
        }
