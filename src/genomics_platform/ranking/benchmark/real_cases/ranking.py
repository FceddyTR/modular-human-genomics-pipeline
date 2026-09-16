"""Truth-blind ranking execution for frozen real-case benchmarks."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from genomics_platform.ranking.borda import (
    BordaResult,
    aggregate_borda,
)
from genomics_platform.ranking.rank_aggregation import (
    RankAggregationMethod,
    RankAggregationResult,
    aggregate_component_ranks,
)
from genomics_platform.ranking.ranking_contract import (
    RankingResult,
)
from genomics_platform.ranking.ranking_engine import (
    rank_candidates,
)
from genomics_platform.ranking.ranking_modes import (
    RankingMode,
)

from .contract import RealCaseInput


REAL_CASE_RANKING_VERSION = "0.1.0"


@dataclass(frozen=True)
class RealCaseRankingRun:
    strategy: str
    result: RankingResult | BordaResult | RankAggregationResult

    def __post_init__(self) -> None:
        if not self.strategy.strip():
            raise ValueError(
                "strategy is required."
            )


@dataclass(frozen=True)
class RealCaseRankingBundle:
    """Frozen truth-blind ranking outputs for one real case."""

    case_id: str
    runs: tuple[RealCaseRankingRun, ...]
    runner_version: str = REAL_CASE_RANKING_VERSION

    def __post_init__(self) -> None:
        if not self.case_id.strip():
            raise ValueError(
                "case_id is required."
            )

        strategies = tuple(
            run.strategy
            for run in self.runs
        )

        if len(strategies) != len(set(strategies)):
            raise ValueError(
                "Ranking strategies must be unique."
            )

    def result_for(
        self,
        strategy: str,
    ) -> RankingResult | BordaResult | RankAggregationResult:
        for run in self.runs:
            if run.strategy == strategy:
                return run.result

        raise KeyError(strategy)

    def to_dict(self) -> dict[str, Any]:
        return {
            "case_id": self.case_id,
            "runner_version": self.runner_version,
            "runs": [
                {
                    "strategy": run.strategy,
                    "result": run.result.to_dict(),
                }
                for run in self.runs
            ],
            "semantic_boundaries": {
                "truth_used_for_ranking": False,
                "benchmark_truth_present": False,
                "relative_prioritization_only": True,
                "pathogenicity_probability": False,
                "acmg_classification": False,
                "diagnosis": False,
                "clinical_validation": False,
            },
        }


def rank_real_case(
    case: RealCaseInput,
    *,
    semantic_engine,
) -> RealCaseRankingBundle:
    """Run all benchmark ranking strategies without benchmark truth."""

    candidates = case.ranking_candidates()

    weighted_variant_first = rank_candidates(
        candidates,
        mode=RankingMode.VARIANT_FIRST,
    )

    weighted_exact = rank_candidates(
        candidates,
        mode=RankingMode.PHENOTYPE,
    )

    weighted_semantic = rank_candidates(
        candidates,
        mode=RankingMode.PHENOTYPE_SEMANTIC,
        semantic_engine=semantic_engine,
    )

    borda_semantic = aggregate_borda(
        weighted_semantic
    )

    rank_geomean_semantic = aggregate_component_ranks(
        weighted_semantic,
        method=RankAggregationMethod.GEOMETRIC_MEAN,
    )

    return RealCaseRankingBundle(
        case_id=case.case_id,
        runs=(
            RealCaseRankingRun(
                strategy="weighted_variant_first",
                result=weighted_variant_first,
            ),
            RealCaseRankingRun(
                strategy="weighted_exact",
                result=weighted_exact,
            ),
            RealCaseRankingRun(
                strategy="weighted_semantic",
                result=weighted_semantic,
            ),
            RealCaseRankingRun(
                strategy="borda_semantic",
                result=borda_semantic,
            ),
            RealCaseRankingRun(
                strategy="rank_geomean_semantic",
                result=rank_geomean_semantic,
            ),
        ),
    )
