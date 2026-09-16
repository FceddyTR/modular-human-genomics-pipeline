"""Machine-readable artifacts for ranking benchmarks.

This module belongs strictly to the benchmark layer.

Benchmark truth may be written to benchmark artifacts for evaluation
and visualization, but it must never flow back into candidate ranking,
interpretation, or evidence aggregation.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from genomics_platform.ranking.borda import (
    BordaResult,
)
from genomics_platform.ranking.rank_aggregation import (
    RankAggregationResult,
)
from genomics_platform.ranking.ranking_contract import (
    RankingResult,
)
from genomics_platform.ranking.benchmark.borda_metrics import (
    evaluate_borda_ranking,
)
from genomics_platform.ranking.benchmark.metrics import (
    evaluate_ranking,
)
from genomics_platform.ranking.benchmark.rank_aggregation_metrics import (
    evaluate_rank_aggregation,
)
from genomics_platform.ranking.benchmark.truth_set import (
    RankingTruth,
)


ARTIFACT_SCHEMA_VERSION = "1.0"

DEFAULT_K_VALUES = (
    1,
    3,
    5,
)


@dataclass(frozen=True)
class BenchmarkRun:
    """One completed ranking strategy for one benchmark case."""

    strategy: str
    result: RankingResult | BordaResult | RankAggregationResult

    def __post_init__(self) -> None:
        strategy = self.strategy.strip()

        if not strategy:
            raise ValueError(
                "strategy is required."
            )

        object.__setattr__(
            self,
            "strategy",
            strategy,
        )


@dataclass(frozen=True)
class BenchmarkCase:
    """Completed rankings plus post-ranking truth for one case."""

    truth: RankingTruth
    runs: tuple[BenchmarkRun, ...]

    def __post_init__(self) -> None:
        if not self.runs:
            raise ValueError(
                "At least one benchmark run is required."
            )

        names = [
            run.strategy
            for run in self.runs
        ]

        if len(names) != len(set(names)):
            raise ValueError(
                "Benchmark strategy names must be unique "
                "within a case."
            )


def _phenotype_mode(
    result: RankingResult | BordaResult | RankAggregationResult,
) -> str:
    if isinstance(result, RankingResult):
        return result.mode.value

    return result.source_mode


def _aggregation(
    result: RankingResult | BordaResult | RankAggregationResult,
) -> str:
    if isinstance(result, RankingResult):
        return "weighted_additive"

    if isinstance(result, BordaResult):
        return "borda"

    if isinstance(result, RankAggregationResult):
        return f"rank_{result.method.value}"

    raise TypeError(
        f"Unsupported ranking result type: {type(result)!r}"
    )


def _score_type(
    result: RankingResult | BordaResult | RankAggregationResult,
) -> str:
    if isinstance(result, RankingResult):
        return "prioritization_score"

    if isinstance(result, BordaResult):
        return "borda_score"

    if isinstance(result, RankAggregationResult):
        return "aggregate_rank_score"

    raise TypeError(
        f"Unsupported ranking result type: {type(result)!r}"
    )


def _metric_row(
    *,
    case_id: str,
    run: BenchmarkRun,
    truth: RankingTruth,
    k_values: tuple[int, ...],
) -> dict[str, Any]:
    result = run.result

    row: dict[str, Any] = {
        "case_id": case_id,
        "strategy": run.strategy,
        "aggregation": _aggregation(result),
        "phenotype_mode": _phenotype_mode(result),
        "best_causal_rank": None,
        "reciprocal_rank": 0.0,
        "top_tie_size": 1,
        "causal_top_tie": False,
    }

    if isinstance(result, RankingResult):
        metrics = evaluate_ranking(
            result,
            truth,
            k_values=k_values,
        )

        row["best_causal_rank"] = (
            metrics.best_causal_rank
        )
        row["reciprocal_rank"] = (
            metrics.reciprocal_rank
        )

        for k in k_values:
            row[f"hit_at_{k}"] = (
                metrics.hit_at_k[k]
            )

            # Weighted ranking currently has no
            # final-score tie contract. Its Hit@K
            # is therefore also the unique Hit@K
            # under the current ranking contract.
            row[f"unique_hit_at_{k}"] = (
                metrics.hit_at_k[k]
            )

        return row

    if isinstance(result, BordaResult):
        metrics = evaluate_borda_ranking(
            result,
            truth,
            k_values=k_values,
        )
    elif isinstance(result, RankAggregationResult):
        metrics = evaluate_rank_aggregation(
            result,
            truth,
            k_values=k_values,
        )
    else:
        raise TypeError(
            f"Unsupported ranking result type: {type(result)!r}"
        )

    row["best_causal_rank"] = (
        metrics.best_causal_rank
    )
    row["reciprocal_rank"] = (
        metrics.reciprocal_rank
    )
    row["top_tie_size"] = (
        metrics.top_tie_size
    )
    row["causal_top_tie"] = (
        metrics.causal_top_tie
    )

    for k in k_values:
        row[f"hit_at_{k}"] = (
            metrics.hit_at_k[k]
        )
        row[f"unique_hit_at_{k}"] = (
            metrics.unique_hit_at_k[k]
        )

    return row


def build_case_metric_rows(
    cases: Iterable[BenchmarkCase],
    *,
    k_values: Iterable[int] = DEFAULT_K_VALUES,
) -> list[dict[str, Any]]:
    """Build one benchmark metric row per case and strategy."""

    k_values = tuple(
        sorted(
            set(
                int(k)
                for k in k_values
            )
        )
    )

    if not k_values:
        raise ValueError(
            "At least one k value is required."
        )

    if any(k < 1 for k in k_values):
        raise ValueError(
            "All k values must be >= 1."
        )

    rows: list[dict[str, Any]] = []

    for case in cases:
        for run in case.runs:
            rows.append(
                _metric_row(
                    case_id=case.truth.case_id,
                    run=run,
                    truth=case.truth,
                    k_values=k_values,
                )
            )

    return rows


def build_candidate_rows(
    cases: Iterable[BenchmarkCase],
) -> list[dict[str, Any]]:
    """Build candidate-level audit rows.

    is_causal exists only in this benchmark artifact layer.
    It is never part of the ranking engine input.
    """

    rows: list[dict[str, Any]] = []

    for case in cases:
        truth_ids = set(
            case.truth.causal_candidate_ids
        )

        for run in case.runs:
            result = run.result
            aggregation = _aggregation(
                result
            )
            phenotype_mode = _phenotype_mode(
                result
            )
            score_type = _score_type(
                result
            )

            for candidate in result.candidates:
                if isinstance(
                    result,
                    RankingResult,
                ):
                    score = (
                        candidate.prioritization_score
                    )
                    participated = None
                    available = None
                    coverage = None
                elif isinstance(
                    result,
                    BordaResult,
                ):
                    score = candidate.borda_score
                    participated = (
                        candidate.ballots_participated
                    )
                    available = (
                        candidate.ballots_available
                    )
                    coverage = (
                        participated / available
                        if available
                        else 0.0
                    )
                elif isinstance(
                    result,
                    RankAggregationResult,
                ):
                    score = (
                        candidate.aggregate_rank_score
                    )
                    participated = (
                        candidate.components_participated
                    )
                    available = (
                        candidate.components_available
                    )
                    coverage = (
                        candidate.evidence_coverage
                    )
                else:
                    raise TypeError(
                        "Unsupported ranking result type: "
                        f"{type(result)!r}"
                    )

                rows.append(
                    {
                        "case_id": (
                            case.truth.case_id
                        ),
                        "strategy": run.strategy,
                        "aggregation": aggregation,
                        "phenotype_mode": (
                            phenotype_mode
                        ),
                        "candidate_id": (
                            candidate.candidate_id
                        ),
                        "variant_key": (
                            candidate.variant_key
                        ),
                        "gene_symbol": (
                            candidate.gene_symbol
                        ),
                        "rank": candidate.rank,
                        "score": score,
                        "score_type": score_type,
                        "is_causal": (
                            candidate.candidate_id
                            in truth_ids
                        ),
                        "ballots_participated": (
                            participated
                        ),
                        "ballots_available": (
                            available
                        ),
                        "evidence_coverage": (
                            coverage
                        ),
                    }
                )

    return rows


def build_strategy_comparison_rows(
    metric_rows: Iterable[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Build within-case pairwise strategy comparisons.

    Positive causal_rank_delta means the comparison strategy moved
    the best causal candidate closer to rank 1 than the baseline.
    """

    rows = list(metric_rows)

    grouped: dict[
        str,
        list[dict[str, Any]],
    ] = {}

    for row in rows:
        grouped.setdefault(
            str(row["case_id"]),
            [],
        ).append(row)

    comparisons: list[dict[str, Any]] = []

    for case_id, case_rows in grouped.items():
        ordered = sorted(
            case_rows,
            key=lambda row: str(
                row["strategy"]
            ),
        )

        for baseline in ordered:
            for comparison in ordered:
                if (
                    baseline["strategy"]
                    == comparison["strategy"]
                ):
                    continue

                baseline_rank = (
                    baseline["best_causal_rank"]
                )
                comparison_rank = (
                    comparison["best_causal_rank"]
                )

                rank_delta = None

                if (
                    baseline_rank is not None
                    and comparison_rank is not None
                ):
                    rank_delta = (
                        baseline_rank
                        - comparison_rank
                    )

                comparisons.append(
                    {
                        "case_id": case_id,
                        "baseline_strategy": (
                            baseline["strategy"]
                        ),
                        "comparison_strategy": (
                            comparison["strategy"]
                        ),
                        "baseline_causal_rank": (
                            baseline_rank
                        ),
                        "comparison_causal_rank": (
                            comparison_rank
                        ),
                        "causal_rank_delta": (
                            rank_delta
                        ),
                        "baseline_mrr": (
                            baseline[
                                "reciprocal_rank"
                            ]
                        ),
                        "comparison_mrr": (
                            comparison[
                                "reciprocal_rank"
                            ]
                        ),
                        "mrr_delta": (
                            comparison[
                                "reciprocal_rank"
                            ]
                            - baseline[
                                "reciprocal_rank"
                            ]
                        ),
                    }
                )

    return comparisons


def _write_tsv(
    path: Path,
    rows: list[dict[str, Any]],
) -> None:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not rows:
        path.write_text(
            "",
            encoding="utf-8",
        )
        return

    fieldnames = list(
        rows[0].keys()
    )

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="raise",
        )

        writer.writeheader()
        writer.writerows(rows)


def export_benchmark_artifacts(
    cases: Iterable[BenchmarkCase],
    output_dir: str | Path,
    *,
    k_values: Iterable[int] = DEFAULT_K_VALUES,
) -> dict[str, Path]:
    """Export machine-readable benchmark artifacts."""

    cases = tuple(cases)

    if not cases:
        raise ValueError(
            "At least one benchmark case is required."
        )

    case_ids = [
        case.truth.case_id
        for case in cases
    ]

    if len(case_ids) != len(set(case_ids)):
        raise ValueError(
            "Benchmark case IDs must be unique."
        )

    k_values = tuple(
        sorted(
            set(
                int(k)
                for k in k_values
            )
        )
    )

    metric_rows = build_case_metric_rows(
        cases,
        k_values=k_values,
    )

    candidate_rows = build_candidate_rows(
        cases
    )

    comparison_rows = (
        build_strategy_comparison_rows(
            metric_rows
        )
    )

    root = Path(
        output_dir
    )

    tables = root / "tables"

    case_metrics_path = (
        tables / "case_metrics.tsv"
    )

    candidate_rankings_path = (
        tables / "candidate_rankings.tsv"
    )

    comparison_path = (
        tables / "strategy_comparison.tsv"
    )

    summary_path = (
        root / "summary.json"
    )

    _write_tsv(
        case_metrics_path,
        metric_rows,
    )

    _write_tsv(
        candidate_rankings_path,
        candidate_rows,
    )

    _write_tsv(
        comparison_path,
        comparison_rows,
    )

    summary = {
        "schema_version": (
            ARTIFACT_SCHEMA_VERSION
        ),
        "benchmark": {
            "case_count": len(cases),
            "run_count": sum(
                len(case.runs)
                for case in cases
            ),
            "k_values": list(
                k_values
            ),
        },
        "semantic_boundaries": {
            "truth_used_for_ranking": False,
            "truth_used_for_benchmark": True,
            "relative_prioritization_only": True,
            "pathogenicity_probability": False,
            "acmg_amp_classification": False,
            "variant_classification": False,
            "diagnosis": False,
        },
        "artifacts": {
            "case_metrics": (
                "tables/case_metrics.tsv"
            ),
            "candidate_rankings": (
                "tables/candidate_rankings.tsv"
            ),
            "strategy_comparison": (
                "tables/strategy_comparison.tsv"
            ),
        },
    }

    root.mkdir(
        parents=True,
        exist_ok=True,
    )

    summary_path.write_text(
        json.dumps(
            summary,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    return {
        "summary": summary_path,
        "case_metrics": (
            case_metrics_path
        ),
        "candidate_rankings": (
            candidate_rankings_path
        ),
        "strategy_comparison": (
            comparison_path
        ),
    }
