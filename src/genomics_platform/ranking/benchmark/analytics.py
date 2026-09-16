"""Suite-level analytics for controlled ranking benchmarks.

This module consumes frozen benchmark artifacts. It does not rerun
ranking and does not introduce benchmark truth into prioritization.
"""

from __future__ import annotations

import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path


ANALYTICS_SCHEMA_VERSION = "1.0"

STRATEGY_ORDER = (
    "weighted_variant_first",
    "weighted_exact",
    "weighted_semantic",
    "borda_semantic",
    "rank_geomean_semantic",
)

PAIRWISE_COMPARISONS = (
    ("weighted_variant_first", "weighted_exact"),
    ("weighted_variant_first", "weighted_semantic"),
    ("weighted_exact", "weighted_semantic"),
    ("weighted_semantic", "borda_semantic"),
    ("weighted_semantic", "rank_geomean_semantic"),
    ("borda_semantic", "rank_geomean_semantic"),
)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )


def _write_tsv(
    path: Path,
    rows: list[dict],
    fieldnames: tuple[str, ...],
) -> None:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
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
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def _as_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def _as_int(value: str) -> int:
    return int(value)


def _as_float(value: str) -> float:
    return float(value)


def _mean(values) -> float:
    values = list(values)
    return (
        statistics.fmean(values)
        if values
        else 0.0
    )


def _median(values) -> float:
    values = list(values)
    return (
        float(statistics.median(values))
        if values
        else 0.0
    )


def _round(value: float) -> float:
    return round(float(value), 6)


def build_suite_metric_rows(
    case_metrics: list[dict[str, str]],
) -> list[dict]:
    grouped = defaultdict(list)

    for row in case_metrics:
        grouped[row["strategy"]].append(row)

    output = []

    for strategy in STRATEGY_ORDER:
        rows = grouped.get(strategy, [])

        if not rows:
            continue

        ranks = [
            _as_int(row["best_causal_rank"])
            for row in rows
        ]

        output.append(
            {
                "strategy": strategy,
                "case_count": len(rows),
                "hit_at_1": _round(
                    _mean(
                        _as_bool(row["hit_at_1"])
                        for row in rows
                    )
                ),
                "hit_at_3": _round(
                    _mean(
                        _as_bool(row["hit_at_3"])
                        for row in rows
                    )
                ),
                "hit_at_5": _round(
                    _mean(
                        _as_bool(row["hit_at_5"])
                        for row in rows
                    )
                ),
                "unique_hit_at_1": _round(
                    _mean(
                        _as_bool(
                            row["unique_hit_at_1"]
                        )
                        for row in rows
                    )
                ),
                "mrr": _round(
                    _mean(
                        _as_float(
                            row["reciprocal_rank"]
                        )
                        for row in rows
                    )
                ),
                "mean_causal_rank": _round(
                    _mean(ranks)
                ),
                "median_causal_rank": _round(
                    _median(ranks)
                ),
                "top_tie_frequency": _round(
                    _mean(
                        _as_int(
                            row["top_tie_size"]
                        ) > 1
                        for row in rows
                    )
                ),
                "causal_top_tie_frequency": _round(
                    _mean(
                        _as_bool(
                            row["causal_top_tie"]
                        )
                        for row in rows
                    )
                ),
            }
        )

    return output


def build_family_metric_rows(
    case_metrics: list[dict[str, str]],
    case_metadata: list[dict[str, str]],
) -> list[dict]:
    family_by_case = {
        row["case_id"]: row["family"]
        for row in case_metadata
    }

    grouped = defaultdict(list)

    for row in case_metrics:
        family = family_by_case[
            row["case_id"]
        ]
        grouped[
            (family, row["strategy"])
        ].append(row)

    output = []

    for family, strategy in sorted(
        grouped,
        key=lambda key: (
            key[0],
            STRATEGY_ORDER.index(key[1])
            if key[1] in STRATEGY_ORDER
            else 999,
            key[1],
        ),
    ):
        rows = grouped[(family, strategy)]

        ranks = [
            _as_int(row["best_causal_rank"])
            for row in rows
        ]

        output.append(
            {
                "family": family,
                "strategy": strategy,
                "case_count": len(rows),
                "hit_at_1": _round(
                    _mean(
                        _as_bool(row["hit_at_1"])
                        for row in rows
                    )
                ),
                "unique_hit_at_1": _round(
                    _mean(
                        _as_bool(
                            row["unique_hit_at_1"]
                        )
                        for row in rows
                    )
                ),
                "mrr": _round(
                    _mean(
                        _as_float(
                            row["reciprocal_rank"]
                        )
                        for row in rows
                    )
                ),
                "mean_causal_rank": _round(
                    _mean(ranks)
                ),
                "median_causal_rank": _round(
                    _median(ranks)
                ),
                "top_tie_frequency": _round(
                    _mean(
                        _as_int(
                            row["top_tie_size"]
                        ) > 1
                        for row in rows
                    )
                ),
            }
        )

    return output


def build_pairwise_rows(
    case_metrics: list[dict[str, str]],
) -> list[dict]:
    by_case = defaultdict(dict)

    for row in case_metrics:
        by_case[row["case_id"]][
            row["strategy"]
        ] = row

    output = []

    for baseline, comparison in PAIRWISE_COMPARISONS:
        wins = 0
        ties = 0
        losses = 0
        rank_deltas = []
        mrr_deltas = []

        for strategies in by_case.values():
            if (
                baseline not in strategies
                or comparison not in strategies
            ):
                continue

            baseline_row = strategies[baseline]
            comparison_row = strategies[
                comparison
            ]

            baseline_rank = _as_int(
                baseline_row["best_causal_rank"]
            )
            comparison_rank = _as_int(
                comparison_row[
                    "best_causal_rank"
                ]
            )

            delta = (
                baseline_rank
                - comparison_rank
            )
            rank_deltas.append(delta)

            mrr_delta = (
                _as_float(
                    comparison_row[
                        "reciprocal_rank"
                    ]
                )
                - _as_float(
                    baseline_row[
                        "reciprocal_rank"
                    ]
                )
            )
            mrr_deltas.append(mrr_delta)

            if delta > 0:
                wins += 1
            elif delta < 0:
                losses += 1
            else:
                ties += 1

        case_count = (
            wins + ties + losses
        )

        output.append(
            {
                "baseline_strategy": baseline,
                "comparison_strategy": comparison,
                "case_count": case_count,
                "comparison_wins": wins,
                "ties": ties,
                "comparison_losses": losses,
                "mean_causal_rank_delta": _round(
                    _mean(rank_deltas)
                ),
                "mean_mrr_delta": _round(
                    _mean(mrr_deltas)
                ),
            }
        )

    return output


def build_failure_rows(
    case_metrics: list[dict[str, str]],
    case_metadata: list[dict[str, str]],
) -> list[dict]:
    family_by_case = {
        row["case_id"]: row["family"]
        for row in case_metadata
    }

    output = []

    for row in case_metrics:
        rank = _as_int(
            row["best_causal_rank"]
        )

        if rank <= 1:
            continue

        output.append(
            {
                "case_id": row["case_id"],
                "family": family_by_case[
                    row["case_id"]
                ],
                "strategy": row["strategy"],
                "best_causal_rank": rank,
                "reciprocal_rank": _as_float(
                    row["reciprocal_rank"]
                ),
                "hit_at_1": _as_bool(
                    row["hit_at_1"]
                ),
                "unique_hit_at_1": _as_bool(
                    row["unique_hit_at_1"]
                ),
                "top_tie_size": _as_int(
                    row["top_tie_size"]
                ),
                "causal_top_tie": _as_bool(
                    row["causal_top_tie"]
                ),
            }
        )

    return sorted(
        output,
        key=lambda row: (
            row["family"],
            row["case_id"],
            STRATEGY_ORDER.index(
                row["strategy"]
            )
            if row["strategy"]
            in STRATEGY_ORDER
            else 999,
        ),
    )


def build_borda_coverage_summary(
    candidate_rankings: list[dict[str, str]],
) -> dict:
    rows = [
        row
        for row in candidate_rankings
        if row["strategy"]
        == "borda_semantic"
    ]

    coverages = [
        _as_float(row["evidence_coverage"])
        for row in rows
        if row.get(
            "evidence_coverage",
            "",
        ).strip()
    ]

    causal_coverages = [
        _as_float(row["evidence_coverage"])
        for row in rows
        if _as_bool(row["is_causal"])
        and row.get(
            "evidence_coverage",
            "",
        ).strip()
    ]

    return {
        "candidate_count": len(rows),
        "mean_evidence_coverage": _round(
            _mean(coverages)
        ),
        "mean_causal_evidence_coverage": _round(
            _mean(causal_coverages)
        ),
    }


def generate_suite_analytics(
    benchmark_dir: str | Path,
) -> dict[str, Path]:
    root = Path(benchmark_dir)

    case_metrics = _read_tsv(
        root / "tables" / "case_metrics.tsv"
    )
    candidate_rankings = _read_tsv(
        root
        / "tables"
        / "candidate_rankings.tsv"
    )
    case_metadata = _read_tsv(
        root / "case_metadata.tsv"
    )

    suite_rows = build_suite_metric_rows(
        case_metrics
    )
    family_rows = build_family_metric_rows(
        case_metrics,
        case_metadata,
    )
    pairwise_rows = build_pairwise_rows(
        case_metrics
    )
    failure_rows = build_failure_rows(
        case_metrics,
        case_metadata,
    )
    borda_coverage = (
        build_borda_coverage_summary(
            candidate_rankings
        )
    )

    analytics_dir = root / "analytics"
    analytics_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    suite_path = (
        analytics_dir / "suite_metrics.tsv"
    )
    family_path = (
        analytics_dir / "family_metrics.tsv"
    )
    pairwise_path = (
        analytics_dir / "strategy_pairwise.tsv"
    )
    failure_path = (
        analytics_dir / "failure_cases.tsv"
    )
    json_path = (
        analytics_dir / "suite_analytics.json"
    )

    _write_tsv(
        suite_path,
        suite_rows,
        (
            "strategy",
            "case_count",
            "hit_at_1",
            "hit_at_3",
            "hit_at_5",
            "unique_hit_at_1",
            "mrr",
            "mean_causal_rank",
            "median_causal_rank",
            "top_tie_frequency",
            "causal_top_tie_frequency",
        ),
    )

    _write_tsv(
        family_path,
        family_rows,
        (
            "family",
            "strategy",
            "case_count",
            "hit_at_1",
            "unique_hit_at_1",
            "mrr",
            "mean_causal_rank",
            "median_causal_rank",
            "top_tie_frequency",
        ),
    )

    _write_tsv(
        pairwise_path,
        pairwise_rows,
        (
            "baseline_strategy",
            "comparison_strategy",
            "case_count",
            "comparison_wins",
            "ties",
            "comparison_losses",
            "mean_causal_rank_delta",
            "mean_mrr_delta",
        ),
    )

    _write_tsv(
        failure_path,
        failure_rows,
        (
            "case_id",
            "family",
            "strategy",
            "best_causal_rank",
            "reciprocal_rank",
            "hit_at_1",
            "unique_hit_at_1",
            "top_tie_size",
            "causal_top_tie",
        ),
    )

    payload = {
        "schema_version": (
            ANALYTICS_SCHEMA_VERSION
        ),
        "case_count": len(
            {
                row["case_id"]
                for row in case_metrics
            }
        ),
        "strategy_count": len(
            {
                row["strategy"]
                for row in case_metrics
            }
        ),
        "suite_metrics": suite_rows,
        "family_metrics": family_rows,
        "pairwise_comparisons": pairwise_rows,
        "borda_coverage": borda_coverage,
        "failure_case_count": len(
            failure_rows
        ),
        "semantic_boundaries": {
            "synthetic_benchmark": True,
            "clinical_validation": False,
            "truth_used_for_ranking": False,
            "truth_used_for_benchmark": True,
            "pathogenicity_probability": False,
            "acmg_classification": False,
            "diagnosis": False,
        },
    }

    json_path.write_text(
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    return {
        "suite_metrics": suite_path,
        "family_metrics": family_path,
        "strategy_pairwise": pairwise_path,
        "failure_cases": failure_path,
        "suite_analytics": json_path,
    }
