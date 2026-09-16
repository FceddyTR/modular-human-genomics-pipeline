import csv
import json

from genomics_platform.ranking.benchmark.analytics import (
    build_failure_rows,
    build_pairwise_rows,
    build_suite_metric_rows,
    generate_suite_analytics,
)


def _metric(
    case_id,
    strategy,
    rank,
    *,
    top_tie_size=1,
    causal_top_tie=False,
    unique_hit_at_1=None,
):
    if unique_hit_at_1 is None:
        unique_hit_at_1 = rank == 1

    return {
        "case_id": case_id,
        "strategy": strategy,
        "aggregation": "test",
        "phenotype_mode": "test",
        "best_causal_rank": str(rank),
        "reciprocal_rank": str(1.0 / rank),
        "top_tie_size": str(top_tie_size),
        "causal_top_tie": str(causal_top_tie),
        "hit_at_1": str(rank <= 1),
        "unique_hit_at_1": str(
            unique_hit_at_1
        ),
        "hit_at_3": str(rank <= 3),
        "unique_hit_at_3": str(rank <= 3),
        "hit_at_5": str(rank <= 5),
        "unique_hit_at_5": str(rank <= 5),
    }


def test_suite_metrics_preserve_borda_ties():
    rows = [
        _metric(
            "case-1",
            "borda_semantic",
            1,
            top_tie_size=3,
            causal_top_tie=True,
            unique_hit_at_1=False,
        ),
        _metric(
            "case-2",
            "borda_semantic",
            2,
        ),
    ]

    result = build_suite_metric_rows(rows)
    borda = result[0]

    assert borda["case_count"] == 2
    assert borda["hit_at_1"] == 0.5
    assert borda["unique_hit_at_1"] == 0.0
    assert borda["mrr"] == 0.75
    assert borda["top_tie_frequency"] == 0.5
    assert (
        borda["causal_top_tie_frequency"]
        == 0.5
    )


def test_pairwise_direction_is_comparison_improvement():
    rows = [
        _metric(
            "case-1",
            "weighted_exact",
            2,
        ),
        _metric(
            "case-1",
            "weighted_semantic",
            1,
        ),
        _metric(
            "case-2",
            "weighted_exact",
            1,
        ),
        _metric(
            "case-2",
            "weighted_semantic",
            1,
        ),
    ]

    pairwise = build_pairwise_rows(rows)

    row = next(
        x
        for x in pairwise
        if x["baseline_strategy"]
        == "weighted_exact"
        and x["comparison_strategy"]
        == "weighted_semantic"
    )

    assert row["comparison_wins"] == 1
    assert row["ties"] == 1
    assert row["comparison_losses"] == 0
    assert row["mean_causal_rank_delta"] == 0.5


def test_failure_rows_only_include_rank_above_one():
    metrics = [
        _metric(
            "case-good",
            "weighted_semantic",
            1,
        ),
        _metric(
            "case-fail",
            "weighted_semantic",
            2,
        ),
    ]

    metadata = [
        {
            "case_id": "case-good",
            "family": "test",
        },
        {
            "case_id": "case-fail",
            "family": "test",
        },
    ]

    rows = build_failure_rows(
        metrics,
        metadata,
    )

    assert len(rows) == 1
    assert rows[0]["case_id"] == "case-fail"


def test_generate_suite_analytics(tmp_path):
    root = tmp_path / "benchmark"
    tables = root / "tables"
    tables.mkdir(parents=True)

    metrics = [
        _metric(
            "case-1",
            strategy,
            1,
            top_tie_size=(
                2
                if strategy == "borda_semantic"
                else 1
            ),
            causal_top_tie=(
                strategy == "borda_semantic"
            ),
            unique_hit_at_1=(
                strategy != "borda_semantic"
            ),
        )
        for strategy in (
            "weighted_variant_first",
            "weighted_exact",
            "weighted_semantic",
            "borda_semantic",
        )
    ]

    with (
        tables / "case_metrics.tsv"
    ).open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=metrics[0].keys(),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(metrics)

    candidates = [
        {
            "case_id": "case-1",
            "strategy": "borda_semantic",
            "aggregation": "borda",
            "phenotype_mode": (
                "phenotype_semantic"
            ),
            "candidate_id": "causal",
            "variant_key": "1:1:A:G",
            "gene_symbol": "GENE",
            "rank": "1",
            "score": "1.0",
            "score_type": "borda_score",
            "is_causal": "True",
            "ballots_participated": "3",
            "ballots_available": "6",
            "evidence_coverage": "0.5",
        }
    ]

    with (
        tables / "candidate_rankings.tsv"
    ).open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=candidates[0].keys(),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(candidates)

    with (
        root / "case_metadata.tsv"
    ).open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "case_id",
                "family",
                "description",
                "expected_behavior",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "case_id": "case-1",
                "family": "test",
                "description": "test",
                "expected_behavior": "test",
            }
        )

    paths = generate_suite_analytics(root)

    for path in paths.values():
        assert path.exists()

    payload = json.loads(
        paths["suite_analytics"].read_text()
    )

    assert payload["case_count"] == 1
    assert payload["strategy_count"] == 4
    assert (
        payload["semantic_boundaries"][
            "truth_used_for_ranking"
        ]
        is False
    )
    assert (
        payload["semantic_boundaries"][
            "clinical_validation"
        ]
        is False
    )
