"""Tests for machine-readable ranking benchmark artifacts."""

from __future__ import annotations

import csv
import json

from genomics_platform.ranking.benchmark.artifacts import (
    BenchmarkCase,
    BenchmarkRun,
    build_candidate_rows,
    build_case_metric_rows,
    export_benchmark_artifacts,
)
from genomics_platform.ranking.benchmark.truth_set import (
    RankingTruth,
)
from genomics_platform.ranking.borda import (
    aggregate_borda,
)
from genomics_platform.ranking.ranking_engine import (
    rank_candidates,
)
from genomics_platform.ranking.ranking_modes import (
    RankingMode,
)

from test_semantic_ranking_controlled_cohort import (
    build_cohort,
)


def build_case():
    cohort = build_cohort()

    weighted = rank_candidates(
        cohort,
        mode=RankingMode.PHENOTYPE,
    )

    borda = aggregate_borda(
        weighted
    )

    # Truth is introduced only after ranking.
    truth = RankingTruth(
        case_id="artifact-case-001",
        causal_candidate_ids=(
            "causal",
        ),
    )

    return BenchmarkCase(
        truth=truth,
        runs=(
            BenchmarkRun(
                strategy="weighted_exact",
                result=weighted,
            ),
            BenchmarkRun(
                strategy="borda_exact",
                result=borda,
            ),
        ),
    )


def test_metric_rows_cover_weighted_and_borda():
    rows = build_case_metric_rows(
        (build_case(),),
        k_values=(1, 3),
    )

    assert len(rows) == 2

    by_strategy = {
        row["strategy"]: row
        for row in rows
    }

    assert (
        by_strategy[
            "weighted_exact"
        ]["aggregation"]
        == "weighted_additive"
    )

    assert (
        by_strategy[
            "borda_exact"
        ]["aggregation"]
        == "borda"
    )

    assert (
        by_strategy[
            "weighted_exact"
        ]["best_causal_rank"]
        == 2
    )

    assert (
        by_strategy[
            "borda_exact"
        ]["best_causal_rank"]
        == 2
    )


def test_truth_exists_only_in_benchmark_rows():
    case = build_case()

    rows = build_candidate_rows(
        (case,)
    )

    causal = [
        row
        for row in rows
        if row["candidate_id"]
        == "causal"
    ]

    assert causal

    assert all(
        row["is_causal"] is True
        for row in causal
    )

    for run in case.runs:
        data = run.result.to_dict()

        assert "truth" not in data
        assert "is_causal" not in data


def test_borda_candidate_rows_include_coverage():
    rows = build_candidate_rows(
        (build_case(),)
    )

    borda_rows = [
        row
        for row in rows
        if row["aggregation"]
        == "borda"
    ]

    assert borda_rows

    for row in borda_rows:
        assert (
            row["ballots_participated"]
            is not None
        )

        assert (
            row["ballots_available"]
            is not None
        )

        assert (
            0.0
            <= row["evidence_coverage"]
            <= 1.0
        )


def test_export_writes_expected_artifacts(
    tmp_path,
):
    paths = export_benchmark_artifacts(
        (build_case(),),
        tmp_path,
        k_values=(1, 3),
    )

    assert paths["summary"].exists()
    assert paths["case_metrics"].exists()
    assert (
        paths["candidate_rankings"]
        .exists()
    )
    assert (
        paths["strategy_comparison"]
        .exists()
    )

    summary = json.loads(
        paths["summary"].read_text(
            encoding="utf-8"
        )
    )

    assert (
        summary["benchmark"][
            "case_count"
        ]
        == 1
    )

    assert (
        summary[
            "semantic_boundaries"
        ]["truth_used_for_ranking"]
        is False
    )

    assert (
        summary[
            "semantic_boundaries"
        ]["truth_used_for_benchmark"]
        is True
    )


def test_case_metrics_tsv_is_tabular(
    tmp_path,
):
    paths = export_benchmark_artifacts(
        (build_case(),),
        tmp_path,
        k_values=(1, 3),
    )

    with paths[
        "case_metrics"
    ].open(
        encoding="utf-8",
        newline="",
    ) as handle:
        rows = list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )

    assert len(rows) == 2

    assert {
        row["strategy"]
        for row in rows
    } == {
        "weighted_exact",
        "borda_exact",
    }


def test_summary_does_not_claim_clinical_meaning(
    tmp_path,
):
    paths = export_benchmark_artifacts(
        (build_case(),),
        tmp_path,
    )

    summary = json.loads(
        paths["summary"].read_text(
            encoding="utf-8"
        )
    )

    boundaries = summary[
        "semantic_boundaries"
    ]

    assert (
        boundaries[
            "pathogenicity_probability"
        ]
        is False
    )

    assert (
        boundaries[
            "acmg_amp_classification"
        ]
        is False
    )

    assert (
        boundaries["diagnosis"]
        is False
    )
