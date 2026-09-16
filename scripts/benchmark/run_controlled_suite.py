"""Run the controlled synthetic ranking benchmark suite.

Ranking is completed before benchmark truth is constructed.

Outputs are machine-readable benchmark artifacts under:
results/ranking_benchmark/controlled_suite_v1/
"""

from __future__ import annotations

import csv
from pathlib import Path

from genomics_platform.phenotype.hpo_semantic_similarity import (
    HPOSemanticSimilarity,
)
from genomics_platform.ranking.benchmark.artifacts import (
    BenchmarkCase,
    BenchmarkRun,
    export_benchmark_artifacts,
)
from genomics_platform.ranking.benchmark.controlled_cases import (
    build_initial_registry,
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


OUTPUT_DIR = Path(
    "results/ranking_benchmark/controlled_suite_v1"
)

SEMANTIC_DB = Path(
    "data/evidence/hpo_ontology/"
    "hpo_ontology.sqlite"
)

IC_DB = Path(
    "data/evidence/hpo_ontology/"
    "hpo_ic.sqlite"
)


def build_benchmark_case(
    controlled_case,
    semantic_engine,
):
    """Rank first; introduce benchmark truth afterward."""

    candidates = (
        controlled_case.ranking_candidates()
    )

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

    # Borda is intentionally derived from the
    # corresponding weighted component results.
    #
    # We expose Borda-semantic here because the goal
    # of this first diagnostic run is four strategies:
    #
    # variant-first
    # exact phenotype
    # semantic phenotype
    # semantic phenotype + Borda aggregation
    borda_semantic = aggregate_borda(
        weighted_semantic
    )

    # IMPORTANT:
    # Truth is created only after every ranking result
    # above has been frozen.
    truth = RankingTruth(
        case_id=controlled_case.case_id,
        causal_candidate_ids=(
            controlled_case.causal_candidate_ids()
        ),
    )

    return BenchmarkCase(
        truth=truth,
        runs=(
            BenchmarkRun(
                strategy="weighted_variant_first",
                result=weighted_variant_first,
            ),
            BenchmarkRun(
                strategy="weighted_exact",
                result=weighted_exact,
            ),
            BenchmarkRun(
                strategy="weighted_semantic",
                result=weighted_semantic,
            ),
            BenchmarkRun(
                strategy="borda_semantic",
                result=borda_semantic,
            ),
        ),
    )


def write_case_metadata(
    controlled_cases,
    output_dir,
):
    """Write presentation metadata separately from ranking artifacts."""

    path = output_dir / "case_metadata.tsv"

    fieldnames = (
        "case_id",
        "family",
        "description",
        "expected_behavior",
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
        )

        writer.writeheader()

        for case in controlled_cases:
            writer.writerow(
                {
                    "case_id": case.case_id,
                    "family": case.family,
                    "description": (
                        case.spec.description
                    ),
                    "expected_behavior": (
                        case.spec.expected_behavior
                    ),
                }
            )

    return path


def print_diagnostic_table(
    benchmark_cases,
    controlled_cases,
):
    family_by_case = {
        case.case_id: case.family
        for case in controlled_cases
    }

    strategies = (
        "weighted_variant_first",
        "weighted_exact",
        "weighted_semantic",
        "borda_semantic",
    )

    print()
    print(
        "CONTROLLED SUITE DIAGNOSTIC"
    )
    print(
        "=" * 100
    )

    header = (
        f"{'case_id':38}"
        f"{'family':12}"
        f"{'VF':>7}"
        f"{'Exact':>8}"
        f"{'Semantic':>10}"
        f"{'Borda-S':>10}"
    )

    print(header)
    print("-" * 100)

    for case in benchmark_cases:
        causal_ids = set(
            case.truth.causal_candidate_ids
        )

        ranks = {}

        for run in case.runs:
            causal_ranks = [
                candidate.rank
                for candidate
                in run.result.candidates
                if candidate.candidate_id
                in causal_ids
            ]

            ranks[run.strategy] = (
                min(causal_ranks)
                if causal_ranks
                else None
            )

        print(
            f"{case.truth.case_id:38}"
            f"{family_by_case[case.truth.case_id]:12}"
            f"{str(ranks[strategies[0]]):>7}"
            f"{str(ranks[strategies[1]]):>8}"
            f"{str(ranks[strategies[2]]):>10}"
            f"{str(ranks[strategies[3]]):>10}"
        )

    print("=" * 100)


def main():
    registry = build_initial_registry()

    controlled_cases = (
        registry.build_all()
    )

    semantic_engine = HPOSemanticSimilarity(
        ontology_database=SEMANTIC_DB,
        information_content_database=IC_DB,
    )

    # Ranking happens here.
    benchmark_cases = tuple(
        build_benchmark_case(
            case,
            semantic_engine,
        )
        for case in controlled_cases
    )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    paths = export_benchmark_artifacts(
        benchmark_cases,
        OUTPUT_DIR,
        k_values=(1, 3, 5),
    )

    metadata_path = write_case_metadata(
        controlled_cases,
        OUTPUT_DIR,
    )

    print_diagnostic_table(
        benchmark_cases,
        controlled_cases,
    )

    print()
    print(
        f"Cases: {len(benchmark_cases)}"
    )

    print(
        "Ranking runs: "
        f"{sum(len(case.runs) for case in benchmark_cases)}"
    )

    print(
        f"Metadata: {metadata_path}"
    )

    for name, path in paths.items():
        print(
            f"{name}: {path}"
        )


if __name__ == "__main__":
    main()
