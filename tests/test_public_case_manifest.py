import pytest

from genomics_platform.ranking.benchmark.public_cases import (
    CohortSelection,
    PublicCaseEligibility,
    PublicCaseRecord,
    PublicCausalVariant,
    build_public_cohort_manifest,
    split_public_benchmark,
)


def _record(i):
    variant = PublicCausalVariant(
        gene_symbol=f"GENE{i}",
        hgnc_id=None,
        chrom="chr1",
        pos=i,
        ref="A",
        alt="G",
        variant_type="SNV",
    )

    eligibility = PublicCaseEligibility(
        case_id=f"case-{i}",
        eligible=True,
        present_hpo_terms=(
            "HP:0000001",
            "HP:0000002",
            "HP:0000003",
        ),
        absent_hpo_terms=(),
        causal_variants=(variant,),
        exclusion_reasons=(),
    )

    return PublicCaseRecord(
        eligibility=eligibility,
        publication_id=f"PMID:{i}",
        source_reference=f"GENE{i}/case-{i}.json",
    )


def test_manifest_freezes_300_200_contract():
    cases = tuple(
        _record(i)
        for i in range(1, 501)
    )

    cohort = CohortSelection(
        cases=cases,
        target_size=500,
    )

    split = split_public_benchmark(
        cohort,
        development_size=300,
        locked_size=200,
    )

    manifest = (
        build_public_cohort_manifest(
            split
        )
    )

    payload = manifest.to_dict()

    assert payload["cohort"] == {
        "total_cases": 500,
        "development_cases": 300,
        "locked_evaluation_cases": 200,
    }

    assert len(
        payload["development"]
    ) == 300

    assert len(
        payload["locked_evaluation"]
    ) == 200

    assert payload[
        "semantic_boundaries"
    ][
        "contains_benchmark_truth"
    ] is True

    assert payload[
        "semantic_boundaries"
    ][
        "ranking_input"
    ] is False


def test_manifest_rejects_wrong_split_size():
    with pytest.raises(
        ValueError,
        match="300 development",
    ):
        from genomics_platform.ranking.benchmark.public_cases import (
            PublicCohortManifest,
        )

        PublicCohortManifest(
            development=(
                _record(1),
            ),
            locked_evaluation=(
                _record(2),
            ),
        )
