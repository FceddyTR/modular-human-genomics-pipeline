import pytest

from genomics_platform.ranking.benchmark.public_cases import (
    CohortSelection,
    PublicBenchmarkSplit,
    PublicCaseEligibility,
    PublicCaseRecord,
    PublicCausalVariant,
    split_public_benchmark,
)


def _record(
    case_id,
    gene,
    publication,
    pos,
):
    variant = PublicCausalVariant(
        gene_symbol=gene,
        hgnc_id=None,
        chrom="chr1",
        pos=pos,
        ref="A",
        alt="G",
        variant_type="SNV",
    )

    eligibility = PublicCaseEligibility(
        case_id=case_id,
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
        publication_id=publication,
        source_reference=case_id,
    )


def test_split_has_no_publication_overlap():
    cases = tuple(
        _record(
            f"c{i}",
            f"GENE{i}",
            f"P{i}",
            i,
        )
        for i in range(1, 6)
    )

    cohort = CohortSelection(
        cases=cases,
        target_size=5,
    )

    result = split_public_benchmark(
        cohort,
        development_size=3,
        locked_size=2,
    )

    development_publications = {
        case.publication_id
        for case in result.development
    }

    locked_publications = {
        case.publication_id
        for case
        in result.locked_evaluation
    }

    assert not (
        development_publications
        & locked_publications
    )


def test_split_is_deterministic():
    cases = tuple(
        _record(
            f"c{i}",
            f"GENE{i}",
            f"P{i}",
            i,
        )
        for i in range(1, 6)
    )

    cohort = CohortSelection(
        cases=cases,
        target_size=5,
    )

    first = split_public_benchmark(
        cohort,
        development_size=3,
        locked_size=2,
    )

    second = split_public_benchmark(
        cohort,
        development_size=3,
        locked_size=2,
    )

    assert first == second


def test_split_rejects_wrong_total_size():
    cases = tuple(
        _record(
            f"c{i}",
            f"GENE{i}",
            f"P{i}",
            i,
        )
        for i in range(1, 5)
    )

    cohort = CohortSelection(
        cases=cases,
        target_size=4,
    )

    with pytest.raises(
        ValueError,
        match="equal",
    ):
        split_public_benchmark(
            cohort,
            development_size=3,
            locked_size=2,
        )


def test_split_contract_rejects_publication_leakage():
    development = (
        _record(
            "dev", "GENE1", "P1", 1
        ),
    )

    locked = (
        _record(
            "locked", "GENE2", "P1", 2
        ),
    )

    with pytest.raises(
        ValueError,
        match="Publication leakage",
    ):
        PublicBenchmarkSplit(
            development=development,
            locked_evaluation=locked,
        )


def test_split_contract_rejects_variant_leakage():
    development = (
        _record(
            "dev", "GENE1", "P1", 1
        ),
    )

    locked = (
        _record(
            "locked", "GENE1", "P2", 1
        ),
    )

    with pytest.raises(
        ValueError,
        match="variant leakage",
    ):
        PublicBenchmarkSplit(
            development=development,
            locked_evaluation=locked,
        )


def test_seen_and_unseen_gene_sets():
    development = (
        _record(
            "dev1", "SHARED", "P1", 1
        ),
        _record(
            "dev2", "DEVONLY", "P2", 2
        ),
    )

    locked = (
        _record(
            "lock1", "SHARED", "P3", 3
        ),
        _record(
            "lock2", "NEWGENE", "P4", 4
        ),
    )

    result = PublicBenchmarkSplit(
        development=development,
        locked_evaluation=locked,
    )

    assert result.locked_seen_genes == {
        "SHARED"
    }

    assert result.locked_unseen_genes == {
        "NEWGENE"
    }
