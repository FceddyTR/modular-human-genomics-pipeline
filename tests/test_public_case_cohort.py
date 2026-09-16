import pytest

from genomics_platform.ranking.benchmark.public_cases import (
    PublicCaseEligibility,
    PublicCaseRecord,
    PublicCausalVariant,
    select_diverse_cohort,
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
        source_reference=(
            f"{publication}/{case_id}"
        ),
    )


def test_selection_is_deterministic():
    records = [
        _record(
            "c3", "GENE2", "P2", 3
        ),
        _record(
            "c1", "GENE1", "P1", 1
        ),
        _record(
            "c2", "GENE1", "P1", 2
        ),
    ]

    first = select_diverse_cohort(
        records,
        target_size=2,
    )

    second = select_diverse_cohort(
        reversed(records),
        target_size=2,
    )

    assert [
        case.case_id
        for case in first.cases
    ] == [
        case.case_id
        for case in second.cases
    ]


def test_selection_prioritizes_gene_diversity():
    records = [
        _record(
            "a1", "GENEA", "P1", 1
        ),
        _record(
            "a2", "GENEA", "P2", 2
        ),
        _record(
            "b1", "GENEB", "P3", 3
        ),
    ]

    cohort = select_diverse_cohort(
        records,
        target_size=2,
    )

    assert cohort.unique_gene_count == 2


def test_selection_rejects_insufficient_pool():
    records = [
        _record(
            "c1", "GENE1", "P1", 1
        ),
    ]

    with pytest.raises(
        ValueError,
        match="Insufficient",
    ):
        select_diverse_cohort(
            records,
            target_size=2,
        )


def test_exact_variant_set_is_deduplicated():
    first = _record(
        "c1", "GENE1", "P1", 1
    )

    duplicate = _record(
        "c2", "GENE1", "P2", 1
    )

    with pytest.raises(
        ValueError,
        match="Insufficient",
    ):
        select_diverse_cohort(
            [first, duplicate],
            target_size=2,
        )
