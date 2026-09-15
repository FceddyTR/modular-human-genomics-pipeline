import pytest

from genomics_platform.phenotype.hpo_semantic_similarity import (
    HPOSemanticSimilarity,
)


ONTOLOGY = (
    "data/evidence/hpo_ontology/"
    "hpo_ontology.sqlite"
)

IC = (
    "data/evidence/hpo_ontology/"
    "hpo_ic.sqlite"
)


@pytest.fixture(scope="module")
def semantic():
    return HPOSemanticSimilarity(
        ontology_database=ONTOLOGY,
        information_content_database=IC,
    )


def test_root_ic_is_zero(semantic):
    assert semantic.information_content(
        "HP:0000118"
    ) == pytest.approx(
        0.0,
        abs=1e-12,
    )


def test_acne_inversa_has_positive_ic(
    semantic,
):
    assert semantic.information_content(
        "HP:0040154"
    ) > 0


def test_inheritance_not_in_ic(
    semantic,
):
    assert (
        semantic.information_content(
            "HP:0000006"
        )
        is None
    )

    assert (
        semantic.information_content(
            "HP:0000007"
        )
        is None
    )


def test_identical_resnik_equals_term_ic(
    semantic,
):
    result = semantic.resnik(
        "HP:0040154",
        "HP:0040154",
    )

    assert result.mica_terms == (
        "HP:0040154",
    )

    assert result.similarity == (
        pytest.approx(
            semantic.information_content(
                "HP:0040154"
            )
        )
    )


def test_parent_child_resnik_uses_parent_ic(
    semantic,
):
    result = semantic.resnik(
        "HP:0040154",
        "HP:0001061",
    )

    assert "HP:0001061" in (
        result.mica_terms
    )

    assert result.similarity == (
        pytest.approx(
            semantic.information_content(
                "HP:0001061"
            )
        )
    )


def test_exact_more_similar_than_parent_child(
    semantic,
):
    exact = semantic.resnik(
        "HP:0040154",
        "HP:0040154",
    )

    parent = semantic.resnik(
        "HP:0040154",
        "HP:0001061",
    )

    assert (
        exact.similarity
        > parent.similarity
    )


def test_related_skin_terms_have_nonzero_similarity(
    semantic,
):
    result = semantic.resnik(
        "HP:0040154",
        "HP:0000987",
    )

    assert result.similarity > 0
    assert result.mica_terms


def test_unknown_term_has_zero_similarity(
    semantic,
):
    result = semantic.resnik(
        "HP:0040154",
        "HP:9999999",
    )

    assert result.similarity == 0.0
    assert result.mica_terms == ()


def test_single_exact_bma_equals_resnik(
    semantic,
):
    result = semantic.bma(
        ["HP:0040154"],
        ["HP:0040154"],
    )

    expected = semantic.resnik(
        "HP:0040154",
        "HP:0040154",
    ).similarity

    assert result.bma_similarity == (
        pytest.approx(expected)
    )


def test_bma_is_symmetric(
    semantic,
):
    first = semantic.bma(
        [
            "HP:0040154",
            "HP:0000987",
        ],
        [
            "HP:0001061",
        ],
    )

    second = semantic.bma(
        [
            "HP:0001061",
        ],
        [
            "HP:0040154",
            "HP:0000987",
        ],
    )

    assert first.bma_similarity == (
        pytest.approx(
            second.bma_similarity
        )
    )


def test_exact_set_scores_above_related_set(
    semantic,
):
    exact = semantic.bma(
        [
            "HP:0040154",
            "HP:0000987",
        ],
        [
            "HP:0040154",
            "HP:0000987",
        ],
    )

    related = semantic.bma(
        [
            "HP:0040154",
            "HP:0000987",
        ],
        [
            "HP:0001061",
        ],
    )

    assert (
        exact.bma_similarity
        > related.bma_similarity
    )


def test_empty_set_is_zero(
    semantic,
):
    result = semantic.bma(
        [],
        ["HP:0040154"],
    )

    assert result.bma_similarity == 0.0


def test_semantic_boundaries(
    semantic,
):
    result = semantic.bma(
        ["HP:0040154"],
        ["HP:0040154"],
    )

    boundaries = (
        result.to_dict()[
            "semantic_boundaries"
        ]
    )

    assert (
        boundaries[
            "phenotype_similarity_only"
        ]
        is True
    )

    assert boundaries["diagnosis"] is False

    assert (
        boundaries[
            "pathogenicity_probability"
        ]
        is False
    )


def test_identical_lin_is_one(
    semantic,
):
    result = semantic.lin(
        "HP:0040154",
        "HP:0040154",
    )

    assert result.similarity == (
        pytest.approx(1.0)
    )


def test_parent_child_lin_is_bounded(
    semantic,
):
    result = semantic.lin(
        "HP:0040154",
        "HP:0001061",
    )

    assert (
        0.0
        < result.similarity
        < 1.0
    )


def test_related_lin_below_parent_child(
    semantic,
):
    parent = semantic.lin(
        "HP:0040154",
        "HP:0001061",
    )

    distant = semantic.lin(
        "HP:0040154",
        "HP:0000987",
    )

    assert (
        parent.similarity
        > distant.similarity
    )


def test_lin_is_symmetric(
    semantic,
):
    first = semantic.lin(
        "HP:0040154",
        "HP:0001061",
    )

    second = semantic.lin(
        "HP:0001061",
        "HP:0040154",
    )

    assert first.similarity == (
        pytest.approx(
            second.similarity
        )
    )


def test_unknown_lin_is_zero(
    semantic,
):
    result = semantic.lin(
        "HP:0040154",
        "HP:9999999",
    )

    assert result.similarity == 0.0


def test_exact_lin_bma_is_one(
    semantic,
):
    result = semantic.lin_bma(
        [
            "HP:0040154",
            "HP:0000987",
        ],
        [
            "HP:0040154",
            "HP:0000987",
        ],
    )

    assert result.bma_similarity == (
        pytest.approx(1.0)
    )


def test_lin_bma_is_symmetric(
    semantic,
):
    first = semantic.lin_bma(
        [
            "HP:0040154",
            "HP:0000987",
        ],
        [
            "HP:0001061",
        ],
    )

    second = semantic.lin_bma(
        [
            "HP:0001061",
        ],
        [
            "HP:0040154",
            "HP:0000987",
        ],
    )

    assert first.bma_similarity == (
        pytest.approx(
            second.bma_similarity
        )
    )


def test_lin_bma_is_bounded(
    semantic,
):
    result = semantic.lin_bma(
        [
            "HP:0040154",
            "HP:0000987",
        ],
        [
            "HP:0001061",
        ],
    )

    assert (
        0.0
        <= result.bma_similarity
        <= 1.0
    )


def test_exact_lin_bma_above_related(
    semantic,
):
    exact = semantic.lin_bma(
        [
            "HP:0040154",
            "HP:0000987",
        ],
        [
            "HP:0040154",
            "HP:0000987",
        ],
    )

    related = semantic.lin_bma(
        [
            "HP:0040154",
            "HP:0000987",
        ],
        [
            "HP:0001061",
        ],
    )

    assert (
        exact.bma_similarity
        > related.bma_similarity
    )


def test_identical_lin_is_one(
    semantic,
):
    result = semantic.lin(
        "HP:0040154",
        "HP:0040154",
    )

    assert result.similarity == (
        pytest.approx(1.0)
    )


def test_parent_child_lin_is_bounded(
    semantic,
):
    result = semantic.lin(
        "HP:0040154",
        "HP:0001061",
    )

    assert (
        0.0
        < result.similarity
        < 1.0
    )


def test_related_lin_below_parent_child(
    semantic,
):
    parent = semantic.lin(
        "HP:0040154",
        "HP:0001061",
    )

    distant = semantic.lin(
        "HP:0040154",
        "HP:0000987",
    )

    assert (
        parent.similarity
        > distant.similarity
    )


def test_lin_is_symmetric(
    semantic,
):
    first = semantic.lin(
        "HP:0040154",
        "HP:0001061",
    )

    second = semantic.lin(
        "HP:0001061",
        "HP:0040154",
    )

    assert first.similarity == (
        pytest.approx(
            second.similarity
        )
    )


def test_unknown_lin_is_zero(
    semantic,
):
    result = semantic.lin(
        "HP:0040154",
        "HP:9999999",
    )

    assert result.similarity == 0.0


def test_exact_lin_bma_is_one(
    semantic,
):
    result = semantic.lin_bma(
        [
            "HP:0040154",
            "HP:0000987",
        ],
        [
            "HP:0040154",
            "HP:0000987",
        ],
    )

    assert result.bma_similarity == (
        pytest.approx(1.0)
    )


def test_lin_bma_is_symmetric(
    semantic,
):
    first = semantic.lin_bma(
        [
            "HP:0040154",
            "HP:0000987",
        ],
        [
            "HP:0001061",
        ],
    )

    second = semantic.lin_bma(
        [
            "HP:0001061",
        ],
        [
            "HP:0040154",
            "HP:0000987",
        ],
    )

    assert first.bma_similarity == (
        pytest.approx(
            second.bma_similarity
        )
    )


def test_lin_bma_is_bounded(
    semantic,
):
    result = semantic.lin_bma(
        [
            "HP:0040154",
            "HP:0000987",
        ],
        [
            "HP:0001061",
        ],
    )

    assert (
        0.0
        <= result.bma_similarity
        <= 1.0
    )


def test_exact_lin_bma_above_related(
    semantic,
):
    exact = semantic.lin_bma(
        [
            "HP:0040154",
            "HP:0000987",
        ],
        [
            "HP:0040154",
            "HP:0000987",
        ],
    )

    related = semantic.lin_bma(
        [
            "HP:0040154",
            "HP:0000987",
        ],
        [
            "HP:0001061",
        ],
    )

    assert (
        exact.bma_similarity
        > related.bma_similarity
    )
