from genomics_platform.phenotype.hpo_ontology import (
    HPOOntology,
)


DB = (
    "data/evidence/hpo_ontology/"
    "hpo_ontology.sqlite"
)


def ontology():
    return HPOOntology(DB)


def test_known_terms_exist():
    graph = ontology()

    assert graph.term_exists(
        "HP:0000001"
    )

    assert graph.term_exists(
        "HP:0000987"
    )

    assert graph.term_exists(
        "HP:0040154"
    )


def test_unknown_term():
    graph = ontology()

    assert not graph.term_exists(
        "HP:9999999"
    )

    assert (
        graph.ancestors(
            "HP:9999999"
        )
        == frozenset()
    )


def test_term_names():
    graph = ontology()

    assert (
        graph.term_name(
            "HP:0040154"
        )
        == "Acne inversa"
    )


def test_scarring_direct_parents():
    graph = ontology()

    parents = graph.parents(
        "HP:0000987"
    )

    assert (
        "HP:0011355"
        in parents
    )

    assert (
        "HP:0100699"
        in parents
    )


def test_acne_inversa_direct_parent():
    graph = ontology()

    assert (
        "HP:0001061"
        in graph.parents(
            "HP:0040154"
        )
    )


def test_self_is_ancestor_when_requested():
    graph = ontology()

    ancestors = graph.ancestors(
        "HP:0040154",
        include_self=True,
    )

    assert (
        "HP:0040154"
        in ancestors
    )


def test_self_can_be_excluded():
    graph = ontology()

    ancestors = graph.ancestors(
        "HP:0040154",
        include_self=False,
    )

    assert (
        "HP:0040154"
        not in ancestors
    )

    assert (
        "HP:0001061"
        in ancestors
    )


def test_root_reachable():
    graph = ontology()

    assert (
        "HP:0000001"
        in graph.ancestors(
            "HP:0040154"
        )
    )


def test_identical_terms_share_themselves():
    graph = ontology()

    common = graph.common_ancestors(
        "HP:0040154",
        "HP:0040154",
    )

    assert (
        "HP:0040154"
        in common
    )


def test_related_terms_have_common_ancestor():
    graph = ontology()

    common = graph.common_ancestors(
        "HP:0040154",
        "HP:0000987",
    )

    assert common

    assert (
        "HP:0000001"
        in common
    )


def test_depth_is_structural_only():
    graph = ontology()

    root_depth = graph.depth(
        "HP:0000001"
    )

    acne_depth = graph.depth(
        "HP:0040154"
    )

    assert root_depth == 0
    assert acne_depth is not None
    assert acne_depth > root_depth


def test_structural_specific_common_ancestor():
    graph = ontology()

    terms = (
        graph.most_specific_common_ancestors(
            "HP:0040154",
            "HP:0000987",
        )
    )

    assert isinstance(
        terms,
        tuple,
    )

    assert terms


def test_identical_term_is_most_specific_common_ancestor():
    graph = ontology()

    assert (
        graph.most_specific_common_ancestors(
            "HP:0040154",
            "HP:0040154",
        )
        == (
            "HP:0040154",
        )
    )
