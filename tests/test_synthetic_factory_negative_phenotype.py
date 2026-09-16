from genomics_platform.ranking.benchmark.synthetic_factory import (
    make_synthetic_candidate,
)


def _phenotype_observation(candidate, code):
    matches = [
        observation
        for observation
        in candidate.interpretation_profile.phenotype.observations
        if observation.code == code
    ]

    assert len(matches) == 1
    return matches[0]


def test_absent_hpo_survives_controlled_phenotype_replacement():
    candidate = make_synthetic_candidate(
        candidate_id="candidate",
        gene="GENE1",
        position=100,
        patient_hpo="HP:0040154",
        associated_hpo="HP:0001061",
        absent_hpo_terms=("HP:0000987",),
    )

    summary = _phenotype_observation(
        candidate,
        "PHENOTYPE_EXACT_MATCH_SUMMARY",
    )

    assert summary.value["unmatched_absent"] == [
        "HP:0000987"
    ]


def test_absent_exact_match_is_preserved_as_negative_context():
    candidate = make_synthetic_candidate(
        candidate_id="candidate",
        gene="GENE1",
        position=101,
        patient_hpo="HP:0040154",
        associated_hpo="HP:0001061",
        absent_hpo_terms=("HP:0001061",),
    )

    association = _phenotype_observation(
        candidate,
        "PHENOTYPE_SOURCE_ASSOCIATION",
    )

    assert association.value[
        "absent_exact_matches"
    ] == ["HP:0001061"]


def test_negative_phenotype_does_not_remove_present_context():
    candidate = make_synthetic_candidate(
        candidate_id="candidate",
        gene="GENE1",
        position=102,
        patient_hpo="HP:0040154",
        associated_hpo="HP:0040154",
        exact_phenotype_match=True,
        absent_hpo_terms=("HP:0000987",),
    )

    summary = _phenotype_observation(
        candidate,
        "PHENOTYPE_EXACT_MATCH_SUMMARY",
    )

    assert summary.value["matched_present"] == [
        "HP:0040154"
    ]

    assert summary.value["unmatched_absent"] == [
        "HP:0000987"
    ]
