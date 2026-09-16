import copy
import json

import pytest

from genomics_platform.interpretation.case_model import (
    CandidateCase,
)
from genomics_platform.ranking.benchmark.real_cases.replay import (
    replay_candidate_case,
)
from test_case_model import make_candidate


def test_candidate_case_round_trip():
    original = make_candidate()

    replayed = replay_candidate_case(
        original.to_dict()
    )

    assert isinstance(
        replayed,
        CandidateCase,
    )

    assert (
        replayed.to_dict()
        == original.to_dict()
    )


def test_candidate_case_survives_json_round_trip():
    original = make_candidate()

    payload = json.loads(
        json.dumps(
            original.to_dict()
        )
    )

    replayed = replay_candidate_case(
        payload
    )

    assert (
        replayed.to_dict()
        == original.to_dict()
    )


def test_candidate_replay_preserves_case_context():
    original = make_candidate()

    replayed = replay_candidate_case(
        original.to_dict()
    )

    assert (
        replayed.case_context
        == original.case_context
    )


def test_candidate_replay_preserves_gene_consistency():
    original = make_candidate()

    replayed = replay_candidate_case(
        original.to_dict()
    )

    assert (
        replayed.gene_consistency
        == original.gene_consistency
    )


def test_candidate_replay_rejects_profile_variant_mismatch():
    payload = copy.deepcopy(
        make_candidate().to_dict()
    )

    payload[
        "interpretation_profile"
    ]["variant_key"] = "1:999:A:T"

    with pytest.raises(
        ValueError,
        match="Variant identity mismatch",
    ):
        replay_candidate_case(
            payload
        )


def test_candidate_replay_rejects_serialized_variant_key_mismatch():
    payload = copy.deepcopy(
        make_candidate().to_dict()
    )

    payload["variant_key"] = (
        "1:999:A:T"
    )

    with pytest.raises(
        ValueError,
        match="variant_key",
    ):
        replay_candidate_case(
            payload
        )


def test_candidate_replay_rejects_missing_required_dimension():
    payload = copy.deepcopy(
        make_candidate().to_dict()
    )

    payload[
        "interpretation_profile"
    ]["dimensions"]["clinical"] = None

    with pytest.raises(
        ValueError,
        match="clinical",
    ):
        replay_candidate_case(
            payload
        )
