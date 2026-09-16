import copy

import pytest

from genomics_platform.evidence.variant_evidence_contract import (
    VariantEvidenceContract,
)
from genomics_platform.ranking.benchmark.real_cases.replay import (
    replay_variant_evidence_contract,
)
from test_case_model import make_contract


def test_variant_evidence_round_trip():
    original = make_contract()

    replayed = replay_variant_evidence_contract(
        original.to_dict()
    )

    assert isinstance(
        replayed,
        VariantEvidenceContract,
    )

    assert (
        replayed.to_dict()
        == original.to_dict()
    )


def test_replay_preserves_variant_identity():
    original = make_contract()

    replayed = replay_variant_evidence_contract(
        original.to_dict()
    )

    assert (
        replayed.identity.variant_key
        == original.identity.variant_key
    )


def test_replay_preserves_availability_semantics():
    original = make_contract()

    replayed = replay_variant_evidence_contract(
        original.to_dict()
    )

    assert (
        replayed.evidence_complete
        == original.evidence_complete
    )

    assert (
        replayed.error_sources
        == original.error_sources
    )


def test_replay_rejects_serialized_variant_key_mismatch():
    payload = copy.deepcopy(
        make_contract().to_dict()
    )

    payload["variant_key"] = (
        "GRCh38:99:999:A:T"
    )

    with pytest.raises(
        ValueError,
        match="variant_key",
    ):
        replay_variant_evidence_contract(
            payload
        )


def test_replay_rejects_invalid_availability_shape():
    payload = copy.deepcopy(
        make_contract().to_dict()
    )

    payload["availability"] = "bad"

    with pytest.raises(
        ValueError,
        match="availability",
    ):
        replay_variant_evidence_contract(
            payload
        )


def test_variant_evidence_survives_json_round_trip():
    import json

    original = make_contract()

    payload = json.loads(
        json.dumps(
            original.to_dict()
        )
    )

    replayed = replay_variant_evidence_contract(
        payload
    )

    assert (
        replayed.to_dict()
        == original.to_dict()
    )
