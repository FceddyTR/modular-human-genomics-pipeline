from scripts.benchmark.build_public_cohort import (
    canonical_json_bytes,
    manifest_sha256,
)


def test_canonical_json_is_order_independent():
    first = {
        "b": 2,
        "a": 1,
    }

    second = {
        "a": 1,
        "b": 2,
    }

    assert (
        canonical_json_bytes(first)
        == canonical_json_bytes(second)
    )

    assert (
        manifest_sha256(first)
        == manifest_sha256(second)
    )


def test_manifest_hash_changes_with_content():
    first = {
        "locked_evaluation": [
            "case-a"
        ]
    }

    second = {
        "locked_evaluation": [
            "case-b"
        ]
    }

    assert (
        manifest_sha256(first)
        != manifest_sha256(second)
    )
