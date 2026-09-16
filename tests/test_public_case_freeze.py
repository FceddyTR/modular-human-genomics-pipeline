import pytest

from genomics_platform.ranking.benchmark.public_cases import (
    PUBLIC_BENCHMARK_V0_1_SHA256,
    verify_public_benchmark_v0_1,
)


def test_frozen_public_benchmark_digest():
    assert PUBLIC_BENCHMARK_V0_1_SHA256 == (
        "a8c51d43b0300e0ea12de95e92388cbb"
        "8a4cc493b2255e0e7ac833ee4a284d29"
    )

    verify_public_benchmark_v0_1(
        PUBLIC_BENCHMARK_V0_1_SHA256
    )


def test_frozen_public_benchmark_rejects_change():
    with pytest.raises(
        ValueError,
        match="frozen SHA-256",
    ):
        verify_public_benchmark_v0_1(
            "0" * 64
        )
