"""Frozen identity contract for public benchmark v0.1."""

from __future__ import annotations


PUBLIC_BENCHMARK_V0_1_SHA256 = (
    "a8c51d43b0300e0ea12de95e92388cbb"
    "8a4cc493b2255e0e7ac833ee4a284d29"
)


def verify_public_benchmark_v0_1(
    digest: str,
) -> None:
    if digest != PUBLIC_BENCHMARK_V0_1_SHA256:
        raise ValueError(
            "Public benchmark v0.1 manifest "
            "does not match the frozen SHA-256."
        )
