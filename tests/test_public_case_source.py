from pathlib import Path

from genomics_platform.ranking.benchmark.public_cases import (
    publication_id_from_path,
)


def test_publication_id_from_filename():
    path = Path(
        "AAGAB/"
        "PMID_24573067_CASE_REPORT.json"
    )

    assert publication_id_from_path(
        path
    ) == "PMID:24573067"


def test_publication_fallback_is_deterministic():
    path = Path(
        "GENE1/case_without_pmid.json"
    )

    assert publication_id_from_path(
        path
    ) == "SOURCE:GENE1"
