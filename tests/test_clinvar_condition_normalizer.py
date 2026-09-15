#!/usr/bin/env python3

import sqlite3

from src.genomics_platform.evidence.clinvar_condition_normalizer import (
    canonicalize_clinvar_identifier,
    normalize_clinvar_condition,
)


def build_fixture_database(path):
    con = sqlite3.connect(path)

    try:
        con.executescript(
            """
            CREATE TABLE disease (
                mondo_id TEXT PRIMARY KEY,
                label TEXT
            );

            CREATE TABLE xref (
                source_id TEXT,
                mondo_id TEXT
            );

            INSERT INTO disease(
                mondo_id,
                label
            )
            VALUES (
                'MONDO:0019200',
                'retinitis pigmentosa'
            );

            INSERT INTO xref(
                source_id,
                mondo_id
            )
            VALUES
                (
                    'OMIM:268000',
                    'MONDO:0019200'
                ),
                (
                    'Orphanet:791',
                    'MONDO:0019200'
                );
            """
        )

        con.commit()

    finally:
        con.close()


def test_clinvar_condition_normalization(tmp_path):
    database = tmp_path / "mondo.sqlite"

    build_fixture_database(database)

    assert canonicalize_clinvar_identifier(
        "MONDO:MONDO:0019200"
    ) == "MONDO:0019200"

    condition = {
        "name": "Retinitis_pigmentosa",
        "identifiers": [
            "MONDO:MONDO:0019200",
            "MeSH:D012174",
            "MedGen:C0035334",
            "OMIM:268000",
            "OMIM:PS268000",
            "Orphanet:791",
        ],
    }

    result = normalize_clinvar_condition(
        str(database),
        condition,
    )

    # Raw ClinVar provenance must be preserved.
    assert result["raw_name"] == (
        "Retinitis_pigmentosa"
    )

    assert (
        "MONDO:MONDO:0019200"
        in result["raw_identifiers"]
    )

    # Three independent identifiers resolve
    # to the same canonical MONDO concept.
    normalization = result["normalization"]

    assert normalization["status"] == "FOUND"
    assert normalization["mondo_id"] == (
        "MONDO:0019200"
    )
    assert normalization["mondo_label"] == (
        "retinitis pigmentosa"
    )

    assert set(normalization["matched_via"]) == {
        "MONDO:MONDO:0019200",
        "OMIM:268000",
        "Orphanet:791",
    }

    by_raw = {
        item["raw_identifier"]: item
        for item in result["identifier_results"]
    }

    # Source-specific ClinVar cleanup is derived;
    # raw identifier remains untouched.
    assert by_raw[
        "MONDO:MONDO:0019200"
    ]["lookup_identifier"] == (
        "MONDO:0019200"
    )

    assert by_raw[
        "MONDO:MONDO:0019200"
    ]["status"] == "FOUND"

    # Unsupported/unmapped identifiers must
    # remain explicitly NOT_FOUND.
    assert by_raw[
        "MeSH:D012174"
    ]["status"] == "NOT_FOUND"

    assert by_raw[
        "MedGen:C0035334"
    ]["status"] == "NOT_FOUND"

    assert by_raw[
        "OMIM:PS268000"
    ]["status"] == "NOT_FOUND"

    # Normalization must not imply causality.
    assert "does not establish" in (
        result["semantics"]
    )


def test_no_identifier_match(tmp_path):
    database = tmp_path / "mondo.sqlite"

    build_fixture_database(database)

    result = normalize_clinvar_condition(
        str(database),
        {
            "name": "Unknown_condition",
            "identifiers": [
                "MedGen:UNKNOWN"
            ],
        },
    )

    assert (
        result["normalization"]["status"]
        == "NOT_FOUND"
    )

    assert (
        result["normalization"]["mondo_id"]
        is None
    )
