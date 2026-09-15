import pytest

from genomics_platform.interpretation.interpretation_profile import (
    EvidenceDimension,
    InterpretationObservation,
    InterpretationProfile,
)


def make_profile():
    population = EvidenceDimension(
        name="population",
        status="AVAILABLE",
        observations=(
            InterpretationObservation(
                code="GNOMAD_EXOME_AF",
                label="gnomAD exome allele frequency",
                value=0.00001,
                source="gnomAD",
            ),
        ),
    )

    clinical = EvidenceDimension(
        name="clinical",
        status="AVAILABLE",
        observations=(
            InterpretationObservation(
                code="CLINVAR_SIGNIFICANCE",
                label="ClinVar submitted significance",
                value=["Uncertain_significance"],
                source="ClinVar",
            ),
        ),
    )

    functional = EvidenceDimension(
        name="functional",
        status="AVAILABLE",
        observations=(
            InterpretationObservation(
                code="VEP_CONSEQUENCE",
                label="VEP consequence",
                value=["intron_variant"],
                source="VEP",
            ),
        ),
    )

    return InterpretationProfile(
        variant_key="1:66926:AG:A",
        population=population,
        clinical=clinical,
        functional=functional,
    )


def test_profile_serializes_dimensions():
    data = make_profile().to_dict()

    assert data["schema_version"] == "1.0"
    assert data["variant_key"] == "1:66926:AG:A"

    assert (
        data["dimensions"]["population"]["status"]
        == "AVAILABLE"
    )

    assert (
        data["dimensions"]["clinical"][
            "observations"
        ][0]["value"]
        == ["Uncertain_significance"]
    )

    assert (
        data["dimensions"]["functional"][
            "observations"
        ][0]["value"]
        == ["intron_variant"]
    )


def test_profile_tracks_available_dimensions():
    profile = make_profile()

    assert profile.available_dimensions == [
        "population",
        "clinical",
        "functional",
    ]


def test_profile_tracks_error_dimensions():
    profile = InterpretationProfile(
        variant_key="1:66926:AG:A",
        population=EvidenceDimension(
            name="population",
            status="ERROR",
        ),
        clinical=EvidenceDimension(
            name="clinical",
            status="AVAILABLE",
        ),
        functional=EvidenceDimension(
            name="functional",
            status="AVAILABLE",
        ),
    )

    assert profile.error_dimensions == [
        "population"
    ]


def test_partial_dimension_is_available():
    profile = InterpretationProfile(
        variant_key="1:66926:AG:A",
        population=EvidenceDimension(
            name="population",
            status="PARTIAL",
        ),
        clinical=EvidenceDimension(
            name="clinical",
            status="UNAVAILABLE",
        ),
        functional=EvidenceDimension(
            name="functional",
            status="AVAILABLE",
        ),
    )

    assert profile.available_dimensions == [
        "population",
        "functional",
    ]


def test_invalid_dimension_status_is_rejected():
    with pytest.raises(ValueError):
        EvidenceDimension(
            name="population",
            status="PATHOGENIC",
        )


def test_future_dimensions_can_be_absent():
    data = make_profile().to_dict()

    assert data["dimensions"]["phenotype"] is None
    assert data["dimensions"]["inheritance"] is None
    assert data["dimensions"]["gene_disease"] is None


def test_profile_has_no_ranking_or_classification_fields():
    data = make_profile().to_dict()

    forbidden_keys = {
        "ranking_score",
        "tier",
        "pathogenicity_probability",
        "acmg_classification",
        "classification",
    }

    assert forbidden_keys.isdisjoint(data.keys())

    for dimension in (
        "population",
        "clinical",
        "functional",
    ):
        assert forbidden_keys.isdisjoint(
            data["dimensions"][dimension].keys()
        )


def test_observation_preserves_details():
    observation = InterpretationObservation(
        code="TEST",
        label="Test observation",
        value="example",
        source="example_source",
        details={
            "source_release": "1.0",
            "note": "descriptive only",
        },
    )

    data = observation.to_dict()

    assert data["code"] == "TEST"
    assert data["details"]["source_release"] == "1.0"
    assert data["details"]["note"] == "descriptive only"
