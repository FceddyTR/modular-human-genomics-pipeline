from genomics_platform.evidence.variant_evidence_contract import (
    ClinicalConditionEvidence,
    ClinicalEvidence,
    ClinicalRecordEvidence,
    FunctionalEvidence,
    PopulationEvidence,
    SourceAvailability,
    VariantEvidenceContract,
    VariantIdentity,
)
from genomics_platform.interpretation.clinical_interpreter import (
    interpret_clinical,
)


def make_contract(
    *,
    status="FOUND",
    records=True,
    error=None,
    significance=("Uncertain_significance",),
):
    clinical_records = ()

    if records:
        clinical_records = (
            ClinicalRecordEvidence(
                variation_id="3385321",
                allele_id="3544463",
                vcv_accession="VCV003385321",
                clinical_significance=significance,
                conflicting_significance=(),
                review_status=(
                    "criteria_provided,_single_submitter"
                ),
                review_stars=1,
                conditions=(
                    ClinicalConditionEvidence(
                        name="Retinitis_pigmentosa",
                        identifiers=(
                            "MONDO:MONDO:0019200",
                            "OMIM:268000",
                            "Orphanet:791",
                        ),
                    ),
                ),
                hgvs=(
                    "NC_000001.11:g.66927del",
                ),
                scv_accessions=(
                    "SCV005419006",
                ),
            ),
        )

    return VariantEvidenceContract(
        identity=VariantIdentity(
            assembly="GRCh38",
            chrom="1",
            pos=66926,
            ref="AG",
            alt="A",
        ),
        functional=FunctionalEvidence(),
        clinical=ClinicalEvidence(
            records=clinical_records
        ),
        population=PopulationEvidence(),
        availability=(
            SourceAvailability(
                source="VEP",
                status="FOUND",
            ),
            SourceAvailability(
                source="ClinVar",
                status=status,
                error=error,
            ),
            SourceAvailability(
                source="gnomAD",
                status="NOT_FOUND",
            ),
        ),
    )


def test_found_clinvar_is_available():
    dimension = interpret_clinical(
        make_contract()
    )

    assert dimension.status == "AVAILABLE"

    codes = {
        observation.code
        for observation in dimension.observations
    }

    assert "CLINVAR_SUBMITTED_SIGNIFICANCE" in codes
    assert "CLINVAR_REVIEW_STATUS" in codes
    assert "CLINVAR_CONDITIONS" in codes
    assert "CLINVAR_HGVS" in codes


def test_uncertain_significance_is_preserved():
    data = interpret_clinical(
        make_contract()
    ).to_dict()

    significance = next(
        observation
        for observation in data["observations"]
        if observation["code"]
        == "CLINVAR_SUBMITTED_SIGNIFICANCE"
    )

    assert significance["value"] == [
        "Uncertain_significance"
    ]

    assert (
        significance["details"][
            "platform_classification"
        ]
        is False
    )


def test_review_stars_are_metadata_not_pathogenicity():
    data = interpret_clinical(
        make_contract()
    ).to_dict()

    review = next(
        observation
        for observation in data["observations"]
        if observation["code"]
        == "CLINVAR_REVIEW_STATUS"
    )

    assert review["value"]["review_stars"] == 1
    assert "not pathogenicity strength" in (
        review["details"]["meaning"]
    )


def test_not_found_does_not_mean_benign():
    data = interpret_clinical(
        make_contract(
            status="NOT_FOUND",
            records=False,
        )
    ).to_dict()

    assert data["status"] == "UNAVAILABLE"

    observation = data["observations"][0]

    assert observation["code"] == "CLINVAR_NOT_FOUND"

    assert "benign" in (
        observation["details"]["does_not_mean"]
    )


def test_error_is_preserved():
    dimension = interpret_clinical(
        make_contract(
            status="ERROR",
            records=False,
            error="Database unavailable",
        )
    )

    assert dimension.status == "ERROR"

    assert "database unavailable" in (
        str(dimension.to_dict()).lower()
    )


def test_found_without_records_is_partial():
    dimension = interpret_clinical(
        make_contract(
            status="FOUND",
            records=False,
        )
    )

    assert dimension.status == "PARTIAL"


def test_pathogenic_source_value_is_not_platform_classification():
    data = interpret_clinical(
        make_contract(
            significance=("Pathogenic",),
        )
    ).to_dict()

    significance = next(
        observation
        for observation in data["observations"]
        if observation["code"]
        == "CLINVAR_SUBMITTED_SIGNIFICANCE"
    )

    assert significance["value"] == ["Pathogenic"]

    assert (
        significance["details"][
            "platform_classification"
        ]
        is False
    )


def test_interpreter_has_no_acmg_or_ranking_output():
    data = interpret_clinical(
        make_contract()
    ).to_dict()

    serialized = str(data).lower()

    forbidden = (
        "'ranking_score':",
        "'tier':",
        "'pathogenicity_probability':",
        "'acmg_classification':",
    )

    for term in forbidden:
        assert term not in serialized


def test_condition_is_preserved_without_patient_causality():
    data = interpret_clinical(
        make_contract()
    ).to_dict()

    condition = next(
        observation
        for observation in data["observations"]
        if observation["code"]
        == "CLINVAR_CONDITIONS"
    )

    assert (
        condition["value"][0]["name"]
        == "Retinitis_pigmentosa"
    )

    assert "does not independently establish causality" in (
        condition["details"]["meaning"]
    )
