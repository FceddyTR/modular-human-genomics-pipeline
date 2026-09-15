from genomics_platform.evidence.variant_evidence_contract import (
    ClinicalEvidence,
    FunctionalEvidence,
    PopulationEvidence,
    SourceAvailability,
    TranscriptEvidence,
    VariantEvidenceContract,
    VariantIdentity,
)
from genomics_platform.interpretation.functional_interpreter import (
    interpret_functional,
)


def make_contract(
    *,
    status="FOUND",
    payload=True,
    consequence="stop_gained",
    impact="HIGH",
    error=None,
):
    functional = FunctionalEvidence()

    if payload:
        functional = FunctionalEvidence(
            genes=("TEST1",),
            gene_ids=("ENSG000001",),
            consequences=(consequence,),
            impacts=(impact,),
            transcripts=(
                TranscriptEvidence(
                    gene_symbol="TEST1",
                    gene_id="ENSG000001",
                    transcript_id="ENST000001",
                    protein_id="ENSP000001",
                    feature_type="Transcript",
                    biotype="protein_coding",
                    consequences=(consequence,),
                    impact=impact,
                    canonical=True,
                    variant_class="SNV",
                ),
            ),
        )

    return VariantEvidenceContract(
        identity=VariantIdentity(
            assembly="GRCh38",
            chrom="1",
            pos=100,
            ref="A",
            alt="T",
        ),
        functional=functional,
        clinical=ClinicalEvidence(),
        population=PopulationEvidence(),
        availability=(
            SourceAvailability(
                source="VEP",
                status=status,
                error=error,
            ),
            SourceAvailability(
                source="ClinVar",
                status="NOT_FOUND",
            ),
            SourceAvailability(
                source="gnomAD",
                status="NOT_FOUND",
            ),
        ),
    )


def test_found_functional_evidence_is_available():
    dimension = interpret_functional(
        make_contract()
    )

    assert dimension.status == "AVAILABLE"

    codes = {
        observation.code
        for observation in dimension.observations
    }

    assert "VEP_GENE_CONTEXT" in codes
    assert "VEP_CONSEQUENCES" in codes
    assert "VEP_IMPACTS" in codes
    assert "VEP_TRANSCRIPT" in codes


def test_consequence_is_preserved():
    data = interpret_functional(
        make_contract(
            consequence="stop_gained",
        )
    ).to_dict()

    observation = next(
        item
        for item in data["observations"]
        if item["code"] == "VEP_CONSEQUENCES"
    )

    assert observation["value"] == [
        "stop_gained"
    ]


def test_vep_high_impact_is_not_clinical_strength():
    data = interpret_functional(
        make_contract(
            impact="HIGH",
        )
    ).to_dict()

    observation = next(
        item
        for item in data["observations"]
        if item["code"] == "VEP_IMPACTS"
    )

    assert observation["value"] == ["HIGH"]

    assert "not this platform's clinical" in (
        observation["details"]["does_not_mean"]
    )


def test_transcript_context_is_preserved():
    data = interpret_functional(
        make_contract()
    ).to_dict()

    transcript = next(
        item
        for item in data["observations"]
        if item["code"] == "VEP_TRANSCRIPT"
    )

    assert (
        transcript["value"]["transcript_id"]
        == "ENST000001"
    )

    assert transcript["value"]["canonical"] is True

    assert (
        transcript["details"][
            "clinical_interpretation"
        ]
        is False
    )


def test_not_found_does_not_mean_functionally_neutral():
    data = interpret_functional(
        make_contract(
            status="NOT_FOUND",
            payload=False,
        )
    ).to_dict()

    assert data["status"] == "UNAVAILABLE"

    observation = data["observations"][0]

    assert observation["code"] == "VEP_NOT_FOUND"

    assert "functionally_neutral" in (
        observation["details"]["does_not_mean"]
    )


def test_error_is_preserved():
    dimension = interpret_functional(
        make_contract(
            status="ERROR",
            payload=False,
            error="VEP unavailable",
        )
    )

    assert dimension.status == "ERROR"

    assert "vep unavailable" in (
        str(dimension.to_dict()).lower()
    )


def test_found_without_payload_is_partial():
    dimension = interpret_functional(
        make_contract(
            status="FOUND",
            payload=False,
        )
    )

    assert dimension.status == "PARTIAL"
    assert dimension.observations == ()


def test_stop_gained_does_not_create_pathogenicity():
    data = interpret_functional(
        make_contract(
            consequence="stop_gained",
            impact="HIGH",
        )
    ).to_dict()

    serialized = str(data).lower()

    forbidden = (
        "'pathogenic': true",
        "'pathogenicity_probability':",
        "'ranking_score':",
        "'tier':",
        "'acmg_classification':",
    )

    for term in forbidden:
        assert term not in serialized


def test_gene_context_is_preserved():
    data = interpret_functional(
        make_contract()
    ).to_dict()

    gene = next(
        item
        for item in data["observations"]
        if item["code"] == "VEP_GENE_CONTEXT"
    )

    assert gene["value"]["genes"] == ["TEST1"]
    assert gene["value"]["gene_ids"] == [
        "ENSG000001"
    ]
