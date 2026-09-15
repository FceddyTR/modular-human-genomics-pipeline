import pytest

from genomics_platform.evidence.variant_evidence_contract import (
    ClinicalConditionEvidence,
    ClinicalEvidence,
    ClinicalRecordEvidence,
    FunctionalEvidence,
    PopulationDatasetEvidence,
    PopulationEvidence,
    SourceAvailability,
    TranscriptEvidence,
    VariantEvidenceContract,
    VariantIdentity,
)
from genomics_platform.interpretation.interpretation_engine import (
    interpret_variant,
)


def make_contract():
    return VariantEvidenceContract(
        identity=VariantIdentity(
            assembly="GRCh38",
            chrom="1",
            pos=66926,
            ref="AG",
            alt="A",
        ),
        functional=FunctionalEvidence(
            genes=("OR4F5",),
            gene_ids=("ENSG00000186092",),
            consequences=("intron_variant",),
            impacts=("MODIFIER",),
            transcripts=(
                TranscriptEvidence(
                    gene_symbol="OR4F5",
                    gene_id="ENSG00000186092",
                    transcript_id="ENST00000641515",
                    protein_id="ENSP00000493376",
                    feature_type="Transcript",
                    biotype="protein_coding",
                    consequences=("intron_variant",),
                    impact="MODIFIER",
                    canonical=True,
                    variant_class="deletion",
                ),
            ),
        ),
        clinical=ClinicalEvidence(
            records=(
                ClinicalRecordEvidence(
                    variation_id="3385321",
                    allele_id="3544463",
                    vcv_accession="VCV003385321",
                    clinical_significance=(
                        "Uncertain_significance",
                    ),
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
            ),
        ),
        population=PopulationEvidence(
            exome=PopulationDatasetEvidence(
                ac=1,
                an=100000,
                af=0.00001,
                homozygote_count=0,
                populations=(),
            ),
            genome=None,
        ),
        availability=(
            SourceAvailability(
                source="VEP",
                status="FOUND",
            ),
            SourceAvailability(
                source="ClinVar",
                status="FOUND",
            ),
            SourceAvailability(
                source="gnomAD",
                status="FOUND",
            ),
        ),
    )


def test_engine_returns_complete_profile():
    profile = interpret_variant(
        make_contract()
    )

    assert (
        profile.variant_key
        == "1:66926:AG:A"
    )

    assert profile.population.status == "PARTIAL"
    assert profile.clinical.status == "AVAILABLE"
    assert profile.functional.status == "AVAILABLE"


def test_engine_runs_all_three_v01_dimensions():
    data = interpret_variant(
        make_contract()
    ).to_dict()

    dimensions = data["dimensions"]

    assert dimensions["population"] is not None
    assert dimensions["clinical"] is not None
    assert dimensions["functional"] is not None

    assert dimensions["phenotype"] is None
    assert dimensions["inheritance"] is None
    assert dimensions["gene_disease"] is None


def test_engine_preserves_clinvar_vus():
    data = interpret_variant(
        make_contract()
    ).to_dict()

    observations = data[
        "dimensions"
    ]["clinical"]["observations"]

    significance = next(
        observation
        for observation in observations
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


def test_engine_preserves_population_frequency():
    data = interpret_variant(
        make_contract()
    ).to_dict()

    observations = data[
        "dimensions"
    ]["population"]["observations"]

    exome = next(
        observation
        for observation in observations
        if observation["code"]
        == "GNOMAD_EXOME_SUMMARY"
    )

    assert exome["value"]["af"] == 0.00001


def test_engine_preserves_functional_annotation():
    data = interpret_variant(
        make_contract()
    ).to_dict()

    observations = data[
        "dimensions"
    ]["functional"]["observations"]

    consequence = next(
        observation
        for observation in observations
        if observation["code"]
        == "VEP_CONSEQUENCES"
    )

    assert consequence["value"] == [
        "intron_variant"
    ]


def test_engine_has_no_ranking_or_classification():
    data = interpret_variant(
        make_contract()
    ).to_dict()

    forbidden_keys = {
        "ranking_score",
        "score",
        "tier",
        "pathogenicity_probability",
        "acmg_classification",
        "diagnosis",
    }

    assert forbidden_keys.isdisjoint(
        data.keys()
    )

    for dimension in (
        "population",
        "clinical",
        "functional",
    ):
        assert forbidden_keys.isdisjoint(
            data["dimensions"][dimension].keys()
        )


def test_engine_reports_error_dimension():
    contract = make_contract()

    contract = VariantEvidenceContract(
        identity=contract.identity,
        functional=contract.functional,
        clinical=contract.clinical,
        population=contract.population,
        availability=(
            SourceAvailability(
                source="VEP",
                status="FOUND",
            ),
            SourceAvailability(
                source="ClinVar",
                status="FOUND",
            ),
            SourceAvailability(
                source="gnomAD",
                status="ERROR",
                error="Service overloaded",
            ),
        ),
    )

    profile = interpret_variant(contract)

    assert profile.population.status == "ERROR"
    assert profile.error_dimensions == [
        "population"
    ]

    assert profile.clinical.status == "AVAILABLE"
    assert profile.functional.status == "AVAILABLE"


def test_engine_rejects_wrong_input_type():
    with pytest.raises(TypeError):
        interpret_variant(
            {"variant": "not-a-contract"}
        )
