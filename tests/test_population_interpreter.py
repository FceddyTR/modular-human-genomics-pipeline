from genomics_platform.evidence.variant_evidence_contract import (
    ClinicalEvidence,
    FunctionalEvidence,
    PopulationDatasetEvidence,
    PopulationEvidence,
    PopulationObservation,
    SourceAvailability,
    VariantEvidenceContract,
    VariantIdentity,
)
from genomics_platform.interpretation.population_interpreter import (
    interpret_population,
)


def make_contract(
    *,
    status="FOUND",
    exome=True,
    genome=True,
    error=None,
):
    exome_data = (
        PopulationDatasetEvidence(
            ac=2,
            an=200000,
            af=0.00001,
            homozygote_count=0,
            populations=(
                PopulationObservation(
                    population_id="example_exome_population",
                    ac=1,
                    an=10000,
                    af=0.0001,
                    homozygote_count=0,
                ),
            ),
        )
        if exome
        else None
    )

    genome_data = (
        PopulationDatasetEvidence(
            ac=3,
            an=100000,
            af=0.00003,
            homozygote_count=0,
            populations=(
                PopulationObservation(
                    population_id="example_genome_population",
                    ac=1,
                    an=5000,
                    af=0.0002,
                    homozygote_count=0,
                ),
            ),
        )
        if genome
        else None
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
        clinical=ClinicalEvidence(),
        population=PopulationEvidence(
            exome=exome_data,
            genome=genome_data,
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
                status=status,
                error=error,
            ),
        ),
    )


def test_found_exome_and_genome_are_available():
    dimension = interpret_population(
        make_contract()
    )

    data = dimension.to_dict()

    assert data["status"] == "AVAILABLE"

    codes = {
        observation["code"]
        for observation in data["observations"]
    }

    assert "GNOMAD_EXOME_SUMMARY" in codes
    assert "GNOMAD_GENOME_SUMMARY" in codes
    assert "GNOMAD_EXOME_POPULATION" in codes
    assert "GNOMAD_GENOME_POPULATION" in codes


def test_exome_and_genome_values_are_preserved():
    data = interpret_population(
        make_contract()
    ).to_dict()

    observations = {
        observation["code"]: observation
        for observation in data["observations"]
        if observation["code"].endswith("_SUMMARY")
    }

    assert (
        observations["GNOMAD_EXOME_SUMMARY"]["value"]["af"]
        == 0.00001
    )

    assert (
        observations["GNOMAD_GENOME_SUMMARY"]["value"]["af"]
        == 0.00003
    )


def test_missing_genome_is_partial_not_error():
    dimension = interpret_population(
        make_contract(
            exome=True,
            genome=False,
        )
    )

    assert dimension.status == "PARTIAL"

    assert any(
        "genome" in limitation
        for limitation in dimension.limitations
    )


def test_not_found_does_not_mean_rare():
    dimension = interpret_population(
        make_contract(
            status="NOT_FOUND",
            exome=False,
            genome=False,
        )
    )

    data = dimension.to_dict()

    assert data["status"] == "UNAVAILABLE"

    assert (
        data["observations"][0]["code"]
        == "GNOMAD_NOT_FOUND"
    )

    does_not_mean = data["observations"][0][
        "details"
    ]["does_not_mean"]

    assert "AF=0" in does_not_mean
    assert "rare" in does_not_mean
    assert "pathogenic" in does_not_mean


def test_error_is_not_converted_to_not_found():
    dimension = interpret_population(
        make_contract(
            status="ERROR",
            exome=False,
            genome=False,
            error="Service overloaded",
        )
    )

    assert dimension.status == "ERROR"

    serialized = str(
        dimension.to_dict()
    ).lower()

    assert "service overloaded" in serialized


def test_found_without_payload_is_partial():
    dimension = interpret_population(
        make_contract(
            status="FOUND",
            exome=False,
            genome=False,
        )
    )

    assert dimension.status == "PARTIAL"
    assert dimension.observations == ()


def test_population_interpreter_does_not_assign_rarity():
    data = interpret_population(
        make_contract()
    ).to_dict()

    serialized = str(data).lower()

    forbidden = (
        "'rarity':",
        "'rare': true",
        "'rare': false",
        "pathogenicity_probability",
        "ranking_score",
        "acmg_classification",
    )

    for term in forbidden:
        assert term not in serialized


def test_population_ids_are_preserved_not_reinterpreted():
    data = interpret_population(
        make_contract()
    ).to_dict()

    population_observations = [
        observation
        for observation in data["observations"]
        if observation["code"].endswith(
            "_POPULATION"
        )
    ]

    ids = {
        observation["value"]["population_id"]
        for observation in population_observations
    }

    assert "example_exome_population" in ids
    assert "example_genome_population" in ids
