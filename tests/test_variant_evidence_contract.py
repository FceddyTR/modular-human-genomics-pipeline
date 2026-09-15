import pytest

from genomics_platform.evidence.variant_evidence_contract import (
    ClinicalConditionEvidence,
    ClinicalEvidence,
    ClinicalRecordEvidence,
    EvidenceProvenance,
    FunctionalEvidence,
    PopulationDatasetEvidence,
    PopulationEvidence,
    PopulationObservation,
    SourceAvailability,
    TranscriptEvidence,
    VariantEvidenceContract,
    VariantIdentity,
)


def make_contract(
    *,
    vep_status="FOUND",
    clinvar_status="FOUND",
    gnomad_status="FOUND",
):
    identity = VariantIdentity(
        assembly="GRCh38",
        chrom="1",
        pos=66926,
        ref="AG",
        alt="A",
    )

    transcript = TranscriptEvidence(
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
    )

    functional = FunctionalEvidence(
        genes=("OR4F5",),
        gene_ids=("ENSG00000186092",),
        consequences=("intron_variant",),
        impacts=("MODIFIER",),
        transcripts=(transcript,),
    )

    condition = ClinicalConditionEvidence(
        name="Retinitis_pigmentosa",
        identifiers=(
            "MONDO:MONDO:0019200",
            "OMIM:268000",
            "Orphanet:791",
        ),
    )

    clinical_record = ClinicalRecordEvidence(
        variation_id="3385321",
        allele_id="3544463",
        vcv_accession="VCV003385321",
        clinical_significance=("Uncertain_significance",),
        review_status="criteria_provided,_single_submitter",
        review_stars=1,
        conditions=(condition,),
        hgvs=("NC_000001.11:g.66927del",),
        scv_accessions=("SCV005419006",),
    )

    clinical = ClinicalEvidence(records=(clinical_record,))

    population = PopulationEvidence(
        exome=PopulationDatasetEvidence(
            ac=2,
            an=200000,
            af=0.00001,
            homozygote_count=0,
            populations=(
                PopulationObservation(
                    population_id="example_population",
                    ac=1,
                    an=10000,
                    af=0.0001,
                    homozygote_count=0,
                ),
            ),
        ),
        genome=None,
    )

    availability = (
        SourceAvailability("VEP", vep_status),
        SourceAvailability("ClinVar", clinvar_status),
        SourceAvailability("gnomAD", gnomad_status),
    )

    provenance = (
        EvidenceProvenance(
            source="VEP",
            source_release="116",
            access_method="local_cache",
        ),
        EvidenceProvenance(
            source="ClinVar",
            source_release="2026-09-13",
            source_dataset="GRCh38",
            access_method="SQLite",
        ),
        EvidenceProvenance(
            source="gnomAD",
            source_dataset="gnomad_r4",
            access_method="GraphQL API",
        ),
    )

    return VariantEvidenceContract(
        identity=identity,
        functional=functional,
        clinical=clinical,
        population=population,
        availability=availability,
        provenance=provenance,
    )


def test_contract_serializes_normalized_evidence():
    contract = make_contract()
    data = contract.to_dict()

    assert data["schema_version"] == "1.0"
    assert data["variant_key"] == "1:66926:AG:A"

    assert data["identity"]["assembly"] == "GRCh38"
    assert data["functional"]["genes"] == ["OR4F5"]
    assert data["functional"]["consequences"] == ["intron_variant"]
    assert data["functional"]["canonical_transcript_count"] == 1

    assert data["clinical"]["record_count"] == 1
    assert (
        data["clinical"]["records"][0]["variation_id"]
        == "3385321"
    )

    assert data["population"]["exome"]["af"] == 0.00001
    assert data["population"]["genome"] is None

    assert data["availability"]["VEP"]["status"] == "FOUND"
    assert data["availability"]["ClinVar"]["status"] == "FOUND"
    assert data["availability"]["gnomAD"]["status"] == "FOUND"

    assert data["evidence_complete"] is True
    assert data["error_sources"] == []


def test_not_found_is_completed_lookup():
    contract = make_contract(
        clinvar_status="NOT_FOUND",
        gnomad_status="NOT_FOUND",
    )

    assert contract.evidence_complete is True
    assert contract.error_sources == []


def test_error_is_not_completed_lookup():
    contract = make_contract(gnomad_status="ERROR")

    assert contract.evidence_complete is False
    assert contract.error_sources == ["gnomAD"]


def test_invalid_availability_status_is_rejected():
    with pytest.raises(ValueError):
        SourceAvailability(
            source="gnomAD",
            status="MISSING",
        )


def test_contract_contains_no_classification_or_ranking_output():
    data = make_contract().to_dict()

    forbidden_keys = {
        "pathogenic",
        "benign",
        "tier",
        "ranking_score",
        "pathogenicity_probability",
        "acmg_classification",
    }

    assert forbidden_keys.isdisjoint(data.keys())
    assert forbidden_keys.isdisjoint(data["functional"].keys())
    assert forbidden_keys.isdisjoint(data["clinical"].keys())
    assert forbidden_keys.isdisjoint(data["population"].keys())


def test_exome_and_genome_are_not_merged():
    data = make_contract().to_dict()

    assert "exome" in data["population"]
    assert "genome" in data["population"]
    assert data["population"]["exome"] is not None
    assert data["population"]["genome"] is None


def test_provenance_is_preserved():
    data = make_contract().to_dict()

    sources = {
        item["source"]: item
        for item in data["provenance"]
    }

    assert sources["VEP"]["source_release"] == "116"
    assert sources["ClinVar"]["source_release"] == "2026-09-13"
    assert sources["gnomAD"]["source_dataset"] == "gnomad_r4"
