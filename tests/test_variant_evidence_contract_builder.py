from genomics_platform.evidence.variant_evidence_bundle import (
    EvidenceComponent,
    VariantEvidenceBundle,
)
from genomics_platform.evidence.variant_evidence_contract_builder import (
    build_variant_evidence_contract,
)


def make_bundle(
    *,
    vep_status="FOUND",
    clinvar_status="FOUND",
    gnomad_status="FOUND",
):
    vep_data = [
        {
            "Allele": "-",
            "Consequence": "intron_variant",
            "IMPACT": "MODIFIER",
            "SYMBOL": "OR4F5",
            "Gene": "ENSG00000186092",
            "Feature_type": "Transcript",
            "Feature": "ENST00000641515",
            "BIOTYPE": "protein_coding",
            "ENSP": "ENSP00000493376",
            "CANONICAL": "YES",
            "VARIANT_CLASS": "deletion",
            "_normalized": {
                "allele": "-",
                "consequence": "intron_variant",
                "impact": "MODIFIER",
                "symbol": "OR4F5",
                "gene_id": "ENSG00000186092",
                "feature_type": "Transcript",
                "transcript_id": "ENST00000641515",
                "biotype": "protein_coding",
                "protein_id": "ENSP00000493376",
                "canonical": True,
                "variant_class": "deletion",
            },
        }
    ]

    clinvar_data = [
        {
            "assembly": "GRCh38",
            "chrom": "1",
            "pos": 66926,
            "ref": "AG",
            "alt": "A",
            "variation_id": "3385321",
            "allele_id": "3544463",
            "vcv_accession": "VCV003385321",
            "gene_symbols": ["OR4F5"],
            "gene_ids": ["79501"],
            "clinical_significance": [
                "Uncertain_significance"
            ],
            "conflicting_significance": [],
            "review_status": (
                "criteria_provided,_single_submitter"
            ),
            "review_stars": 1,
            "conditions": [
                {
                    "name": "Retinitis_pigmentosa",
                    "identifiers": [
                        "MONDO:MONDO:0019200",
                        "OMIM:268000",
                        "Orphanet:791",
                    ],
                }
            ],
            "hgvs": [
                "NC_000001.11:g.66927del"
            ],
            "molecular_consequences": [
                {
                    "so_id": "SO:0001627",
                    "consequence": "intron_variant",
                }
            ],
            "scv_accessions": [
                "SCV005419006"
            ],
            "variant_type": "Deletion",
            "variant_type_so": "SO:0000159",
            "origin": "0",
            "source": "ClinVar",
            "source_release": "2026-09-13",
            "source_reference": "GRCh38",
            "schema_version": "1.1",
            "variant_key": "1:66926:AG:A",
        }
    ]

    gnomad_data = {
        "schema_version": "1.0",
        "source": "gnomAD",
        "source_dataset": "gnomad_r4",
        "source_release": "4.1.1",
        "access_method": "gnomAD Browser GraphQL API",
        "assembly": "GRCh38",
        "variant_key": "1:66926:AG:A",
        "gnomad_variant_id": "1-66926-AG-A",
        "chrom": "1",
        "pos": 66926,
        "ref": "AG",
        "alt": "A",
        "lookup_status": "FOUND",
        "exome": {
            "ac": 2,
            "an": 200000,
            "af": 0.00001,
            "homozygote_count": 0,
            "populations": [
                {
                    "id": "example_exome_population",
                    "ac": 1,
                    "an": 10000,
                    "af": 0.0001,
                    "homozygote_count": 0,
                }
            ],
        },
        "genome": {
            "ac": 3,
            "an": 100000,
            "af": 0.00003,
            "homozygote_count": 0,
            "populations": [
                {
                    "id": "example_genome_population",
                    "ac": 1,
                    "an": 5000,
                    "af": 0.0002,
                    "homozygote_count": 0,
                }
            ],
        },
        "error": None,
    }

    return VariantEvidenceBundle(
        assembly="GRCh38",
        chrom="1",
        pos=66926,
        ref="AG",
        alt="A",
        vep=EvidenceComponent(
            source="VEP",
            status=vep_status,
            data=vep_data if vep_status == "FOUND" else [],
            error=(
                "VEP unavailable"
                if vep_status == "ERROR"
                else None
            ),
            provenance={
                "source_release": "116",
                "access_method": "local_cache",
                "parser_version": "0.1.0",
            },
        ),
        clinvar=EvidenceComponent(
            source="ClinVar",
            status=clinvar_status,
            data=(
                clinvar_data
                if clinvar_status == "FOUND"
                else []
            ),
            error=(
                "ClinVar unavailable"
                if clinvar_status == "ERROR"
                else None
            ),
            provenance={
                "source_release": "2026-09-13",
                "source_dataset": "GRCh38",
                "access_method": "SQLite",
            },
        ),
        gnomad=EvidenceComponent(
            source="gnomAD",
            status=gnomad_status,
            data=(
                gnomad_data
                if gnomad_status == "FOUND"
                else None
            ),
            error=(
                "gnomAD unavailable"
                if gnomad_status == "ERROR"
                else None
            ),
            provenance={
                "source_release": "4.1.1",
                "source_dataset": "gnomad_r4",
                "access_method": (
                    "gnomAD Browser GraphQL API"
                ),
            },
        ),
    )


def test_bundle_to_contract_identity():
    contract = build_variant_evidence_contract(
        make_bundle()
    )

    assert contract.identity.variant_key == "1:66926:AG:A"
    assert contract.identity.assembly == "GRCh38"
    assert contract.identity.chrom == "1"
    assert contract.identity.pos == 66926
    assert contract.identity.ref == "AG"
    assert contract.identity.alt == "A"


def test_bundle_to_contract_preserves_vep_evidence():
    data = build_variant_evidence_contract(
        make_bundle()
    ).to_dict()

    functional = data["functional"]

    assert functional["genes"] == ["OR4F5"]
    assert functional["gene_ids"] == [
        "ENSG00000186092"
    ]
    assert functional["consequences"] == [
        "intron_variant"
    ]
    assert functional["impacts"] == ["MODIFIER"]
    assert functional["canonical_transcript_count"] == 1

    transcript = functional["transcripts"][0]

    assert transcript["gene_symbol"] == "OR4F5"
    assert (
        transcript["transcript_id"]
        == "ENST00000641515"
    )
    assert transcript["canonical"] is True


def test_bundle_to_contract_preserves_clinvar_evidence():
    data = build_variant_evidence_contract(
        make_bundle()
    ).to_dict()

    clinical = data["clinical"]

    assert clinical["record_count"] == 1

    record = clinical["records"][0]

    assert record["variation_id"] == "3385321"
    assert record["allele_id"] == "3544463"
    assert record["vcv_accession"] == "VCV003385321"
    assert record["clinical_significance"] == [
        "Uncertain_significance"
    ]
    assert record["review_stars"] == 1

    assert (
        record["conditions"][0]["name"]
        == "Retinitis_pigmentosa"
    )


def test_bundle_to_contract_preserves_exome_and_genome_separately():
    data = build_variant_evidence_contract(
        make_bundle()
    ).to_dict()

    population = data["population"]

    assert population["exome"]["af"] == 0.00001
    assert population["genome"]["af"] == 0.00003

    assert (
        population["exome"]["populations"][0][
            "population_id"
        ]
        == "example_exome_population"
    )

    assert (
        population["genome"]["populations"][0][
            "population_id"
        ]
        == "example_genome_population"
    )


def test_not_found_does_not_create_fake_evidence():
    contract = build_variant_evidence_contract(
        make_bundle(
            clinvar_status="NOT_FOUND",
            gnomad_status="NOT_FOUND",
        )
    )

    data = contract.to_dict()

    assert data["clinical"]["record_count"] == 0
    assert data["population"]["exome"] is None
    assert data["population"]["genome"] is None

    assert (
        data["availability"]["ClinVar"]["status"]
        == "NOT_FOUND"
    )
    assert (
        data["availability"]["gnomAD"]["status"]
        == "NOT_FOUND"
    )

    assert data["evidence_complete"] is True


def test_error_is_preserved_and_marks_contract_incomplete():
    contract = build_variant_evidence_contract(
        make_bundle(
            gnomad_status="ERROR",
        )
    )

    data = contract.to_dict()

    assert (
        data["availability"]["gnomAD"]["status"]
        == "ERROR"
    )
    assert (
        data["availability"]["gnomAD"]["error"]
        == "gnomAD unavailable"
    )

    assert data["evidence_complete"] is False
    assert data["error_sources"] == ["gnomAD"]


def test_provenance_survives_bundle_to_contract():
    data = build_variant_evidence_contract(
        make_bundle()
    ).to_dict()

    provenance = {
        item["source"]: item
        for item in data["provenance"]
    }

    assert provenance["VEP"]["source_release"] == "116"
    assert (
        provenance["VEP"]["metadata"]["parser_version"]
        == "0.1.0"
    )

    assert (
        provenance["ClinVar"]["source_release"]
        == "2026-09-13"
    )

    assert (
        provenance["gnomAD"]["source_dataset"]
        == "gnomad_r4"
    )


def test_contract_does_not_invent_interpretation():
    data = build_variant_evidence_contract(
        make_bundle()
    ).to_dict()

    serialized = str(data).lower()

    forbidden = (
        "ranking_score",
        "pathogenicity_probability",
        "acmg_classification",
        "'tier':",
    )

    for term in forbidden:
        assert term not in serialized
