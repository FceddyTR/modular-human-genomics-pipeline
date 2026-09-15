from genomics_platform.evidence.gnomad_record import (
    DatasetFrequency,
    GnomADRecord,
)
from genomics_platform.evidence.variant_evidence_builder import (
    build_variant_evidence_bundle,
)


def make_gnomad_found():
    return GnomADRecord(
        assembly="GRCh38",
        chrom="1",
        pos=12345,
        ref="A",
        alt="G",
        lookup_status="FOUND",
        exome=DatasetFrequency(
            ac=1,
            an=1000,
            homozygote_count=0,
        ),
        genome=DatasetFrequency(
            ac=2,
            an=2000,
            homozygote_count=0,
        ),
    )


def test_complete_bundle():
    bundle = build_variant_evidence_bundle(
        assembly="GRCh38",
        chrom="chr1",
        pos=12345,
        ref="a",
        alt="g",
        vep_annotations=[
            {
                "Consequence": "missense_variant",
                "SYMBOL": "TEST1",
            }
        ],
        clinvar_records=[
            {
                "clinical_significance": ["Uncertain_significance"],
            }
        ],
        gnomad_record=make_gnomad_found(),
        vep_release="116",
    )

    data = bundle.to_dict()

    assert data["variant_key"] == "1:12345:A:G"
    assert data["evidence_complete"] is True
    assert data["error_sources"] == []

    assert data["evidence"]["functional"]["status"] == "FOUND"
    assert data["evidence"]["clinical"]["status"] == "FOUND"
    assert data["evidence"]["population"]["status"] == "FOUND"

    assert (
        data["evidence"]["population"]["data"]["exome"]["af"]
        == 0.001
    )


def test_not_found_is_completed_lookup():
    gnomad = GnomADRecord(
        assembly="GRCh38",
        chrom="1",
        pos=12345,
        ref="A",
        alt="G",
        lookup_status="NOT_FOUND",
    )

    bundle = build_variant_evidence_bundle(
        assembly="GRCh38",
        chrom="1",
        pos=12345,
        ref="A",
        alt="G",
        vep_annotations=[],
        clinvar_records=[],
        gnomad_record=gnomad,
    )

    data = bundle.to_dict()

    assert data["evidence_complete"] is True
    assert data["evidence"]["functional"]["status"] == "NOT_FOUND"
    assert data["evidence"]["clinical"]["status"] == "NOT_FOUND"
    assert data["evidence"]["population"]["status"] == "NOT_FOUND"

    assert data["evidence"]["population"]["data"]["exome"] is None
    assert data["evidence"]["population"]["data"]["genome"] is None


def test_error_is_distinct_from_not_found():
    gnomad = GnomADRecord(
        assembly="GRCh38",
        chrom="1",
        pos=12345,
        ref="A",
        alt="G",
        lookup_status="ERROR",
        error="Synthetic network failure",
    )

    bundle = build_variant_evidence_bundle(
        assembly="GRCh38",
        chrom="1",
        pos=12345,
        ref="A",
        alt="G",
        vep_annotations=[
            {"Consequence": "missense_variant"}
        ],
        clinvar_records=[],
        gnomad_record=gnomad,
    )

    data = bundle.to_dict()

    assert data["evidence_complete"] is False
    assert data["error_sources"] == ["gnomAD"]

    assert data["evidence"]["clinical"]["status"] == "NOT_FOUND"
    assert data["evidence"]["population"]["status"] == "ERROR"

    assert (
        data["evidence"]["population"]["error"]
        == "Synthetic network failure"
    )


def test_missing_source_result_becomes_error():
    bundle = build_variant_evidence_bundle(
        assembly="GRCh38",
        chrom="1",
        pos=12345,
        ref="A",
        alt="G",
        vep_annotations=None,
        clinvar_records=None,
        gnomad_record=None,
    )

    data = bundle.to_dict()

    assert data["evidence_complete"] is False
    assert set(data["error_sources"]) == {
        "VEP",
        "ClinVar",
        "gnomAD",
    }


def test_bundle_does_not_classify_variant():
    bundle = build_variant_evidence_bundle(
        assembly="GRCh38",
        chrom="1",
        pos=12345,
        ref="A",
        alt="G",
        vep_annotations=[
            {"Consequence": "stop_gained"}
        ],
        clinvar_records=[
            {"clinical_significance": ["Pathogenic"]}
        ],
        gnomad_record=make_gnomad_found(),
    )

    data = bundle.to_dict()

    assert "classification" not in data
    assert "tier" not in data
    assert "score" not in data
    assert "pathogenicity_probability" not in data
