from genomics_platform.evidence.gnomad_record import (
    DatasetFrequency,
    GnomADRecord,
    PopulationFrequency,
    safe_af,
)


def test_safe_af():
    assert safe_af(10, 100) == 0.1
    assert safe_af(0, 100) == 0.0
    assert safe_af(None, 100) is None
    assert safe_af(1, None) is None
    assert safe_af(1, 0) is None


def test_record_semantics():
    pop = PopulationFrequency(
        id="example",
        ac=1,
        an=100,
        homozygote_count=0,
    )

    exome = DatasetFrequency(
        ac=2,
        an=200,
        homozygote_count=0,
        populations=(pop,),
    )

    record = GnomADRecord(
        assembly="GRCh38",
        chrom="1",
        pos=12345,
        ref="A",
        alt="G",
        lookup_status="FOUND",
        exome=exome,
    )

    data = record.to_dict()

    assert data["variant_key"] == "1:12345:A:G"
    assert data["gnomad_variant_id"] == "1-12345-A-G"
    assert data["exome"]["af"] == 0.01
    assert data["exome"]["populations"][0]["af"] == 0.01
    assert data["genome"] is None


def test_not_found_is_not_zero_frequency():
    record = GnomADRecord(
        assembly="GRCh38",
        chrom="1",
        pos=999,
        ref="A",
        alt="T",
        lookup_status="NOT_FOUND",
    )

    data = record.to_dict()

    assert data["lookup_status"] == "NOT_FOUND"
    assert data["exome"] is None
    assert data["genome"] is None
