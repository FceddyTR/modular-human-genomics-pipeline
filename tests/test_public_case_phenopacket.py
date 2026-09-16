from copy import deepcopy

from genomics_platform.ranking.benchmark.public_cases import (
    PHENOPACKET_STORE_RELEASE,
    PUBLIC_CASE_PROTOCOL_VERSION,
    evaluate_phenopacket,
)


def _packet():
    return {
        "id": "public-case-001",
        "phenotypicFeatures": [
            {
                "type": {
                    "id": "HP:0001250"
                }
            },
            {
                "type": {
                    "id": "HP:0001263"
                }
            },
            {
                "type": {
                    "id": "HP:0004322"
                }
            },
            {
                "type": {
                    "id": "HP:0002014"
                },
                "excluded": True,
            },
        ],
        "interpretations": [
            {
                "progressStatus": "SOLVED",
                "diagnosis": {
                    "genomicInterpretations": [
                        {
                            "interpretationStatus": (
                                "CAUSATIVE"
                            ),
                            "variantInterpretation": {
                                "variationDescriptor": {
                                    "allelicState": {
                                        "id": (
                                            "GENO:0000135"
                                        ),
                                        "label": (
                                            "heterozygous"
                                        ),
                                    },
                                    "geneContext": {
                                        "valueId": (
                                            "HGNC:25662"
                                        ),
                                        "symbol": "AAGAB",
                                    },
                                    "vcfRecord": {
                                        "genomeAssembly": (
                                            "hg38"
                                        ),
                                        "chrom": "chr15",
                                        "pos": 67236719,
                                        "ref": "TG",
                                        "alt": "T",
                                    },
                                }
                            },
                        }
                    ]
                },
            }
        ],
    }


def test_protocol_versions_are_frozen():
    assert PUBLIC_CASE_PROTOCOL_VERSION == "0.1"
    assert PHENOPACKET_STORE_RELEASE == "0.1.27"


def test_eligible_public_case():
    result = evaluate_phenopacket(
        _packet()
    )

    assert result.eligible is True
    assert result.exclusion_reasons == ()

    assert result.present_hpo_terms == (
        "HP:0001250",
        "HP:0001263",
        "HP:0004322",
    )

    assert result.absent_hpo_terms == (
        "HP:0002014",
    )

    assert result.causal_genes == (
        "AAGAB",
    )

    assert len(
        result.causal_variants
    ) == 1

    variant = result.causal_variants[0]

    assert variant.variant_type == "INDEL"

    assert variant.variant_key == (
        "chr15:67236719:TG:T"
    )


def test_requires_three_present_hpo_terms():
    packet = _packet()

    packet["phenotypicFeatures"] = (
        packet["phenotypicFeatures"][:2]
    )

    result = evaluate_phenopacket(
        packet
    )

    assert result.eligible is False

    assert "present_hpo_lt_3" in (
        result.exclusion_reasons
    )


def test_rejects_non_hg38_causal_variant():
    packet = _packet()

    descriptor = (
        packet["interpretations"][0]
        ["diagnosis"]
        ["genomicInterpretations"][0]
        ["variantInterpretation"]
        ["variationDescriptor"]
    )

    descriptor["vcfRecord"][
        "genomeAssembly"
    ] = "hg19"

    result = evaluate_phenopacket(
        packet
    )

    assert result.eligible is False

    assert (
        "no_usable_causative_small_variant"
        in result.exclusion_reasons
    )


def test_rejects_structural_only_case():
    packet = _packet()

    descriptor = (
        packet["interpretations"][0]
        ["diagnosis"]
        ["genomicInterpretations"][0]
        ["variantInterpretation"]
        ["variationDescriptor"]
    )

    descriptor.pop("vcfRecord")

    descriptor["structuralType"] = {
        "id": "SO:1000029",
        "label": "chromosomal_deletion",
    }

    result = evaluate_phenopacket(
        packet
    )

    assert result.eligible is False


def test_requires_explicit_causative_status():
    packet = _packet()

    genomic = (
        packet["interpretations"][0]
        ["diagnosis"]
        ["genomicInterpretations"][0]
    )

    genomic[
        "interpretationStatus"
    ] = "CONTRIBUTORY"

    result = evaluate_phenopacket(
        packet
    )

    assert result.eligible is False


def test_retains_two_causative_variants():
    packet = _packet()

    genomic = (
        packet["interpretations"][0]
        ["diagnosis"]
        ["genomicInterpretations"][0]
    )

    second = deepcopy(genomic)

    descriptor = (
        second["variantInterpretation"]
        ["variationDescriptor"]
    )

    descriptor["vcfRecord"] = {
        "genomeAssembly": "hg38",
        "chrom": "chr15",
        "pos": 67231842,
        "ref": "A",
        "alt": "ATT",
    }

    packet["interpretations"][0][
        "diagnosis"
    ][
        "genomicInterpretations"
    ].append(second)

    result = evaluate_phenopacket(
        packet
    )

    assert result.eligible is True

    assert len(
        result.causal_variants
    ) == 2


def test_noncausative_variant_does_not_become_truth():
    packet = _packet()

    genomic = (
        packet["interpretations"][0]
        ["diagnosis"]
        ["genomicInterpretations"][0]
    )

    decoy = deepcopy(genomic)

    decoy[
        "interpretationStatus"
    ] = "CONTRIBUTORY"

    packet["interpretations"][0][
        "diagnosis"
    ][
        "genomicInterpretations"
    ].append(decoy)

    result = evaluate_phenopacket(
        packet
    )

    assert len(
        result.causal_variants
    ) == 1
