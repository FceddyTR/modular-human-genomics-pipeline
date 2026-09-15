#!/usr/bin/env bash

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

REF="${ROOT}/reference/grch38/GRCh38_no_alt_analysis_set.fasta"

FIXTURE_DIR="${ROOT}/tests/fixtures/variant_normalization"

INPUT_VCF="${FIXTURE_DIR}/input.vcf"
INPUT_VCF_GZ="${FIXTURE_DIR}/input.vcf.gz"

NORMALIZED_VCF="${FIXTURE_DIR}/observed.normalized.vcf.gz"

QC_JSON="${FIXTURE_DIR}/observed.variant_qc.json"

SAMPLE="SYNTHETIC_NORMALIZATION_TEST"


echo
echo "=============================================="
echo " Variant Normalization Regression Test"
echo "=============================================="
echo


if [[ ! -f "${REF}" ]]; then
    echo "ERROR: reference FASTA not found:"
    echo "${REF}"
    exit 1
fi


if [[ ! -f "${REF}.fai" ]]; then
    echo "ERROR: reference FASTA index not found:"
    echo "${REF}.fai"
    exit 1
fi


command -v samtools >/dev/null 2>&1 || {
    echo "ERROR: samtools not found"
    exit 1
}

command -v bcftools >/dev/null 2>&1 || {
    echo "ERROR: bcftools not found"
    exit 1
}


mkdir -p "${FIXTURE_DIR}"


get_base() {

    local pos="$1"

    samtools faidx \
        "${REF}" \
        "chr20:${pos}-${pos}" \
        | grep -v '^>' \
        | tr -d '\n' \
        | tr '[:lower:]' '[:upper:]'
}


get_sequence() {

    local start="$1"
    local end="$2"

    samtools faidx \
        "${REF}" \
        "chr20:${start}-${end}" \
        | grep -v '^>' \
        | tr -d '\n' \
        | tr '[:lower:]' '[:upper:]'
}


transition_alt() {

    case "$1" in

        A)
            echo "G"
            ;;

        G)
            echo "A"
            ;;

        C)
            echo "T"
            ;;

        T)
            echo "C"
            ;;

        *)
            echo ""
            ;;

    esac
}


transversion_alt() {

    case "$1" in

        A)
            echo "C"
            ;;

        G)
            echo "T"
            ;;

        C)
            echo "A"
            ;;

        T)
            echo "G"
            ;;

        *)
            echo ""
            ;;

    esac
}


echo "[1/8] Reading GRCh38 reference bases..."


POS_SNP=10000000
POS_INS=10000010
POS_DEL=10000020
POS_MULTI=10000030


REF_SNP="$(get_base "${POS_SNP}")"

REF_INS="$(get_base "${POS_INS}")"

REF_DEL="$(get_sequence \
    "${POS_DEL}" \
    "$((POS_DEL + 1))"
)"

REF_MULTI="$(get_base "${POS_MULTI}")"


for BASE in \
    "${REF_SNP}" \
    "${REF_INS}" \
    "${REF_MULTI}"
do

    if [[ ! "${BASE}" =~ ^[ACGT]$ ]]; then

        echo "ERROR: unexpected reference base:"
        echo "${BASE}"

        exit 1

    fi

done


if [[ ! "${REF_DEL}" =~ ^[ACGT]{2}$ ]]; then

    echo "ERROR: unexpected deletion reference sequence:"
    echo "${REF_DEL}"

    exit 1

fi


ALT_SNP="$(transition_alt "${REF_SNP}")"

ALT_INS="${REF_INS}A"

if [[ "${REF_INS}" == "A" ]]; then
    ALT_INS="${REF_INS}C"
fi

ALT_DEL="${REF_DEL:0:1}"

ALT_MULTI_1="$(transition_alt "${REF_MULTI}")"

ALT_MULTI_2="$(transversion_alt "${REF_MULTI}")"


echo "SNP:"
echo "  chr20:${POS_SNP} ${REF_SNP}>${ALT_SNP}"

echo "Insertion:"
echo "  chr20:${POS_INS} ${REF_INS}>${ALT_INS}"

echo "Deletion:"
echo "  chr20:${POS_DEL} ${REF_DEL}>${ALT_DEL}"

echo "Multiallelic:"
echo "  chr20:${POS_MULTI} ${REF_MULTI}>${ALT_MULTI_1},${ALT_MULTI_2}"

echo


echo "[2/8] Creating synthetic VCF..."


cat > "${INPUT_VCF}" <<EOF
##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=chr20,length=64444167>
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	${SAMPLE}
chr20	${POS_SNP}	test_snp	${REF_SNP}	${ALT_SNP}	60	PASS	.	GT:DP:GQ	0/1:35:99
chr20	${POS_INS}	test_insertion	${REF_INS}	${ALT_INS}	70	PASS	.	GT:DP:GQ	1/1:42:90
chr20	${POS_DEL}	test_deletion	${REF_DEL}	${ALT_DEL}	80	PASS	.	GT:DP:GQ	0/1:38:95
chr20	${POS_MULTI}	test_multiallelic	${REF_MULTI}	${ALT_MULTI_1},${ALT_MULTI_2}	90	PASS	.	GT:DP:GQ	1/2:50:99
EOF


echo "[3/8] Validating original VCF..."


bcftools view \
    "${INPUT_VCF}" \
    >/dev/null


echo "[4/8] Compressing/indexing input VCF..."


rm -f \
    "${INPUT_VCF_GZ}" \
    "${INPUT_VCF_GZ}.tbi"


bgzip \
    -c "${INPUT_VCF}" \
    > "${INPUT_VCF_GZ}"


bcftools index \
    --tbi \
    "${INPUT_VCF_GZ}"


echo "[5/8] Normalizing VCF..."


rm -f \
    "${NORMALIZED_VCF}" \
    "${NORMALIZED_VCF}.tbi"


bcftools norm \
    --fasta-ref "${REF}" \
    --multiallelics -any \
    --output-type z \
    --output "${NORMALIZED_VCF}" \
    "${INPUT_VCF_GZ}"


bcftools index \
    --tbi \
    "${NORMALIZED_VCF}"


echo "[6/8] Checking normalized records..."


RECORDS="$(
    bcftools view \
        -H \
        "${NORMALIZED_VCF}" \
        | wc -l
)"


SNP_COUNT="$(
    bcftools view \
        --types snps \
        -H \
        "${NORMALIZED_VCF}" \
        | wc -l
)"


INDEL_COUNT="$(
    bcftools view \
        --types indels \
        -H \
        "${NORMALIZED_VCF}" \
        | wc -l
)"


echo "Normalized records : ${RECORDS}"
echo "SNP records        : ${SNP_COUNT}"
echo "INDEL records      : ${INDEL_COUNT}"


if [[ "${RECORDS}" -ne 5 ]]; then

    echo
    echo "FAIL: expected 5 normalized records"
    echo "Observed: ${RECORDS}"

    exit 1

fi


if [[ "${SNP_COUNT}" -ne 3 ]]; then

    echo
    echo "FAIL: expected 3 SNP records"
    echo "Observed: ${SNP_COUNT}"

    exit 1

fi


if [[ "${INDEL_COUNT}" -ne 2 ]]; then

    echo
    echo "FAIL: expected 2 INDEL records"
    echo "Observed: ${INDEL_COUNT}"

    exit 1

fi


echo "[7/8] Running Variant QC engine..."


python \
    "${ROOT}/src/genomics_platform/qc/variant_qc.py" \
    --vcf "${NORMALIZED_VCF}" \
    --sample-id "${SAMPLE}" \
    --profile smoke \
    --pipeline-version regression-test \
    --output "${QC_JSON}"


echo "[8/8] Validating Variant QC JSON..."


python - \
    "${QC_JSON}" <<'PY'

import json
import math
import sys


path = sys.argv[1]


with open(
    path,
    "r",
    encoding="utf-8",
) as handle:

    data = json.load(handle)


metrics = data["metrics"]


expected = {
    "records": 5,
    "variant_alleles": 5,
    "transitions": 2,
    "transversions": 1,
}


errors = []


for key, value in expected.items():

    observed = metrics.get(key)

    if observed != value:

        errors.append(
            f"{key}: expected {value}, observed {observed}"
        )


types = metrics.get(
    "variant_types",
    {},
)


if types.get("SNP") != 3:

    errors.append(
        f"SNP count: expected 3, observed {types.get('SNP')}"
    )


if types.get("INDEL") != 2:

    errors.append(
        f"INDEL count: expected 2, observed {types.get('INDEL')}"
    )


titv = metrics.get(
    "ti_tv"
)


if titv is None or not math.isclose(
    titv,
    2.0,
    rel_tol=1e-9,
):

    errors.append(
        f"Ti/Tv: expected 2.0, observed {titv}"
    )


genotypes = metrics.get(
    "genotypes",
    {},
)


if genotypes.get("HET") != 4:

    errors.append(
        "HET count: expected 4, "
        f"observed {genotypes.get('HET')}"
    )


if genotypes.get("HOM_ALT") != 1:

    errors.append(
        "HOM_ALT count: expected 1, "
        f"observed {genotypes.get('HOM_ALT')}"
    )


if data.get("overall_status") != "PASS":

    errors.append(
        "overall_status: expected PASS, "
        f"observed {data.get('overall_status')}"
    )


if errors:

    print()
    print("VARIANT REGRESSION TEST FAILED")

    for error in errors:
        print(" -", error)

    sys.exit(1)


print()
print("Variant QC assertions:")
print("  records       = 5")
print("  SNPs          = 3")
print("  INDELs        = 2")
print("  transitions   = 2")
print("  transversions = 1")
print("  Ti/Tv         = 2.0")
print("  HET           = 4")
print("  HOM_ALT       = 1")

PY


echo
echo "=============================================="
echo " PASS: VARIANT NORMALIZATION + QC"
echo "=============================================="
echo

echo "Normalized VCF:"
echo "${NORMALIZED_VCF}"

echo
echo "Variant QC JSON:"
echo "${QC_JSON}"

echo
