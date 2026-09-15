process VARIANT_QC {

    tag "${sample_id}"

    publishDir "${params.outdir}/variant_qc",
        mode: 'copy'

    cpus 1
    memory '1 GB'

    input:
    tuple val(sample_id),
          path(vcf),
          path(vcf_index)

    output:
    tuple val(sample_id),
          path("${sample_id}.variant_qc.json"),
          emit: json

    path "${sample_id}.variant_qc.log",
         emit: log

    script:
    """
    set -euo pipefail

    python ${projectDir}/../src/genomics_platform/qc/variant_qc.py \
        --vcf ${vcf} \
        --sample-id ${sample_id} \
        --profile ${params.variant_qc_profile} \
        --pipeline-version ${params.pipeline_version} \
        --output ${sample_id}.variant_qc.json \
        > ${sample_id}.variant_qc.log 2>&1

    test -s ${sample_id}.variant_qc.json
    """
}
