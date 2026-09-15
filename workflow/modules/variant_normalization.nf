process VARIANT_NORMALIZATION {

    tag "${sample_id}"

    publishDir "${params.outdir}/variant_normalization",
        mode: 'copy'

    cpus 1
    memory '2 GB'

    input:
    tuple val(sample_id),
          path(vcf),
          path(vcf_index)

    path reference_fasta
    path reference_fai

    output:
    tuple val(sample_id),
          path("${sample_id}.normalized.vcf.gz"),
          path("${sample_id}.normalized.vcf.gz.tbi"),
          emit: vcf

    path "${sample_id}.normalization.log",
         emit: log

    script:
    """
    set -euo pipefail

    echo "Sample: ${sample_id}" > ${sample_id}.normalization.log
    echo "Input VCF: ${vcf}" >> ${sample_id}.normalization.log
    echo "Reference: ${reference_fasta}" >> ${sample_id}.normalization.log
    echo "---" >> ${sample_id}.normalization.log

    bcftools norm \
        --fasta-ref ${reference_fasta} \
        --multiallelics -any \
        --output-type z \
        --output ${sample_id}.normalized.vcf.gz \
        ${vcf} \
        >> ${sample_id}.normalization.log 2>&1

    bcftools index \
        --tbi \
        ${sample_id}.normalized.vcf.gz \
        >> ${sample_id}.normalization.log 2>&1

    bcftools view \
        --header-only \
        ${sample_id}.normalized.vcf.gz \
        > /dev/null

    test -s ${sample_id}.normalized.vcf.gz
    test -s ${sample_id}.normalized.vcf.gz.tbi

    echo "---" >> ${sample_id}.normalization.log
    echo "Variant normalization completed successfully." \
        >> ${sample_id}.normalization.log
    """
}
