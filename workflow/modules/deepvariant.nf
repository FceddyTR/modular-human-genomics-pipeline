process DEEPVARIANT {

    tag "${sample_id}"

    publishDir "${params.outdir}/variant_calling",
        mode: 'copy'

    container "${params.deepvariant_image}"

    cpus 2
    memory '6 GB'

    input:
    tuple val(sample_id),
          path(bam),
          path(bai)

    path reference_fasta
    path reference_fai

    output:
    tuple val(sample_id),
          path("${sample_id}.deepvariant.vcf.gz"),
          path("${sample_id}.deepvariant.vcf.gz.tbi"),
          emit: vcf

    tuple val(sample_id),
          path("${sample_id}.deepvariant.g.vcf.gz"),
          path("${sample_id}.deepvariant.g.vcf.gz.tbi"),
          emit: gvcf

    path "${sample_id}.deepvariant.log",
         emit: log

    script:
    def region_arg = params.deepvariant_region
        ? "--regions=${params.deepvariant_region}"
        : ""

    """
    set -euo pipefail

    echo "Sample: ${sample_id}" > ${sample_id}.deepvariant.log
    echo "DeepVariant image: ${params.deepvariant_image}" >> ${sample_id}.deepvariant.log
    echo "Model type: ${params.deepvariant_model_type}" >> ${sample_id}.deepvariant.log
    echo "Region: ${params.deepvariant_region ?: 'whole genome'}" >> ${sample_id}.deepvariant.log
    echo "Shards: ${task.cpus}" >> ${sample_id}.deepvariant.log
    echo "---" >> ${sample_id}.deepvariant.log

    /opt/deepvariant/bin/run_deepvariant \
        --model_type=${params.deepvariant_model_type} \
        --ref=${reference_fasta} \
        --reads=${bam} \
        --output_vcf=${sample_id}.deepvariant.vcf.gz \
        --output_gvcf=${sample_id}.deepvariant.g.vcf.gz \
        ${region_arg} \
        --num_shards=${task.cpus} \
        --intermediate_results_dir=${sample_id}.deepvariant_intermediate \
        >> ${sample_id}.deepvariant.log 2>&1

    test -s ${sample_id}.deepvariant.vcf.gz
    test -s ${sample_id}.deepvariant.vcf.gz.tbi
    test -s ${sample_id}.deepvariant.g.vcf.gz
    test -s ${sample_id}.deepvariant.g.vcf.gz.tbi

    echo "---" >> ${sample_id}.deepvariant.log
    echo "DeepVariant completed successfully." >> ${sample_id}.deepvariant.log
    """
}
