process RAW_QC {

    tag "${sample_id}"

    publishDir "${params.outdir}/raw_qc", mode: 'copy'

    cpus 1
    memory '2 GB'

    input:
    tuple val(sample_id), path(fastqc_zips)
    path thresholds
    path qc_script

    output:
    tuple val(sample_id), path("${sample_id}.raw_qc.json"), emit: json
    tuple val(sample_id), path("${sample_id}.raw_qc.html"), emit: html
    tuple val(sample_id), path("${sample_id}.raw_qc.log"), emit: log

    script:

    def zip_list = fastqc_zips instanceof List
        ? fastqc_zips
        : [fastqc_zips]

    def zip_args = zip_list
        .collect { it.toString() }
        .join(' ')

    """
    python ${qc_script} \
        --sample-id "${sample_id}" \
        --fastqc-zip ${zip_args} \
        --thresholds ${thresholds} \
        --output-json ${sample_id}.raw_qc.json \
        --output-html ${sample_id}.raw_qc.html \
        --log ${sample_id}.raw_qc.log
    """
}
