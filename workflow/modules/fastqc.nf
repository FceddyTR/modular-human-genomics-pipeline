process FASTQC {

    tag "${sample_id}"

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    cpus 4

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.html"), emit: html
    tuple val(sample_id), path("*_fastqc.zip"), emit: zip

    script:
    """
    fastqc \
        --threads ${task.cpus} \
        ${reads}
    """
}
