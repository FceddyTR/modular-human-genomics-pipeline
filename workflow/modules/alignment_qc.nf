process ALIGNMENT_QC {

    tag "${sample_id}"

    publishDir "${params.outdir}/alignment_qc",
        mode: 'copy'

    cpus 2
    memory '2 GB'

    input:
    tuple val(sample_id),
          path(bam),
          path(bai)

    path thresholds
    path qc_script

    output:
    tuple val(sample_id),
          path("${sample_id}.alignment_qc.json"),
          emit: json

    tuple val(sample_id),
          path("${sample_id}.alignment_qc.html"),
          emit: html

    path "${sample_id}.alignment_qc.log",
         emit: log

    script:
    """
    set -euo pipefail

    samtools flagstat \
        -@ ${task.cpus} \
        ${bam} \
        > ${sample_id}.flagstat.txt

    samtools stats \
        -@ ${task.cpus} \
        ${bam} \
        > ${sample_id}.stats.txt

    samtools idxstats \
        ${bam} \
        > ${sample_id}.idxstats.txt

    mosdepth \
        -n \
        --threads ${task.cpus} \
        ${sample_id} \
        ${bam}

    python ${qc_script} \
        --sample-id "${sample_id}" \
        --flagstat ${sample_id}.flagstat.txt \
        --stats ${sample_id}.stats.txt \
        --idxstats ${sample_id}.idxstats.txt \
        --mosdepth-summary ${sample_id}.mosdepth.summary.txt \
        --mosdepth-global-dist ${sample_id}.mosdepth.global.dist.txt \
        --thresholds ${thresholds} \
        --profile "${params.alignment_qc_profile}" \
        --output-json ${sample_id}.alignment_qc.json \
        --output-html ${sample_id}.alignment_qc.html \
        --log ${sample_id}.alignment_qc.log
    """
}
