process ALIGNMENT {

    tag "${sample_id}"

    publishDir "${params.outdir}/alignment", mode: 'copy'

    cpus 8
    memory '6 GB'

    input:
    tuple val(sample_id), path(reads)
    path reference_fasta
    path bwa_index

    output:
    tuple val(sample_id),
          path("${sample_id}.markdup.bam"),
          path("${sample_id}.markdup.bam.bai"),
          emit: bam

    path "${sample_id}.alignment.log",
         emit: log

    script:

    def r1 = reads[0]
    def r2 = reads[1]

    """
    set -euo pipefail

    echo "Sample: ${sample_id}" > ${sample_id}.alignment.log
    echo "Reference: ${reference_fasta}" >> ${sample_id}.alignment.log
    echo "BWA threads: 6" >> ${sample_id}.alignment.log
    echo "samtools sort threads: 2" >> ${sample_id}.alignment.log
    echo "samtools sort memory/thread: 384M" >> ${sample_id}.alignment.log
    echo "---" >> ${sample_id}.alignment.log

    bwa mem \
        -t 6 \
        -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}_PCRFREE' \
        ${reference_fasta} \
        ${r1} \
        ${r2} \
        2>> ${sample_id}.alignment.log \
    | samblaster \
        2>> ${sample_id}.alignment.log \
    | samtools sort \
        -@ 2 \
        -m 384M \
        -o ${sample_id}.markdup.bam \
        - \
        2>> ${sample_id}.alignment.log

    samtools quickcheck -v \
        ${sample_id}.markdup.bam \
        >> ${sample_id}.alignment.log 2>&1

    samtools index \
        -@ 2 \
        ${sample_id}.markdup.bam \
        2>> ${sample_id}.alignment.log

    echo "---" >> ${sample_id}.alignment.log
    echo "Alignment completed successfully." >> ${sample_id}.alignment.log
    """
}
