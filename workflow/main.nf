nextflow.enable.dsl = 2


include { FASTQC }       from './modules/fastqc'
include { RAW_QC }       from './modules/raw_qc'
include { QC_GATE }      from './modules/qc_gate'
include { ALIGNMENT }    from './modules/alignment'
include { ALIGNMENT_QC } from './modules/alignment_qc'


workflow {

    /*
     * ---------------------------------------------------------
     * Samples
     * ---------------------------------------------------------
     */

    raw_samples_ch = Channel
        .fromPath(
            params.input,
            checkIfExists: true
        )
        .splitCsv(
            header: true
        )
        .map { row ->

            def r1 = file(
                row.fastq_1,
                checkIfExists: true
            )

            def r2 = file(
                row.fastq_2,
                checkIfExists: true
            )

            tuple(
                row.sample_id,
                [r1, r2]
            )
        }


    /*
     * Create a second deterministic sample channel for
     * downstream join with the QC gate.
     */

    alignment_samples_ch = Channel
        .fromPath(
            params.input,
            checkIfExists: true
        )
        .splitCsv(
            header: true
        )
        .map { row ->

            def r1 = file(
                row.fastq_1,
                checkIfExists: true
            )

            def r2 = file(
                row.fastq_2,
                checkIfExists: true
            )

            tuple(
                row.sample_id,
                [r1, r2]
            )
        }


    /*
     * ---------------------------------------------------------
     * Static resources
     * ---------------------------------------------------------
     */

    raw_qc_thresholds_ch = Channel.value(
        file(
            params.raw_qc_thresholds,
            checkIfExists: true
        )
    )


    raw_qc_script_ch = Channel.value(
        file(
            "${baseDir}/../src/genomics_platform/qc/fastqc_parser.py",
            checkIfExists: true
        )
    )


    gate_script_ch = Channel.value(
        file(
            "${baseDir}/../src/genomics_platform/qc/qc_gate.py",
            checkIfExists: true
        )
    )


    alignment_qc_thresholds_ch = Channel.value(
        file(
            params.alignment_qc_thresholds,
            checkIfExists: true
        )
    )


    alignment_qc_script_ch = Channel.value(
        file(
            "${baseDir}/../src/genomics_platform/qc/alignment_qc.py",
            checkIfExists: true
        )
    )


    reference_fasta_ch = Channel.value(
        file(
            params.reference,
            checkIfExists: true
        )
    )


    bwa_index_ch = Channel.value(
        [
            file(
                "${params.reference}.amb",
                checkIfExists: true
            ),
            file(
                "${params.reference}.ann",
                checkIfExists: true
            ),
            file(
                "${params.reference}.bwt",
                checkIfExists: true
            ),
            file(
                "${params.reference}.pac",
                checkIfExists: true
            ),
            file(
                "${params.reference}.sa",
                checkIfExists: true
            )
        ]
    )


    /*
     * ---------------------------------------------------------
     * Gate 1 — Raw FASTQ QC
     * ---------------------------------------------------------
     */

    FASTQC(
        raw_samples_ch
    )


    RAW_QC(
        FASTQC.out.zip,
        raw_qc_thresholds_ch,
        raw_qc_script_ch
    )


    QC_GATE(
        RAW_QC.out.json,
        gate_script_ch
    )


    /*
     * ---------------------------------------------------------
     * Alignment input only becomes available after QC gate.
     * ---------------------------------------------------------
     */

    alignment_ready_ch =
        alignment_samples_ch
            .join(
                QC_GATE.out.ok
            )
            .map {
                sample_id,
                reads,
                gate_token
                ->

                tuple(
                    sample_id,
                    reads
                )
            }


    /*
     * ---------------------------------------------------------
     * Alignment
     * ---------------------------------------------------------
     */

    ALIGNMENT(
        alignment_ready_ch,
        reference_fasta_ch,
        bwa_index_ch
    )


    /*
     * ---------------------------------------------------------
     * Gate 2 — Alignment QC
     * ---------------------------------------------------------
     */

    ALIGNMENT_QC(
        ALIGNMENT.out.bam,
        alignment_qc_thresholds_ch,
        alignment_qc_script_ch
    )
}
