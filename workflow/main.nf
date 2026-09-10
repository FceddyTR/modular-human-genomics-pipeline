nextflow.enable.dsl=2

include { FASTQC } from './modules/fastqc'

params.input  = null
params.outdir = 'results'

workflow {

    if (!params.input) {
        error "Please provide an input samplesheet using --input"
    }

    samples_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample_id,
                file(row.fastq_1),
                file(row.fastq_2)
            )
        }

    reads_ch = samples_ch.map { sample_id, fastq_1, fastq_2 ->
        tuple(sample_id, [fastq_1, fastq_2])
    }

    FASTQC(reads_ch)
}
