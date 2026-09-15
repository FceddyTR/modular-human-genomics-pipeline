nextflow.enable.dsl=2

include { VEP_ANNOTATION } from './modules/vep_annotation'

params.input_vcf     = null
params.sample_id     = 'SAMPLE'
params.vep_cache_dir = null
params.docker_uid    = 1000
params.docker_gid    = 1000
params.outdir        = 'results/vep'

workflow {

    if (!params.input_vcf) {
        error "Missing --input_vcf"
    }

    if (!params.vep_cache_dir) {
        error "Missing --vep_cache_dir"
    }

    input_vcf = file(params.input_vcf, checkIfExists: true)
    input_tbi = file("${params.input_vcf}.tbi", checkIfExists: true)

    input_ch = Channel.of(
        tuple(
            params.sample_id,
            input_vcf,
            input_tbi
        )
    )

    VEP_ANNOTATION(input_ch)
}
