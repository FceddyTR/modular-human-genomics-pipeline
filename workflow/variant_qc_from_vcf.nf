nextflow.enable.dsl=2

include {
    VARIANT_QC
} from './modules/variant_qc'


params.vcf = null
params.vcf_index = null
params.sample_id = null

params.outdir = 'results'

params.variant_qc_profile = 'smoke'
params.pipeline_version = 'development'


workflow {

    if (!params.vcf) {
        error "Missing --vcf"
    }

    if (!params.vcf_index) {
        error "Missing --vcf_index"
    }

    if (!params.sample_id) {
        error "Missing --sample_id"
    }

    vcf_ch = Channel.of(
        tuple(
            params.sample_id,
            file(params.vcf),
            file(params.vcf_index)
        )
    )

    VARIANT_QC(vcf_ch)
}
