nextflow.enable.dsl=2

include {
    VARIANT_NORMALIZATION
} from './modules/variant_normalization'


params.vcf = null
params.vcf_index = null
params.sample_id = null

params.reference =
    '../reference/grch38/GRCh38_no_alt_analysis_set.fasta'

params.reference_fai =
    '../reference/grch38/GRCh38_no_alt_analysis_set.fasta.fai'

params.outdir = '../results'


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

    reference_fasta_ch = Channel.value(
        file(params.reference)
    )

    reference_fai_ch = Channel.value(
        file(params.reference_fai)
    )

    VARIANT_NORMALIZATION(
        vcf_ch,
        reference_fasta_ch,
        reference_fai_ch
    )
}
