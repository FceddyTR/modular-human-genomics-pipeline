process QC_GATE {

    tag "${sample_id}"

    cpus 1
    memory '256 MB'

    input:
    tuple val(sample_id), path(qc_json)
    path gate_script

    output:
    tuple val(sample_id), path("${sample_id}.raw_qc.ok"), emit: ok

    script:
    """
    python ${gate_script} \
        --sample-id "${sample_id}" \
        --qc-json ${qc_json} \
        --output ${sample_id}.raw_qc.ok
    """
}
