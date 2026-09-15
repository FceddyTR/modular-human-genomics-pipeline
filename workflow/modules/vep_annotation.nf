process VEP_ANNOTATION {

    tag "${sample_id}"

    cpus 2
    memory '4 GB'

    container 'ensemblorg/ensembl-vep:release_116.1'

    /*
     * Run the container with the host user's UID/GID so VEP can write
     * into the Nextflow work directory.
     *
     * The large VEP cache is mounted read-only instead of being staged
     * into each Nextflow work directory.
     */
    containerOptions "-u ${params.docker_uid}:${params.docker_gid} -v ${params.vep_cache_dir}:/vep_cache:ro"

    input:
    tuple val(sample_id), path(vcf), path(tbi)

    output:
    tuple val(sample_id),
          path("${sample_id}.vep.vcf.gz"),
          path("${sample_id}.vep.vcf.gz.tbi"),
          path("${sample_id}.vep.summary.html"),
          path("${sample_id}.vep.log"),
          emit: annotated_vcf

    script:
    """
    set -euo pipefail

    vep \
      --input_file ${vcf} \
      --output_file STDOUT \
      --format vcf \
      --vcf \
      --compress_output bgzip \
      --cache \
      --offline \
      --species homo_sapiens \
      --assembly GRCh38 \
      --dir_cache /vep_cache \
      --cache_version 116 \
      --symbol \
      --canonical \
      --protein \
      --numbers \
      --biotype \
      --variant_class \
      --total_length \
      --fork ${task.cpus} \
      --stats_file ${sample_id}.vep.summary.html \
      --force_overwrite \
      2> ${sample_id}.vep.log \
      > ${sample_id}.vep.vcf.gz

    tabix -f -p vcf ${sample_id}.vep.vcf.gz
    """
}
