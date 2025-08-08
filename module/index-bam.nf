/*
    Nextflow module for indexing BAM files

    input:
        (smaple_id, bam): Tuple of sample_id and path to BAM

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_samtools: string
*/
process run_index_SAMtools {
    container params.docker_image_samtools
    publishDir path: "${params.output_dir_base}/output",
        mode: "copy",
        pattern: "*.bai"

    ext log_dir_suffix: { "-${sample_id}" }

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("${bam}.bai"), emit: indexed_out

    script:
    """
    set -euo pipefail

    samtools index ${bam}
    """
}
