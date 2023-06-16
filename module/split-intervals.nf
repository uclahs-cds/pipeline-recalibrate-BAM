/*
    Nextflow module for splitting input genomic intervals by chromosome

    input:
        intervals: path to set of target intervals to split
        reference: path to reference genome fasta file
        reference_index: path to index for reference fasta
        reference_dict: path to dictionary for reference fasta

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
*/
process run_SplitIntervals_GATK {
    container params.docker_image_gatk

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*-scattered.interval_list",
               enabled: params.save_intermediate_files
    publishDir "${params.log_output_dir}/process-log",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
    path intervals
    path reference
    path reference_index
    path reference_dict

    output:
    path "*-scattered.interval_list", emit: interval_list
    path ".command.*"

    script:
    """
    set -euo pipefail
    if ${params.parallelize_by_chromosome}
    then
        assembled_chr_to_exclude=''
        for i in `grep -E '^(chr|)([0-9]+|X|Y|M)\$' ${intervals}`
        do
            gatk SplitIntervals \
                -R ${reference} \
                -L ${intervals} \
                -L \$i \
                --interval-set-rule INTERSECTION \
                --scatter-count 1 \
                -O ./

            mv 0000-scattered.interval_list \$i-scattered.interval_list
            assembled_chr_to_exclude="\$assembled_chr_to_exclude -XL \$i"
        done

        gatk SplitIntervals \
            -R ${reference} \
            -L ${intervals} \
            \$assembled_chr_to_exclude \
            --interval-set-rule INTERSECTION \
            --scatter-count 1 \
            -O ./

        mv 0000-scattered.interval_list nonassembled-scattered.interval_list
    else
        gatk SplitIntervals \
            -R ${reference} \
            -L ${intervals} \
            --scatter-count ${params.scatter_count} \
            ${params.split_intervals_extra_args} \
            -O ./
    fi
    """
}
