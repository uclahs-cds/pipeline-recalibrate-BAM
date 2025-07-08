/*
    Nextflow module for splitting input genomic intervals by chromosome

    input:
        intervals: path to set of target intervals to split
        reference_fasta: path to reference genome fasta file
        reference_fasta_index: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta

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
               pattern: "*-contig.interval_list",
               enabled: params.save_intermediate_files

    input:
    path intervals
    path reference_fasta
    path reference_fasta_index
    path reference_fasta_dict

    output:
    path "*-contig.interval_list", emit: interval_list

    script:
    """
    set -euo pipefail
    if ${params.parallelize_by_chromosome}
    then
        assembled_chr_to_exclude=''
        for i in `grep -E '^(chr|)([0-9]+|X|Y|M)\$' ${intervals}`
        do
            gatk SplitIntervals \
                -R ${reference_fasta} \
                -L ${intervals} \
                -L \$i \
                --interval-set-rule INTERSECTION \
                --scatter-count 1 \
                -O ./

            mv 0000-scattered.interval_list \$i-contig.interval_list
            assembled_chr_to_exclude="\$assembled_chr_to_exclude -XL \$i"
        done

        gatk SplitIntervals \
            -R ${reference_fasta} \
            -L ${intervals} \
            \$assembled_chr_to_exclude \
            --interval-set-rule INTERSECTION \
            --scatter-count 1 \
            -O ./

        mv 0000-scattered.interval_list noncanonical-contig.interval_list
    else
        gatk SplitIntervals \
            -R ${reference_fasta} \
            -L ${intervals} \
            --scatter-count ${params.scatter_count} \
            ${params.split_intervals_extra_args} \
            -O ./

        for interval in `ls *-scattered.interval_list`
        do
            mv \$interval `echo \$interval | cut -d '-' -f 1`-contig.interval_list
        done
    fi
    """
}
