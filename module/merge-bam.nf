include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
/*
    Nextflow module for merging BAM files

    input:
        (sample_id, bams): tuple of sample ID and list of BAMs to merge

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_picard: string
        params.gatk_command_mem_diff: float(memory)
        params.parallelize_by_chromosome: bool.
*/
process run_MergeSamFiles_Picard {
    container params.docker_image_picard
    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        mode: "copy",
        enabled: !params.parallelize_by_chromosome && params.save_intermediate_files,
        pattern: "${output_file_name}"

    publishDir path: "${params.output_dir_base}/output",
        mode: "copy",
        enabled: params.parallelize_by_chromosome,
        pattern: "${output_file_name}"

    publishDir path: "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}-${sample_id}/log${file(it).getName()}" }

    input:
    tuple val(sample_id), path(bams)

    output:
    path(".command.*")
    tuple val(sample_id), path(output_file_name), emit: merged_bam

    script:
    all_bams = bams.collect{ "-INPUT '$it'" }.join(' ')
    additional_information = (params.parallelize_by_chromosome) ?
        ".bam" :
        "realigned_recalibrated_merged.bam"
    output_file_name = generate_standard_filename(
        params.aligner,
        params.dataset_id,
        sample_id,
        [
            'additional_tools': ["GATK-${params.gatk_version}"],
            'additional_information': additional_information
        ]
    )
    """
    set -euo pipefail
    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=${workDir} \
        -jar /usr/local/share/picard-slim-2.26.10-0/picard.jar MergeSamFiles \
        ${all_bams} \
        -OUTPUT ${output_file_name} \
        -SORT_ORDER coordinate \
        -ASSUME_SORTED false \
        -USE_THREADING true \
        -VALIDATION_STRINGENCY LENIENT
    """
}

/*
    Nextflow module for removing duplicated identical records from interval processing.
    Due to parallelization via interval splitting, reads that overlap two intervals end up in
    the BAMs for both intervals. When merged, these records get duplicated, causing potential
    issues in downstream pipelines and analysis. Since records are sorted, uniq is able to de-deuplicate
    these records.

    See Issue #79 (https://github.com/uclahs-cds/pipeline-call-gSNP/issues/79)

    input:
        (sample_id, bam): Tuple of sample ID and BAM to deduplicate

    params:
        params.parallelize_by_chromosome: bool.
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.docker_image_samtools: string
*/
process deduplicate_records_SAMtools {
    container params.docker_image_samtools
    publishDir path: "${params.output_dir_base}/output",
        mode: "copy",
        pattern: "${output_file_name}.bam"

    publishDir path: "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}-${sample_id}/log${file(it).getName()}" }

    when:
    !params.parallelize_by_chromosome

    input:
    tuple val(sample_id), path(bam)

    output:
    path(".command.*")
    tuple val(sample_id), path(output_file_name), emit: dedup_bam
    path(bam), emit: bam_for_deletion

    script:
    output_file_name = generate_standard_filename(
        params.aligner,
        params.dataset_id,
        sample_id,
        [
            'additional_tools': ["GATK-${params.gatk_version}"],
            'additional_information': ".bam"
        ]
    )
    """
    samtools view \
        -h \
        ${bam} | \
        awk '(\$1 \$2 \$3 \$4 \$5 \$6 \$7 \$8 \$9 \$10 \$11)!=f_p && NR>1 {print f} {f=\$0} {f_p=(\$1 \$2 \$3 \$4 \$5 \$6 \$7 \$8 \$9 \$10 \$11)} END {print f}' | \
        samtools view \
        -b \
        -o ${output_file_name}
    """
}

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
    publishDir path: "${params.output_dir}/output",
        mode: "copy",
        pattern: "*.bai"

    publishDir path: "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}-${sample_id}/log${file(it).getName()}" }

    input:
    tuple val(sample_id), path(bam)

    output:
    path(".command.*")
    tuple val(sample_id), path(bam), path("${bam}.bai"), emit: indexed_out

    script:
    """
    set -euo pipefail

    samtools index ${bam}
    """
}

workflow merge_bams {
    take:
    bams_to_merge

    main:
    bams_to_merge
        .map{ [it.sample, it.bam] }
        .groupTuple()
        .set{ input_ch_merge }

    run_MergeSamFiles_Picard(input_ch_merge)

    deduplicate_records_SAMtools(run_MergeSamFiles_Picard.out.merged_bam)

    input_ch_index = (params.parallelize_by_chromosome) ?
        run_MergeSamFiles_Picard.out.merged_bam :
        deduplicate_records_SAMtools.out.dedup_bam

    run_index_SAMtools(input_ch_index)
}