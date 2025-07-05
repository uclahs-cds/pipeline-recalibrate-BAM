include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
include {
    remove_intermediate_files as remove_unmerged_BAMs
    remove_intermediate_files as remove_merged_BAM
    } from '../external/pipeline-Nextflow-module/modules/common/intermediate_file_removal/main.nf' addParams(
        options: [
            save_intermediate_files: params.save_intermediate_files,
            output_dir: params.output_dir_base,
            log_output_dir: "${params.log_output_dir}/process-log"
            ]
        )
include { calculate_sha512 } from './checksum.nf'

/*
    Function to create chromosome ordering from reference dictionary
    Input: Path to reference dictionary file
    Output: List of chromosome names in reference order
*/
def getChromosomeOrderFromDict(dictPath) {
    def chromosomes = []
    def dictFile = file(dictPath)
    
    if (dictFile.exists()) {
        dictFile.eachLine { line ->
            if (line.startsWith('@SQ')) {
                def matcher = line =~ /SN:(\S+)/
                if (matcher.find()) {
                    chromosomes.add(matcher.group(1))
                }
            }
        }
    }
    
    // Add noncanonical as last if not already present
    if (!chromosomes.contains('noncanonical')) {
        chromosomes.add('noncanonical')
    }
    
    return chromosomes
}

/*
    Function to sort BAM files in natural chromosome order for GatherBamFiles
    Input: List of tuples containing [bam_file, interval_id], reference dictionary path
    Output: List of BAM files sorted in natural chromosome order
*/
def sortBamsInNaturalOrder(bams, dictPath) {
    // Get chromosome order from reference dictionary
    def chromosomeOrder = getChromosomeOrderFromDict(dictPath)
    
    def bamMap = [:]
    bams.each { bam ->
        def intervalId
        if (bam instanceof List && bam.size() == 2) {
            // New format: [bam_file, interval_id]
            intervalId = bam[1]
            bamMap[intervalId] = bam[0]
        } else {
            // Legacy format: extract from filename (fallback)
            intervalId = bam.name.find(/_(recalibrated|indelrealigned)-([^.]+)\.bam$/) { match, type, id -> id }
            if (intervalId) {
                bamMap[intervalId] = bam
            }
        }
    }

    def sortedBams = []
    chromosomeOrder.each { chrId ->
        if (bamMap.containsKey(chrId)) {
            sortedBams.add(bamMap[chrId])
        }
    }

    return sortedBams
}

/*
    Nextflow module for gathering BAM files using GatherBamFiles (for chromosome-based parallelization)

    input:
        (sample_id, bams_with_intervals): tuple of sample ID and list of [bam, interval_id] tuples

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.docker_image_picard: string
        params.gatk_command_mem_diff: float(memory)
        params.reference_fasta_dict: string(path)
*/
process run_GatherBamFiles_Picard {
    container params.docker_image_picard
    publishDir path: "${params.output_dir_base}/output",
        mode: "copy",
        pattern: "${output_file_name}"

    ext log_dir_suffix: { "-${sample_id}" }

    when:
    params.parallelize_by_chromosome

    input:
    tuple val(sample_id), val(bams_with_intervals)

    output:
    tuple val(sample_id), path(output_file_name), emit: gathered_bam
    path("*.bam"), emit: output_ch_deletion

    script:
    // Sort BAMs in natural chromosome order using reference dictionary
    sorted_bams = sortBamsInNaturalOrder(bams_with_intervals, params.reference_fasta_dict)
    bam_list = sorted_bams.collect{ "-INPUT '$it'" }.join(' ')
    output_file_name = generate_standard_filename(
        params.aligner,
        params.dataset_id,
        sample_id,
        [
            'additional_tools': ["GATK-${params.gatk_version}"]
        ]
    )
    output_file_name = "${output_file_name}.bam"
    
    // Create symlinks for BAM files to make them available in work directory
    bam_symlinks = bams_with_intervals.collect { bam_info ->
        def bam_file = bam_info instanceof List ? bam_info[0] : bam_info
        return "ln -sf '${bam_file}' ./"
    }.join(' && ')
    
    """
    set -euo pipefail
    ${bam_symlinks}
    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=${workDir} \
        -jar /usr/local/share/picard-slim-2.26.10-0/picard.jar GatherBamFiles \
        ${bam_list} \
        -OUTPUT ${output_file_name}
    """
}

/*
    Nextflow module for merging BAM files (for scatter-count based parallelization)

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
        enabled: false,  // Never publish directly; always go through deduplication when not using chromosome parallelization
        pattern: "${output_file_name}"

    ext log_dir_suffix: { "-${sample_id}" }

    when:
    !params.parallelize_by_chromosome

    input:
    tuple val(sample_id), path(bams)

    output:
    tuple val(sample_id), path(output_file_name), emit: merged_bam
    path(bams), emit: output_ch_deletion

    script:
    all_bams = bams.collect{ "-INPUT '$it'" }.toSorted().join(' ')
    output_file_name = generate_standard_filename(
        params.aligner,
        params.dataset_id,
        sample_id,
        [
            'additional_tools': ["GATK-${params.gatk_version}"],
            'additional_information': "realigned_recalibrated_merged"
        ]
    )
    output_file_name = "${output_file_name}.bam"
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
        pattern: "${output_file_name}"

    ext log_dir_suffix: { "-${sample_id}" }

    when:
    !params.parallelize_by_chromosome

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(output_file_name), emit: dedup_bam
    path(bam), emit: bam_for_deletion

    script:
    output_file_name = generate_standard_filename(
        params.aligner,
        params.dataset_id,
        sample_id,
        [
            'additional_tools': ["GATK-${params.gatk_version}"]
        ]
    )
    output_file_name = "${output_file_name}.bam"
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

workflow merge_bams {
    take:
    bams_to_merge  // Format: [sample, bam, interval_id] or [sample, bam] (legacy)

    main:
    // Handle both new format with interval_id and legacy format
    bams_to_merge
        .map{ input_data ->
            if (input_data.containsKey('interval_id')) {
                // New format: group by sample and include interval_id
                [input_data.sample, [input_data.bam, input_data.interval_id]]
            } else {
                // Legacy format: group by sample only
                [input_data.sample, input_data.bam]
            }
        }
        .groupTuple()
        .set{ input_ch_merge }

    if (params.parallelize_by_chromosome) {
        // Use GatherBamFiles for chromosome-based parallelization
        run_GatherBamFiles_Picard(input_ch_merge)

        // Don't delete immediately - return for coordinated deletion
        run_GatherBamFiles_Picard.out.output_ch_deletion
            .flatten()
            .set{ interval_bams_for_deletion }

        // No deduplication needed for chromosome-based parallelization
        input_ch_index = run_GatherBamFiles_Picard.out.gathered_bam

    } else {
        // Use MergeSamFiles for scatter-count based parallelization
        // Extract just the BAM files for merging (legacy compatibility)
        input_ch_merge
            .map{ sample_id, bam_data ->
                def bams = bam_data.collect { item ->
                    item instanceof List ? item[0] : item
                }
                [sample_id, bams]
            }
            .set{ input_ch_merge_files_only }

        run_MergeSamFiles_Picard(input_ch_merge_files_only)

        // Don't delete immediately - return for coordinated deletion
        run_MergeSamFiles_Picard.out.output_ch_deletion
            .flatten()
            .set{ interval_bams_for_deletion }

        deduplicate_records_SAMtools(run_MergeSamFiles_Picard.out.merged_bam)

        remove_merged_BAM(
            deduplicate_records_SAMtools.out.bam_for_deletion.flatten(),
            "deduplication_complete"
        )

        input_ch_index = deduplicate_records_SAMtools.out.dedup_bam
    }

    run_index_SAMtools(input_ch_index)

    run_index_SAMtools.out.indexed_out
        .map{
            [
                'sample': it[0],
                'bam': it[1],
                'bam_index': it[2]
            ]
        }
        .set{ output_ch_merge_bams }

    output_ch_merge_bams
        .map{ [it.bam, it.bam_index] }
        .flatten()
        .set{ input_ch_calculate_sha512 }

    calculate_sha512(input_ch_calculate_sha512)

    emit:
    output_ch_merge_bams = output_ch_merge_bams
    interval_bams_for_deletion = interval_bams_for_deletion
}
