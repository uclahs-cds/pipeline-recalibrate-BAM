include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
include {
    remove_intermediate_files as remove_unmerged_BAMs
    remove_intermediate_files as remove_dups_BAM
} from '../external/pipeline-Nextflow-module/modules/common/intermediate_file_removal/main.nf' addParams(
    options: [
        save_intermediate_files: params.save_intermediate_files,
        output_dir: params.output_dir_base,
        log_output_dir: "${params.log_output_dir}/process-log"
    ]
)
include { calculate_sha512 } from './checksum.nf'

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

    // Add nonassembled as last if not already present
    if (!chromosomes.contains('nonassembled')) {
        chromosomes.add('nonassembled')
    }

    return chromosomes
}

/*
    Nextflow module for gathering BAM files using GatherBamFiles (for chromosome-based parallelization)

    input:
        (sample_id, bams): tuple of sample ID and list of BAMs to gather

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
        pattern: "${output_file_name}.bam*"

    ext log_dir_suffix: { "-${sample_id}" }

    input:
    tuple val(sample_id), path(bams), val(interval_ids)

    output:
    tuple val(sample_id), path("${output_file_name}.bam"), path("${output_file_name}.bam.bai"), emit: gathered_bam
    path(bams), emit: bams_for_deletion

    script:
    // Get chromosome order outside the script context
    def dict_order = getChromosomeOrderFromDict(params.reference_fasta_dict)

    // Create map of interval to BAM file
    def bamMap = [:]
    for (int i=0; i<bams.size(); i++) {
        bamMap[interval_ids[i]] = bams[i]
    }

    // Sort BAMs according to dictionary order
    def sorted_bams = []
    dict_order.each { chrId ->
        if (bamMap.containsKey(chrId)) {
            sorted_bams.add(bamMap[chrId])
        }
    }

    // Validate that all input BAMs are being included
    if (sorted_bams.size() != bams.size()) {
        def missing_intervals = []
        for (int i = 0; i < interval_ids.size(); i++) {
            if (!dict_order.contains(interval_ids[i])) {
                missing_intervals.add(interval_ids[i])
            }
        }
        println "WARNING: ${bams.size() - sorted_bams.size()} BAM files excluded from gathering for sample ${sample_id}"
        println "WARNING: Excluded interval IDs: ${missing_intervals}"
        println "WARNING: This may result in data loss. Check reference dictionary and interval naming."
    }

    def bam_list = sorted_bams.collect{ "-INPUT '$it'" }.join(' ')

    output_file_name = generate_standard_filename(
        params.aligner,
        params.dataset_id,
        sample_id,
        [
            'additional_tools': ["GATK-${params.gatk_version}"]
        ]
    )

    """
    set -euo pipefail
    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=${workDir} \
        -jar /usr/local/share/picard-slim-2.26.10-0/picard.jar GatherBamFiles \
        ${bam_list} \
        -OUTPUT ${output_file_name}.bam \
        --CREATE_INDEX true

    mv ${output_file_name}.bai ${output_file_name}.bam.bai
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
        enabled: params.save_intermediate_files,
        pattern: "${output_file_name}"

    publishDir path: "${params.output_dir_base}/output",
        mode: "copy",
        enabled: params.save_intermediate_files,
        pattern: "${output_file_name}"

    ext log_dir_suffix: { "-${sample_id}" }

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
    bams_to_merge  // Format: [id, bam, interval_id]

    main:
    bams_to_merge
        .groupTuple()
        .map { id, bams, interval_ids ->
            [id, bams, interval_ids]
        }
        .set{ input_ch_merge_struct }

    if (params.parallelize_by_chromosome) {

        run_GatherBamFiles_Picard(input_ch_merge_struct)

        run_GatherBamFiles_Picard.out.gathered_bam
            .map{
                [
                    'sample_id': it[0],
                    'bam': it[1],
                    'bam_index': it[2]
                ]
            }
            .set{ output_ch_merge_bams }

        // No deduplication needed for chromosome-based parallelization

    } else {

        run_MergeSamFiles_Picard(input_ch_merge_struct.map{ merge_struct -> [merge_struct[0], merge_struct[1]] })

        remove_unmerged_BAMs(
            run_MergeSamFiles_Picard.out.bams_for_deletion.flatten(),
            "merge_complete"
        )

        deduplicate_records_SAMtools(run_MergeSamFiles_Picard.out.merged_bam)

        remove_dups_BAM(
            deduplicate_records_SAMtools.out.bam_for_deletion.flatten(),
            "deduplication_complete"
        )

        input_ch_index = deduplicate_records_SAMtools.out.dedup_bam

        run_index_SAMtools(input_ch_index)

        run_index_SAMtools.out.indexed_out
            .map{
                [
                    'sample_id': it[0],
                    'bam': it[1],
                    'bam_index': it[2]
                ]
            }
            .set{ output_ch_merge_bams }
    }

    output_ch_merge_bams
        .map{ [it.bam, it.bam_index] }
        .flatten()
        .set{ input_ch_calculate_sha512 }

    calculate_sha512(input_ch_calculate_sha512)

    emit:
    output_ch_merge_bams = output_ch_merge_bams
}
