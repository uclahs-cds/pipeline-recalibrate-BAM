include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
include {
    remove_intermediate_files as remove_ungathered_BAMs
    } from '../external/pipeline-Nextflow-module/modules/common/intermediate_file_removal/main.nf' addParams(
        options: [
            save_intermediate_files: params.save_intermediate_files,
            output_dir: params.output_dir_base,
            log_output_dir: "${params.log_output_dir}/process-log"
            ]
        )
include { calculate_sha512 } from './checksum.nf'

List getChromosomeOrderFromDict(dictPath) {
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
    tuple val(sample_id), path(bams)

    output:
    tuple val(sample_id), path("${output_file_name}.bam"), path("${output_file_name}.bam.bai"), emit: gathered_bam
    path(bams), emit: bams_for_deletion

    script:
    def bam_list = bams.collect{ "-INPUT '$it'" }.join(' ')

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

workflow gatherbams {
    take:
    bams_to_merge

    main:
    List chromosome_list = getChromosomeOrderFromDict(params.reference_fasta_dict)
    Map sort_key_map = [:]
    chromosome_list.eachWithIndex { item, index ->
        sort_key_map[item] = index;
    }

    bams_to_merge
        .map { interval_bam ->
            [
                interval_bam.id,
                [
                    'key': sort_key_map.getOrDefault(interval_bam.interval_id, 1000),
                    'bam': interval_bam.bam
                ]
            ]
        }
        .groupTuple()
        .map { grouped_interval_bams ->
            [
                grouped_interval_bams[0],
                grouped_interval_bams[1].sort{ bam_data -> bam_data.key }
            ]
        }
        .map{ sorted_interval_bams ->
            [
                sorted_interval_bams[0],
                sorted_interval_bams[1].collect{ element -> element.bam }
            ]
        }
        .set{ input_ch_gatherbams }

    run_GatherBamFiles_Picard(input_ch_gatherbams)

    remove_ungathered_BAMs(
        run_GatherBamFiles_Picard.out.bams_for_deletion.flatten(),
        "gather_complete"
    )

    run_GatherBamFiles_Picard.out.gathered_bam
        .map{ a_gathered_bam ->
            [
                'sample': a_gathered_bam[0],
                'bam': a_gathered_bam[1],
                'bam_index': a_gathered_bam[2]
            ]
        }
        .set{ output_ch_gathered_bams }

    output_ch_gathered_bams
        .map{ [it.bam, it.bam_index] }
        .flatten()
        .set{ input_ch_calculate_sha512 }

    calculate_sha512(input_ch_calculate_sha512)

    emit:
    gathered_bams = output_ch_gathered_bams
}
