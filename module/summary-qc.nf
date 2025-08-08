include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
include { sanitize_string } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
include { remove_intermediate_files } from '../external/pipeline-Nextflow-module/modules/common/intermediate_file_removal/main.nf' addParams(options: [save_intermediate_files: params.save_intermediate_files, output_dir: params.output_dir_base, log_output_dir: "${params.log_output_dir}/process-log"])
/*
    Nextflow module for getting pileup summaries of BAMs

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        bundle_contest_hapmap_3p3_vcf_gz: path to contamination estimate variants
        bundle_contest_hapmap_3p3_vcf_gz_tbi: path to index of contamination estimate VCFs
        split_intervals: list or tuple of paths to all split intervals (targeted or chromosomes)
        tuple (id, interval_bam, interval_id, interval_path): tuple of sample name, interval BAM path, interval ID, and interval path

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.docker_image_gatk: string
        params.gatk_command_mem_diff: float(memory)
*/
process run_GetPileupSummaries_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: '*.table'

    ext log_dir_suffix: { "-${id}-${interval_id}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(bundle_contest_hapmap_3p3_vcf_gz)
    path(bundle_contest_hapmap_3p3_vcf_gz_tbi)
    path(input_intervals)
    tuple val(sample_id), path(bam), path(bam_index), val(interval_id), path(split_interval), val(has_unmapped)

    output:
    tuple val(sample_id), val(interval_id), path("*getpileupsummaries*.table"), emit: interval_pileupsummaries

    script:
    base_interval_option = "--intervals ${split_interval}"
    interval_padding = params.is_targeted ? "--interval-padding 100" : ""
    combined_interval_options = params.is_targeted ?
        "--intervals \"\$(realpath ${input_intervals})\" ${interval_padding} --interval-set-rule INTERSECTION" :
        ""
    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        sample_id,
        [
            'additional_information': "getpileupsummaries-${interval_id}"
        ]
    )
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
        GetPileupSummaries \
        --input ${bam} \
        --reference ${reference_fasta} \
        --variant ${bundle_contest_hapmap_3p3_vcf_gz} \
        ${base_interval_option} \
        ${combined_interval_options} \
        --output ${output_filename}.table || touch ${output_filename}.table
    """
}

process concatenate_PileupSummaries {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir_base}/QC/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: '*.table'

    ext log_dir_suffix: { "-${id}" }

    input:
    tuple val(id), path(pileup_tables)

    output:
    tuple val(id), path("*getpileupsummaries.table"), emit: concatenated_pileupsummaries

    script:
    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        id,
        [
            'additional_information': "getpileupsummaries"
        ]
    )
    """
    set -euo pipefail

    # Sort the input files by interval ID for consistent ordering
    sorted_files=\$(ls ${pileup_tables.join(' ')} | sort -V)

    # Get the first file to extract the header
    first_file=\$(echo \$sorted_files | cut -d' ' -f1)

    # Write header from first file
    head -n 2 \$first_file > ${output_filename}.table

    # Concatenate all files, skipping headers (first 2 lines) for all files
    for file in \$sorted_files; do
        tail -n +3 \$file >> ${output_filename}.table
    done
    """
}

/*
    Nextflow module for calculating contamination

    input:
        tuple(normal_id, normal_pileup, tumor_id, tumor_pileup, sample_type):
            ID and paths for normal samples and pileup summary files

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.docker_image_gatk: string
        params.gatk_command_mem_diff: float(memory)
*/
process run_CalculateContamination_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir_base}/QC/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: '*.table'

    ext log_dir_suffix: { "-${id}" }

    input:
    tuple val(normal_id), path(normal_pileup), val(tumor_id), path(tumor_pileup), val(sample_type)

    output:
    path("*-tumor-segmentation.table")
    path("*_alone.table"), emit: contamination
    path("*_with-matched-normal.table"), emit: tumor_normal_matched_contamination optional true

    script:
    id = (sample_type == 'normal') ? normal_id : tumor_id
    calc_matched = (sample_type == 'tumor_paired')
    pileupsummaries = (sample_type == 'normal') ? normal_pileup : tumor_pileup
    single_output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        id,
        [:]
    )
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
        CalculateContamination \
        --input ${pileupsummaries} \
        --output ${single_output_filename}_alone.table \
        --tumor-segmentation ${single_output_filename}_alone-tumor-segmentation.table

    if ${calc_matched}
    then
      gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
          CalculateContamination \
          --input ${pileupsummaries} \
          --matched-normal ${normal_pileup} \
          --output ${single_output_filename}_with-matched-normal.table \
          --tumor-segmentation ${single_output_filename}_with-matched-normal-tumor-segmentation.table
    fi
    """
}

/*
    Nextflow module for calculating depth of coverage

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        targeted_intervals: path to set of full target intervals
        tuple(id, bam, bam_index): tuple of sample ID, path to BAM and BAM index

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.docker_image_gatk: string
        params.gatk_command_mem_diff: float(memory)
*/
process run_DepthOfCoverage_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir_base}/QC/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: '*_DOC*'

    ext log_dir_suffix: { "-${id}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(input_intervals)
    tuple val(sample_id), path(bam), path(bam_index)

    output:
    path("*_DOC*")

    script:
    interval_options = input_intervals.collect{ "--intervals \"\$(realpath ${it})\"" }.join(' ')
    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        sample_id,
        [
            'additional_information': "DOC"
        ]
    )
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
        DepthOfCoverage \
        --input ${bam} \
        --output ${output_filename} \
        --output-format TABLE \
        --reference ${reference_fasta} \
        --omit-depth-output-at-each-base \
        --omit-interval-statistics \
        --omit-locus-table \
        --partition-type sample \
        --partition-type readgroup \
        --partition-type library \
        ${interval_options}
    """
}

workflow summary_qc {
    take:
    bams_for_qc

    main:
    bams_for_qc.map{ [it.sample_id, it.bam, it.bam_index, it.interval_id, it.split_interval, it.has_unmapped] }
        .set{ input_ch_pileupsummaries }

    input_input_intervals = (params.is_targeted) ? params.intervals : "/scratch/NO_FILE.interval_list"

    run_GetPileupSummaries_GATK(
        params.reference_fasta,
        params.reference_fasta_fai,
        params.reference_fasta_dict,
        params.bundle_contest_hapmap_3p3_vcf_gz,
        params.bundle_contest_hapmap_3p3_vcf_gz_tbi,
        input_input_intervals,
        input_ch_pileupsummaries
    )

    // Group interval results by sample and concatenate
    run_GetPileupSummaries_GATK.out.interval_pileupsummaries
        .map{ sample_id, interval_id, pileup_file -> [sample_id, pileup_file] }
        .groupTuple()
        .set{ grouped_interval_pileups }

    concatenate_PileupSummaries(grouped_interval_pileups)

    def sample_type_map = [:]
    params.samples_to_process.each { sample_data ->
        sample_type_map["${sample_data.id}"] = sample_data.sample_type
    }

    // Each entry: [sample_id, pileup_path, sample_type]
    concatenate_PileupSummaries.out.concatenated_pileupsummaries
        .map { tup ->
            [ tup[0], tup[1], sample_type_map[tup[0]] ]
        }
        .set{ pileup_summaries_output }

    // Split into normal and tumor lists retaining id and pileup path
    pileup_summaries_output
        .filter{ it[2] == 'normal' }
        .map{ it -> [sanitize_string(it[0]), it[1]] }
        .set{ normal_pileupsummaries }

    pileup_summaries_output
        .filter{ it[2] == 'tumor' }
        .map{ it -> [sanitize_string(it[0]), it[1]] }
        .set{ tumor_pileupsummaries }

    normal_pileupsummaries.combine(tumor_pileupsummaries)
        .map{ it -> it.flatten() + ['tumor_paired'] }
        .set{ paired_pileups }

    normal_pileupsummaries.map{ it -> it + ['NO_ID', '/scratch/NO_FILE.table', 'normal'] }
        .set{ normal_pileups }

    tumor_pileupsummaries.map{ ['NO_ID', '/scratch/NO_FILE.table'] + it + 'tumor' }
        .set{ tumor_pileups }

    input_ch_calculate_contamination = normal_pileups

    def sample_types = params.samples_to_process.collect{ sample -> sample.sample_type }

    if (sample_types.contains('normal') && sample_types.contains('tumor')) {
        input_ch_calculate_contamination
            .mix(paired_pileups)
            .set{ input_ch_calculate_contamination }
    } else {
        input_ch_calculate_contamination
            .mix(tumor_pileups)
            .set{ input_ch_calculate_contamination }
    }

    run_CalculateContamination_GATK(input_ch_calculate_contamination)

    emit:
    contamination = run_CalculateContamination_GATK.out.contamination
    matched_contamination = run_CalculateContamination_GATK.out.tumor_normal_matched_contamination
}
