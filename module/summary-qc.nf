include { generate_standard_filename; sanitize_string } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

/*
    Nextflow module for getting pileup summaries of BAMs

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        bundle_contest_hapmap_3p3_vcf_gz: path to contamination estimate variants
        bundle_contest_hapmap_3p3_vcf_gz_tbi: path to index of contamination estimate VCFs
        all_intervals: path to set of full target intervals
        (sample_id, bam, bam_index): tuple of sample name, BAM path, and BAM index

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
      when: params.save_intermediate_files,
      pattern: '*.table'

    tag "${sample_id}-${interval_id}"

    ext log_dir_suffix: { "-${sample_id}-${interval_id}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(bundle_contest_hapmap_3p3_vcf_gz)
    path(bundle_contest_hapmap_3p3_vcf_gz_tbi)
    path(all_intervals)
    tuple val(sample_id), path(bam), path(bam_index), val(interval_id), path(split_interval)

    output:
    tuple val(sample_id), path("*getpileupsummaries.table"), emit: pileupsummaries

    script:
    base_interval_option = "--intervals ${split_interval}"

    combined_interval_options = (params.is_targeted) ? "--intervals \"\$(realpath ${all_intervals})\" --interval-set-rule INTERSECTION" : ""

    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        sample_id,
        [
            'additional_information': "${interval_id}-getpileupsummaries"
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

/*
    Nextflow module for gathering pileup summaries per sample

    input:
        reference_fasta_dict: path to dictionary for reference fasta
        (sample_id, pileupsummaries): tuple of sample name and set of pileupsummaries

    params:
        params.output_dir_base: string(path)
        params.docker_image_gatk: string
        params.gatk_command_mem_diff: float(memory)
*/
process run_GatherPileupSummaries_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir_base}/QC/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: '*.table'

    tag "${sample_id}"

    ext log_dir_suffix: { "-${sample_id}" }

    input:
    path(reference_fasta_dict)
    tuple val(sample_id), path(pileupsummaries)

    output:
    tuple val(sample_id), path("*getpileupsummaries.table"), emit: gatheredsummaries

    script:
    input_summaries = pileupsummaries.collect{ "--I ${it}" }.join(' ')

    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        sample_id,
        [
            'additional_information': "getpileupsummaries"
        ]
    )
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
        GatherPileupSummaries \
        ${input_summaries} \
        --O ${output_filename}.table \
        --sequence-dictionary ${reference_fasta_dict}
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

    ext log_dir_suffix: { "-${sample_id}" }

    input:
    tuple val(normal_id), path(normal_pileup), val(tumor_id), path(tumor_pileup), val(sample_type)

    output:
    path("*-tumor-segmentation.table")
    path("*_alone.table"), emit: contamination
    path("*_with-matched-normal.table"), emit: tumor_normal_matched_contamination optional true

    script:
    sample_id = (sample_type == 'normal') ? normal_id : tumor_id
    calc_matched = (sample_type == 'tumor_paired')
    pileupsummaries = (sample_type == 'normal') ? normal_pileup : tumor_pileup
    single_output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        sample_id,
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
        all_intervals: path to set of full target intervals
        tuple(sample_id, bam, bam_index): tuple of sample ID, path to BAM and BAM index

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

    ext log_dir_suffix: { "-${sample_id}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(all_intervals)
    tuple val(sample_id), path(bam), path(bam_index)

    output:
    path("*_DOC*")

    when:
    params.is_DOC_run

    script:
    interval_options = all_intervals.collect{ "--intervals \"\$(realpath ${it})\"" }.join(' ')
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

workflow contamination_qc {
    take:
    bams_for_qc
    samples_with_index

    main:
    bams_for_qc
        .map{ input_bam ->
            [
                input_bam.sample_id,
                input_bam.bam,
                input_bam.bam_index,
                input_bam.interval_id,
                input_bam.split_interval
            ]
        }
        .set{ input_ch_pileupsummaries }

    input_intervals = (params.is_targeted) ? params.intervals : "${params.work_dir}/NO_FILE.interval_list"

    run_GetPileupSummaries_GATK(
        params.reference_fasta,
        params.reference_fasta_fai,
        params.reference_fasta_dict,
        params.bundle_contest_hapmap_3p3_vcf_gz,
        params.bundle_contest_hapmap_3p3_vcf_gz_tbi,
        input_intervals,
        input_ch_pileupsummaries
    )

    run_GetPileupSummaries_GATK.out.pileupsummaries
        .filter{ summary -> (new File(summary[1].toString())).length() != 0 } // Filter out non-overlapping summary files
        .groupTuple()
        .set{ input_ch_gathersummaries }

    run_GatherPileupSummaries_GATK(
        params.reference_fasta_dict,
        input_ch_gathersummaries
    )

    run_GatherPileupSummaries_GATK.out.gatheredsummaries
        .map{ gatheredpileup -> gatheredpileup[0] }
        .set{ pileupsgathered }

    samples_with_index
        .filter{ it.sample_type == 'normal' }
        .map{ it -> [sanitize_string(it.sample_id)] }
        .join(run_GatherPileupSummaries_GATK.out.gatheredsummaries)
        .set{ normal_pileupsummaries }

    samples_with_index
        .filter{ it.sample_type == 'tumor' }
        .map{ it -> [sanitize_string(it.sample_id)] }
        .join(run_GatherPileupSummaries_GATK.out.gatheredsummaries)
        .set{ tumor_pileupsummaries }

    normal_pileupsummaries.combine(tumor_pileupsummaries)
        .map{ it -> it.flatten() + ['tumor_paired'] }
        .set{ paired_pileups }

    normal_pileupsummaries.map{ it -> it + ['NO_ID', params.work_dir + '/NO_FILE.table', 'normal'] }
        .set{ normal_pileups }

    tumor_pileupsummaries.map{ ['NO_ID', params.work_dir + '/NO_FILE.table'] + it + 'tumor' }
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
    pileupsgathered = pileupsgathered
}
