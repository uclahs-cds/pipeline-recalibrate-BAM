include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
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
    publishDir path: "${params.output_dir_base}/QC/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: '*.table'

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${sample_id}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(bundle_contest_hapmap_3p3_vcf_gz)
    path(bundle_contest_hapmap_3p3_vcf_gz_tbi)
    path(all_intervals)
    tuple val(sample_id), path(bam), path(bam_index)

    output:
    path(".command.*")
    tuple val(sample_id), path("*getpileupsummaries.table"), emit: pileupsummaries

    script:
    interval_options = all_intervals.collect{ "--intervals \"\$(realpath ${it})\"" }.join(' ')
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
        GetPileupSummaries \
        --input ${bam} \
        --reference ${reference_fasta} \
        --variant ${bundle_contest_hapmap_3p3_vcf_gz} \
        ${interval_options} \
        --output ${output_filename}.table
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

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${sample_id}/log${file(it).getName()}" }

    input:
    tuple val(normal_id), path(normal_pileup), val(tumor_id), path(tumor_pileup), val(sample_type)

    output:
    path(".command.*")
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

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${sample_id}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(all_intervals)
    tuple val(sample_id), path(bam), path(bam_index)

    output:
    path(".command.*")
    path("*_DOC*")

    when:
    params.is_DOC_run

    script:
    interval_options = all_intervals.collect{ "--intervals '\$(realpath ${it})'" }.join(' ')
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
