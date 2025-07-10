include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
include {
    remove_intermediate_files
    } from '../external/pipeline-Nextflow-module/modules/common/intermediate_file_removal/main.nf' addParams(
        options: [
            save_intermediate_files: params.save_intermediate_files,
            output_dir: params.output_dir_base,
            log_output_dir: "${params.log_output_dir}/process-log"
            ]
        )
/*
    Nextflow module for generating base recalibration table

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        bundle_mills_and_1000g_gold_standard_indels_vcf_gz: path to standard Mills and 1000 genomes variants
        bundle_mills_and_1000g_gold_standard_indels_vcf_gz_tbi: path to index file for Mills and 1000g variants
        bundle_known_indels_vcf_gz: path to set of known indels
        bundle_known_indels_vcf_gz_tbi: path to index of known indels VCF
        bundle_v0_dbsnp138_vcf_gz: path to dbSNP variants
        bundle_v0_dbsnp138_vcf_gz_tbi: path to index of dbSNP variants
        all_intervals: list or tuple of paths to all split intervals
        input_bams: list or tuple of paths to input BAMs (raw or indel realigned)
        input_bams_bai: list or tuple of paths to input BAM indices
        sample_id: identifier for the sample

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
        params.is_targeted: bool. Indicator of whether in targeted exome mode or in WGS mode
        params.gatk_command_mem_diff: float(memory)
        params.run_indel_realignment: bool. Indicator of whether to run indel realignment
*/
process run_BaseRecalibrator_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir_base}/QC/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: "*.grp"

    ext log_dir_suffix: { "-${sample_id}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(bundle_mills_and_1000g_gold_standard_indels_vcf_gz)
    path(bundle_mills_and_1000g_gold_standard_indels_vcf_gz_tbi)
    path(bundle_known_indels_vcf_gz)
    path(bundle_known_indels_vcf_gz_tbi)
    path(bundle_v0_dbsnp138_vcf_gz)
    path(bundle_v0_dbsnp138_vcf_gz_tbi)
    path(intervals)
    path(recal_tables)
    tuple path(input_bams), path(input_bams_bai), val(sample_id), val(interval_id)

    output:
    tuple val(sample_id), path("${sample_id}_${interval_id}_recalibration_table.grp"), emit: recalibration_table

    script:
    all_input_bams = input_bams.collect{ "--input '${it}'" }.join(' ')
    targeted_options = params.is_targeted ? "--interval-padding 100" : ""
    """
    set -euo pipefail
    if [ ! -f ${sample_id}_${interval_id}_recalibration_table.grp ]
    then
        gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
            BaseRecalibrator \
            ${all_input_bams} \
            --reference ${reference_fasta} \
            --verbosity INFO \
            --known-sites ${bundle_mills_and_1000g_gold_standard_indels_vcf_gz} \
            --known-sites ${bundle_known_indels_vcf_gz} \
            --known-sites ${bundle_v0_dbsnp138_vcf_gz} \
            --output ${sample_id}_${interval_id}_recalibration_table.grp \
            --intervals ${interval} \
            ${targeted_options} \
            --read-filter SampleReadFilter \
            --sample ${sample_id}
    fi
    """
}

/*
    Nextflow module for gathering BQSR reports from multiple intervals into a single report per sample

    input:
        sample_id: identifier for the sample
        recalibration_tables: list of paths to interval-specific recalibration tables for this sample

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.docker_image_gatk: string
        params.gatk_command_mem_diff: float(memory)
*/
process run_GatherBQSRReports_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir_base}/QC/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: "*.grp"

    ext log_dir_suffix: { "-${sample_id}" }

    input:
    tuple val(sample_id), path(recalibration_tables)

    output:
    tuple val(sample_id), path("${sample_id}_recalibration_table.grp"), emit: gathered_recalibration_table

    script:
    input_args = recalibration_tables.collect{ "--input ${it}" }.join(' ')
    """
    set -euo pipefail

    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
        GatherBQSRReports \
        ${input_args} \
        --output ${sample_id}_recalibration_table.grp
    """
}

/*
    Nextflow module for recalibrating base quality scores in BAM file

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        input_bam: path to input BAM (raw or indel realigned)
        input_bam_index: path to input BAM index
        interval: path to specific intervals file associated with input BAM
        includes_unmapped: boolean to indicate if unmapped reads are included in input BAM
        sample_id: Identifier for sample being processed
        recalibration_table: path to gathered base recalibration table

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
        params.gatk_command_mem_diff: float(memory)
        params.is_emit_original_quals: bool. Indicator of whether to keep original base quality scores
*/
process run_ApplyBQSR_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "*_recalibrated-*"

    ext log_dir_suffix: { "-${interval_id}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    tuple path(input_bam),
          path(input_bam_index),
          val(interval_id),
          path(interval),
          val(includes_unmapped),
          val(sample_id),
          path(recalibration_table)

    output:
    path("*_recalibrated-*"), emit: output_ch_apply_bqsr
    tuple path(input_bam), path(input_bam_index), emit: output_ch_deletion

    script:
    unmapped_interval_option = (includes_unmapped) ? "--intervals unmapped" : ""
    combined_interval_options = "--intervals ${interval} ${unmapped_interval_option}"
    all_commands = sample_id.collect{
        """
        gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \\
            ApplyBQSR \\
            --input ${input_bam} \\
            --bqsr-recal-file ${it}_recalibration_table.grp \\
            --reference ${reference_fasta} \\
            --read-filter SampleReadFilter \\
            ${combined_interval_options} \\
            --emit-original-quals ${params.is_emit_original_quals} \\
            --output /dev/stdout \\
            --sample ${it} 2> .command.err | \\
            samtools view -h | \\
            awk '(/^@RG/ && /SM:${it}/) || ! /^@RG/' | \\
            samtools view -b -o ${it}_recalibrated-${interval_id}.bam
        """
        }
        .join("\n")
    """
    set -euo pipefail
    ${all_commands}
    """
}

workflow recalibrate_base {
    take:
    input_samples  // Channel of sample data with BAM, index, interval info

    main:
    // Group by sample_id and interval_id to process each sample's intervals together
    input_samples
        .map { sample ->
            [
                sample.sample_id,
                sample.interval_id,
                sample.bam,
                sample.bam_index,
                sample.interval,
                sample.has_unmapped
            ]
        }
        .groupTuple(by: [0,1])  // Group by sample_id and interval_id
        .map { sample_id, interval_id, bams, bam_indices, intervals, has_unmapped ->
            [
                bams.flatten(),
                bam_indices.flatten(),
                sample_id,
                interval_id,
                intervals[0]  // Take first interval since they're all the same after grouping
            ]
        }
        .set{ input_ch_base_recalibrator }

    base_recalibrator_intervals = (params.is_targeted) ? params.intervals : "/scratch/NO_FILE.interval_list"

    run_BaseRecalibrator_GATK(
        params.reference_fasta,
        params.reference_fasta_fai,
        params.reference_fasta_dict,
        params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz,
        params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz_tbi,
        params.bundle_known_indels_vcf_gz,
        params.bundle_known_indels_vcf_gz_tbi,
        params.bundle_v0_dbsnp138_vcf_gz,
        params.bundle_v0_dbsnp138_vcf_gz_tbi,
        base_recalibrator_intervals,
        params.input.recalibration_table,
        input_ch_base_recalibrator
    )

    // Gather BQSR reports per sample
    run_BaseRecalibrator_GATK.out.recalibration_table
        .map { sample_id, recal_table -> 
            [sample_id, recal_table]
        }
        .groupTuple()
        .set{ input_ch_gather_bqsr }

    run_GatherBQSRReports_GATK(input_ch_gather_bqsr)

    // Combine input samples with their gathered recalibration tables for ApplyBQSR
    input_samples
        .map { sample ->
            [
                sample.sample_id,
                [
                    'bam': sample.bam,
                    'bam_index': sample.bam_index,
                    'interval_id': sample.interval_id,
                    'interval': sample.interval,
                    'has_unmapped': sample.has_unmapped
                ]
            ]
        }
        .combine(run_GatherBQSRReports_GATK.out.gathered_recalibration_table, by: 0)  // Join by sample_id
        .map{ sample_id, sample_data, recal_table ->
            [
                sample_data.bam,
                sample_data.bam_index,
                sample_data.interval_id,
                sample_data.interval,
                sample_data.has_unmapped,
                [sample_id],  // Wrap in list since ApplyBQSR expects a list of sample IDs
                recal_table
            ]
        }
        .set{ input_ch_apply_bqsr }

    run_ApplyBQSR_GATK(
        params.reference_fasta,
        params.reference_fasta_fai,
        params.reference_fasta_dict,
        input_ch_apply_bqsr
    )

    run_ApplyBQSR_GATK.out.output_ch_apply_bqsr
        .flatten()
        .map{ bam ->
            def sample_id = bam.getName().split('_recalibrated-')[0]
            [
                'sample': sample_id,
                'bam': bam
            ]
        }
        .set{ recalibrated_samples }

    emit:
    recalibrated_samples
}
