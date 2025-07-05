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
    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "*.grp"

    ext log_dir_suffix: { "-${sample_id}-${interval_id}" }

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
    tuple path(input_bam),
          path(input_bam_index),
          val(interval_id),
          path(interval),
          val(includes_unmapped),
          val(sample_id)

    output:
    tuple val(sample_id), val(interval_id), path("${sample_id}_${interval_id}_recalibration_table.grp"), emit: recalibration_table

    script:
    unmapped_interval_option = (includes_unmapped) ? "--intervals unmapped" : ""
    combined_interval_options = "--intervals ${interval} ${unmapped_interval_option}"
    targeted_options = params.is_targeted ? "--interval-padding 100" : ""
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
        BaseRecalibrator \
        --input ${input_bam} \
        --reference ${reference_fasta} \
        --verbosity INFO \
        --known-sites ${bundle_mills_and_1000g_gold_standard_indels_vcf_gz} \
        --known-sites ${bundle_known_indels_vcf_gz} \
        --known-sites ${bundle_v0_dbsnp138_vcf_gz} \
        --output ${sample_id}_${interval_id}_recalibration_table.grp \
        ${combined_interval_options} \
        ${targeted_options} \
        --read-filter SampleReadFilter \
        --sample ${sample_id}
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
        includes_noncanonical_contigs: boolean to indicate if noncanonical contigs are included in input BAM
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

    ext log_dir_suffix: { "-${sample_id}-${interval_id}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    tuple path(input_bam),
          path(input_bam_index),
          val(interval_id),
          path(interval),
          val(includes_noncanonical_contigs),
          val(sample_id),
          path(recalibration_table)

    output:
    tuple path("${output_filename}"), path("${output_filename}.bai"), val(interval_id), path(interval), val(includes_noncanonical_contigs), val(sample_id), emit: output_ch_apply_bqsr
    tuple path(input_bam), path(input_bam_index), emit: output_ch_deletion

    script:
    unmapped_interval_option = (includes_noncanonical_contigs) ? "--intervals unmapped" : ""
    combined_interval_options = "--intervals ${interval} ${unmapped_interval_option}"
    output_filename = "${generate_standard_filename(params.aligner, params.dataset_id, sample_id, ['additional_tools': ["GATK-${params.gatk_version}"]])}_recalibrated-${interval_id}.bam"
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
        ApplyBQSR \
        --input ${input_bam} \
        --bqsr-recal-file ${recalibration_table} \
        --reference ${reference_fasta} \
        --read-filter SampleReadFilter \
        ${combined_interval_options} \
        --emit-original-quals ${params.is_emit_original_quals} \
        --output /dev/stdout \
        --sample ${sample_id} 2> .command.err | \
        samtools view -h | \
        awk '(/^@RG/ && /SM:${sample_id}/) || ! /^@RG/' | \
        samtools view -b -o ${output_filename}
    
    # Index the output BAM file for GetPileupSummaries
    samtools index ${output_filename}
    """
}

workflow recalibrate_base {
    take:
    input_samples  // Channel of sample data with BAM, index, interval info, and sample_id

    main:
    // Extract sample ID from BAM for each sample - handle both modes
    input_samples
        .map{ sample_data ->
            def sample_id

            if (sample_data.containsKey('sample_id')) {
                // Both base recalibration only mode and indel realignment mode now provide sample_id
                sample_id = sample_data.sample_id
            } else {
                // Fallback - this should not happen with the new implementation
                throw new Exception("Sample ID not found in input data: ${sample_data}")
            }

            return [
                sample_data.bam,
                sample_data.bam_index,
                sample_data.interval_id,
                sample_data.interval,
                sample_data.is_noncanonical_contig,
                sample_id
            ]
        }
        .set{ input_ch_base_recalibrator }

    // Run BaseRecalibrator in parallel for each interval and sample
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
        input_ch_base_recalibrator
    )

    // Group recalibration tables by sample and gather them
    run_BaseRecalibrator_GATK.out.recalibration_table
        .map{ sample_id, interval_id, table -> [sample_id, table] }
        .groupTuple()
        .set{ grouped_recal_tables }

    run_GatherBQSRReports_GATK(grouped_recal_tables)

    // Combine input samples with their gathered recalibration tables
    input_ch_base_recalibrator
        .map{ bam, bam_index, interval_id, interval, is_noncanonical_contig, sample_id ->
            [sample_id, [bam, bam_index, interval_id, interval, is_noncanonical_contig, sample_id]]
        }
        .combine(
            run_GatherBQSRReports_GATK.out.gathered_recalibration_table,
            by: 0  // Join by sample_id
        )
        .map{ sample_id, sample_data, gathered_table ->
            sample_data + [gathered_table]  // Add the gathered table path (gathered_table is already just the path)
        }
        .set{ input_ch_apply_bqsr }

    run_ApplyBQSR_GATK(
        params.reference_fasta,
        params.reference_fasta_fai,
        params.reference_fasta_dict,
        input_ch_apply_bqsr
    )

    run_ApplyBQSR_GATK.out.output_ch_apply_bqsr
        .map{ bam_file, bai_file, interval_id, interval, is_noncanonical_contig, sample_id ->
            [
                'sample': sample_id,
                'bam': bam_file
            ]
        }
        .set{ output_ch_base_recalibration }

    // Also emit interval BAMs with their interval information for pileup summaries
    run_ApplyBQSR_GATK.out.output_ch_apply_bqsr
        .map{ bam_file, bai_file, interval_id, interval, is_noncanonical_contig, sample_id ->
            [
                'sample': sample_id,
                'bam': bam_file,
                'bam_index': bai_file,
                'interval_id': interval_id,
                'interval_path': interval
            ]
        }
        .set{ output_ch_interval_bams }

    // Handle intermediate file deletion
    // Only delete files if they came from indel realignment (not original input BAMs)
    if (params.run_indel_realignment) {
        // Indel realignment was run, safe to delete intermediate files
        run_ApplyBQSR_GATK.out.output_ch_deletion
            .flatten()
            .set{ files_to_delete }
    } else {
        // Base recalibration only mode - don't delete original input BAMs
        Channel.empty().set{ files_to_delete }
    }

    remove_intermediate_files(
        files_to_delete,
        "bqsr_complete"
    )

    emit:
    recalibrated_samples = output_ch_base_recalibration
    interval_bams = output_ch_interval_bams
}
