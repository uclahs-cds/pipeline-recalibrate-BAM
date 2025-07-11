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
include { delete_input } from './delete-input.nf'
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
        input_bams: list or tuple of paths to input  BAMs
        input_bams_bai: list or tuple of paths to input BAM indices
        sample_id: identifier for the sample
        
    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
        params.is_targeted: bool. Indicator of whether in targeted exome mode or in WGS mode
        params.gatk_command_mem_diff: float(memory)
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
    tuple path(input_bams), path(input_bams_bai), val(sample_id)

    output:
    tuple val(sample_id), path("${sample_id}_recalibration_table.grp"), emit: recalibration_table

    script:
    all_input_bams = input_bams.collect{ "--input '${it}'" }.join(' ')
    targeted_options = params.is_targeted ? "--intervals \"\$(realpath ${intervals})\" --interval-padding 100" : ""
    """
    set -euo pipefail
    if [ ! -f ${sample_id}_recalibration_table.grp ]
    then
        gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
            BaseRecalibrator \
            ${all_input_bams} \
            --reference ${reference_fasta} \
            --verbosity INFO \
            --known-sites ${bundle_mills_and_1000g_gold_standard_indels_vcf_gz} \
            --known-sites ${bundle_known_indels_vcf_gz} \
            --known-sites ${bundle_v0_dbsnp138_vcf_gz} \
            --output ${sample_id}_recalibration_table.grp \
            ${targeted_options} \
            --read-filter SampleReadFilter \
            --sample ${sample_id}
    fi
    """
}

/*
    Nextflow module for recalibrating base quality scores in BAM file

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        input_bam: path to input BAM
        input_bam_index: path to input BAM index
        interval: path to specific intervals file associated with input BAM
        includes_unmapped: boolean to indicate if unmapped reads are included in input BAM
        sample_id: Identifier for sample being processed
        recalibration_table: path to base recalibration table

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
          val(includes_unmapped),
          val(sample_id),
          path(recalibration_table)

    output:
    tuple val(sample_id), path(output_filename), path("${output_filename}.bai"), val(interval_id), val(interval), val(includes_unmapped), emit: output_ch_apply_bqsr
    tuple path(input_bam), path(input_bam_index), emit: output_ch_deletion

    script:
    unmapped_interval_option = (includes_unmapped) ? "--intervals unmapped" : ""
    combined_interval_options = "--intervals ${interval} ${unmapped_interval_option}"
    output_filename = generate_standard_filename(
        params.aligner,
        params.dataset_id,
        sample_id,
        [
            'additional_tools': ["GATK-${params.gatk_version}"],
            'additional_information': "recalibrated-${interval_id}.bam"
        ]
    )
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
    """
}

workflow recalibrate_base {
    take:
    input_samples

    main:
    input_samples
        .map{ ir_sample ->
            [it_sample.id, ir_sample.bam, it_sample.bam_index]
        }
        .groupTuple(by: 0)
        .map{ ir_grouped ->
            [
                it[1].unique{ bam_path -> file(bam_path).getFileName() },
                it[2].unique{ bam_index -> file(bam_index).getFileName() },
                it[0]
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

    run_BaseRecalibrator_GATK.out.recalibration_table
        .map{ recal_table ->
            [
                recal_table[0], // Sample ID
                ['recal_table': recal_table[1]]
            ]
        }
        .set{ ided_recal_tables }

    input_samples
        .map{ input_sample ->
            [
                input_sample.id,
                input_sample
            ]
        }
        .set{ ided_input_samples }

    ided_recal_tables
        .join(ided_input_samples)
        .map{ joined_input ->
            [
                joined_input[2].bam,
                joined_input[2].bam_index,
                joined_input[2].interval_id,
                joined_input[2].interval_path,
                joined_input[2].has_unmapped,
                joined_input[2].id,
                joined_input[1].recal_table

            ]
        }
        .set{ input_ch_apply_bqsr }

    run_ApplyBQSR_GATK(
        params.reference_fasta,
        params.reference_fasta_fai,
        params.reference_fasta_dict,
        input_ch_apply_bqsr
    )

    // val(sample_id), path(output_filename), path("${output_filename}.bai"), val(interval_id), val(interval), val(includes_unmapped)
    run_ApplyBQSR_GATK.out.output_ch_apply_bqsr
        .map{ bqsred_sample ->
            [
                'id': bqsred_sample[0],
                'bam': bqsred_sample[1],
                'bam_index': bqsred_sample[2],
                'interval_id': bqsred_sample[3],
                'interval_path': bqsred_sample[4],
                'has_unmapped': bqsred_sample[5]
            ]
        }
        .set{ output_ch_base_recalibration }

    /**
    *   Handle deletion
    */
    if (params.run_indelrealignment) {
        // Delete input to BQSR as intermediate files
        remove_intermediate_files(
            run_ApplyBQSR_GATK.out.output_ch_deletion.flatten(),
            "bqsr_complete"
        )
    } else {
        // Pipeline inputs directly sent to BQSR, handle deletion based on metapipeline params
        Channel.from(params.samples_to_process)
            .filter{ params.metapipeline_states_to_delete.contains(it.sample_type) }
            .map{ sample -> sample.path }
            .flatten()
            .unique()
            .set{ input_bams_to_delete }

        run_ApplyBQSR_GATK.out.output_ch_apply_bqsr
            .map{ 'done' }
            .collect()
            .map{ 'done' }
            .set{ bqsr_complete_signal }

        input_bams_to_delete
            .combine(bqsr_complete_signal)
            .map{ it[0] }
            .set{ input_ch_bams_to_delete }

        delete_input(input_ch_bams_to_delete)
    }

    emit:
    recalibrated_samples = output_ch_base_recalibration
}
