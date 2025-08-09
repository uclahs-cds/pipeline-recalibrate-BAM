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
    path(intervals)
    tuple val(sample_id), path(input_bam), path(input_bam_bai), val(interval_id), path(interval_path), val(has_unmapped)

    output:
    tuple val(sample_id), path("${sample_id}_${interval_id}_recalibration_table.grp"), emit: recalibration_table

    script:
    base_interval_option = "--intervals ${interval_path}"
    unmapped_interval_option = has_unmapped ? "--intervals unmapped" : ""
    interval_padding = params.is_targeted ? "--interval-padding 100" : ""
    combined_interval_options = params.is_targeted ?
        "--intervals \"\$(realpath ${intervals})\" ${interval_padding} --interval-set-rule INTERSECTION" :
        "${unmapped_interval_option} --interval-set-rule UNION"
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
        --read-filter SampleReadFilter \
        --sample ${sample_id} \
        ${base_interval_option} \
        ${combined_interval_options} || touch ${sample_id}_${interval_id}_recalibration_table.grp
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
    tuple val(sample_id), path("${output_filename}.bam"), path("${output_filename}.bam.bai"), val(interval_id), val(interval), val(includes_unmapped), emit: output_ch_apply_bqsr
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
            'additional_information': "recalibrated-${interval_id}"
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
        --output ${output_filename}.bam \
        --create-output-bam-index true \
        --sample ${sample_id}

    mv ${output_filename}.bai ${output_filename}.bam.bai
    """
}

workflow recalibrate_base {
    take:
    input_samples

    main:
    /**
    *   BaseRecalibrator
    */
    def given_recal_tables = [:]
    params.input.recalibration_table.each { recal_table ->
        given_recal_tables["${file(recal_table).baseName.replace('_recalibration_table.grp', '')}"] = recal_table
    }

    input_samples
        .filter{ s -> (!given_recal_tables.containsKey(s.id)) }
        .map{ ir_sample ->
            [
                ir_sample.id,
                ir_sample.bam,
                ir_sample.bam_index,
                ir_sample.interval_id,
                ir_sample.interval_path,
                ir_sample.has_unmapped
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
        input_ch_base_recalibrator
    )

    /**
    *   GatherBQSRReports
    */
    run_BaseRecalibrator_GATK.out.recalibration_table
        .filter{ recal_table ->
            def f = new File(recal_table[1].toString());
            f.length() != 0 // Filter out empty files caused by disjoint intervals
        }
        .map{ recal_table ->
            [
                recal_table[0], // Sample ID
                recal_table[1]
            ]
        }
        .groupTuple(by: 0)
        .set{ input_ch_gatherbqsr }

    run_GatherBQSRReports_GATK(input_ch_gatherbqsr)


    /**
    *   ApplyBQSR
    */
    def recal_tables = []
    given_recal_tables.each{ s, t ->
        recal_tables << [s, t]
    }

    run_GatherBQSRReports_GATK.out.gathered_recalibration_table
        .mix(Channel.fromList(recal_tables))
        .map{ tables ->
            [
                tables[0],
                ['recal_table': tables[1]]
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
        .combine(ided_input_samples, by: 0)
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

    // sample_id, bam, bam_index, interval_id, interval, includes_unmapped
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
