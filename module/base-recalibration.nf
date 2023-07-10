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
        bundle_mills_and_1000g_gold_standards_vcf_gz: path to standard Mills and 1000 genomes variants
        bundle_mills_and_1000g_gold_standards_vcf_gz_tbi: path to index file for Mills and 1000g variants
        bundle_known_indels_vcf_gz: path to set of known indels
        bundle_known_indels_vcf_gz_tbi: path to index of known indels VCF
        bundle_v0_dbsnp138_vcf_gz: path to dbSNP variants
        bundle_v0_dbsnp138_vcf_gz_tbi: path to index of dbSNP variants
        all_intervals: list or tuple of paths to all split intervals
        indelrealigned_bams: list or tuple of paths to indel realigned BAMs
        indelrealigned_bams_bai: list or tuple of paths to indel realigned BAM indices
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

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${sample_id}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(bundle_mills_and_1000g_gold_standards_vcf_gz)
    path(bundle_mills_and_1000g_gold_standards_vcf_gz_tbi)
    path(bundle_known_indels_vcf_gz)
    path(bundle_known_indels_vcf_gz_tbi)
    path(bundle_v0_dbsnp138_vcf_gz)
    path(bundle_v0_dbsnp138_vcf_gz_tbi)
    tuple path(indelrealigned_bams), path(indelrealigned_bams_bai), val(sample_id)

    output:
    path(".command.*")
    tuple val(sample_id), path("${sample_id}_recalibration_table.grp"), emit: recalibration_table

    script:
    all_ir_bams = indelrealigned_bams.collect{ "--input '${it}'" }.join(' ')
    targeted_options = params.is_targeted ? "--intervals ${params.intervals} --interval-padding 100" : ""
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
        BaseRecalibrator \
        ${all_ir_bams} \
        --reference ${reference_fasta} \
        --verbosity INFO \
        --known-sites ${bundle_mills_and_1000g_gold_standards_vcf_gz} \
        --known-sites ${bundle_known_indels_vcf_gz} \
        --known-sites ${bundle_v0_dbsnp138_vcf_gz} \
        --output ${sample_id}_recalibration_table.grp \
        ${targeted_options} \
        --read-filter SampleReadFilter \
        --sample ${sample_id}
    """
}

/*
    Nextflow module for recalibrating base quality scores in BAM file

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        indelrealigned_bam: list or tuple of paths to indel realigned BAMs
        indelrealigned_bam_index: list or tuple of paths to indel realigned BAM indices
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

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${interval_id}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    tuple path(indelrealigned_bam),
          path(indelrealigned_bam_index),
          val(interval_id),
          path(interval),
          val(includes_unmapped),
          val(sample_id),
          path(recalibration_table)

    output:
    path(".command.*")
    path("*_recalibrated-*"), emit: output_ch_apply_bqsr
    tuple path(indelrealigned_bam), path(indelrealigned_bam_index), emit: output_ch_deletion

    script:
    unmapped_interval_option = (includes_unmapped) ? "--intervals unmapped" : ""
    combined_interval_options = "--intervals ${interval} ${unmapped_interval_option}"
    all_commands = sample_id.collect{
        """
        gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \\
            ApplyBQSR \\
            --input ${indelrealigned_bam} \\
            --bqsr-recal-file ${it}_recalibration_table.grp \\
            --reference ${reference_fasta} \\
            --read-filter SampleReadFilter \\
            ${combined_interval_options} \\
            --emit-original-quals ${params.is_emit_original_quals} \\
            --output /dev/stdout \\
            --sample ${it} 2> .command.err | \\
            samtools view -h | \\
            awk '(/^@RG/ && /SM:${it}/) || ! /^@RG/' | \\
            samtools view -b -o ${generate_standard_filename(params.aligner, params.dataset_id, it, ['additional_tools': ["GATK-${params.gatk_version}"]])}_recalibrated-${interval_id}.bam
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
    ir_samples
    sample_ids

    main:
    ir_samples
        .reduce( ['bams': [], 'indices': []] ){ a, b ->
            a.bams.add(b.bam);
            a.indices.add(b.bam_index);
            return a
        }
        .combine(sample_ids)
        .map{ it ->
            [
                it[0].bams,
                it[0].indices,
                it[1]
            ]
        }
        .set{ input_ch_base_recalibrator }

    run_BaseRecalibrator_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz,
        "${params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz}.tbi",
        params.bundle_known_indels_vcf_gz,
        "${params.bundle_known_indels_vcf_gz}.tbi",
        params.bundle_v0_dbsnp138_vcf_gz,
        "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
        input_ch_base_recalibrator
    )

    run_BaseRecalibrator_GATK.out.recalibration_table
        .reduce( ['ids': [], 'tables': []] ){ a, b ->
            a.ids.add(b[0]);
            a.tables.add(b[1]);
            return a
        }
        .set{ all_recal_tables }

    ir_samples
        .combine(all_recal_tables)
        .map{ it ->
            [
                it[0].bam,
                it[0].bam_index,
                it[0].interval_id,
                it[0].interval,
                it[0].has_unmapped,
                it[1].ids,
                it[1].tables
            ]
        }
        .set{ input_ch_apply_bqsr }

    run_ApplyBQSR_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        input_ch_apply_bqsr
    )

    run_ApplyBQSR_GATK.out.output_ch_apply_bqsr
        .flatten()
        .map{
            [
                'sample': it.baseName.split("_")[3], // sample ID
                'bam': it
            ]
        }
        .set{ output_ch_base_recalibration }

    remove_intermediate_files(
        run_ApplyBQSR_GATK.out.output_ch_deletion.flatten(),
        "bqsr_complete"
    )

    emit:
    recalibrated_samples = output_ch_base_recalibration
}
