include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
/*
    Nextflow module for generating realignment targets

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        bundle_mills_and_1000g_gold_standard_indels_vcf_gz: path to standard Mills and 1000 genomes variants
        bundle_mills_and_1000g_gold_standard_indels_vcf_gz_tbi: path to index file for Mills and 1000g variants
        bundle_known_indels_vcf_gz: path to set of known indels
        bundle_known_indels_vcf_gz_tbi: path to index of known indels VCF
        (bam, bam_index, interval_id, interval):  
            tuple of input BAM and index files, interval ID, and interval
        
    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk3: string
        params.is_targeted: bool.
        params.gatk_command_mem_diff: float(memory)
*/
process run_RealignerTargetCreator_GATK {
    container params.docker_image_gatk3
    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "*_RTC_*.intervals"

    ext log_dir_suffix: { "-${interval_id}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(bundle_mills_and_1000g_gold_standard_indels_vcf_gz)
    path(bundle_mills_and_1000g_gold_standard_indels_vcf_gz_tbi)
    path(bundle_known_indels_vcf_gz)
    path(bundle_known_indels_vcf_gz_tbi)
    path(original_intervals)
    tuple path(bam), path(bam_index), val(interval_id), path(interval), val(has_unmapped)

    output:
    tuple path(bam), path(bam_index), val(interval_id), path(interval), val(has_unmapped), path("${output_rtc_intervals}"), emit: ir_targets

    script:
    arg_bam = bam.collect{ "--input_file '${it}'" }.join(' ')
    interval_padding = params.is_targeted ? "--interval_padding 100" : ""
    targeted_interval_params = params.is_targeted ? "--intervals ${params.intervals} --interval_set_rule INTERSECTION" : ""
    output_rtc_intervals = "${params.patient_id}_RTC_${interval_id}.intervals"
    """
    set -euo pipefail
    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir} \
        -jar /GenomeAnalysisTK.jar \
        --analysis_type RealignerTargetCreator \
        ${arg_bam} \
        --reference_sequence ${reference_fasta} \
        --known ${bundle_mills_and_1000g_gold_standard_indels_vcf_gz} \
        --known ${bundle_known_indels_vcf_gz} \
        --intervals ${interval} \
        --out ${output_rtc_intervals} \
        --allow_potentially_misencoded_quality_scores \
        --num_threads ${task.cpus} \
        ${interval_padding} \
        ${targeted_interval_params} || touch ${output_rtc_intervals}
    """
}

/*
    Nextflow module for realigning indels

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        bundle_mills_and_1000g_gold_standard_indels_vcf_gz: path to standard Mills and 1000 genomes variants
        bundle_mills_and_1000g_gold_standard_indels_vcf_gz_tbi: path to index file for Mills and 1000g variants
        bundle_known_indels_vcf_gz: path to set of known indels
        bundle_known_indels_vcf_gz_tbi: path to index of known indels VCF
        (bam, bam_index, interval_id, interval, RTC_interval): 
            tuple of input BAM and index files, interval ID, target interval, and realignment targets

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk3: string
        params.is_targeted: bool. Indicator of whether in targeted exome mode or in WGS mode
        params.gatk_command_mem_diff: float(memory)
*/
process run_IndelRealigner_GATK {
    container params.docker_image_gatk3
    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "*indelrealigned*",
      saveAs: { filename -> (file(filename).getExtension() == "bai") ? "${file(filename).baseName}.bam.bai" : "${filename}" }

    ext log_dir_suffix: { "-${interval_id}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(bundle_mills_and_1000g_gold_standard_indels_vcf_gz)
    path(bundle_mills_and_1000g_gold_standard_indels_vcf_gz_tbi)
    path(bundle_known_indels_vcf_gz)
    path(bundle_known_indels_vcf_gz_tbi)
    tuple path(bam), path(bam_index), val(interval_id), path(scatter_intervals), val(has_unmapped), path(target_intervals_RTC)

    output:
    tuple path("*${output_extension}.bam"), path("*${output_extension}.bai"), val(interval_id), path(scatter_intervals), val(has_unmapped), val(output_extension), emit: output_ch_indel_realignment

    script:
    arg_bam = bam.collect{ "--input_file '$it'" }.join(' ')
    unmapped_interval_option = (has_unmapped) ? "--intervals unmapped" : ""
    combined_interval_options = "--intervals ${scatter_intervals} ${unmapped_interval_option}"
    output_extension = "indelrealigned-${interval_id}"
    """
    set -euo pipefail
    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir} \
        -jar /GenomeAnalysisTK.jar \
        --analysis_type IndelRealigner \
        ${arg_bam} \
        --reference_sequence ${reference_fasta} \
        --bam_compression ${params.gatk_ir_compression} \
        --knownAlleles ${bundle_mills_and_1000g_gold_standard_indels_vcf_gz} \
        --knownAlleles ${bundle_known_indels_vcf_gz} \
        --allow_potentially_misencoded_quality_scores \
        --targetIntervals ${target_intervals_RTC} \
        --nWayOut ${output_extension}.bam
        ${combined_interval_options}
    """
}

workflow realign_indels {
    take:
    ir_input

    main:
    ir_input
        .map{ it ->
            [
                it.bams,
                it.indices,
                it.interval_id,
                it.interval_path,
                it.has_unmapped
            ]
        }
        .set{ input_ch_rtc }

    run_RealignerTargetCreator_GATK(
        params.reference_fasta,
        params.reference_fasta_fai,
        params.reference_fasta_dict,
        params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz,
        params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz_tbi,
        params.bundle_known_indels_vcf_gz,
        params.bundle_known_indels_vcf_gz_tbi,
        "${params.getOrDefault('intervals', null) ?: params.work_dir + '/NO_FILE.bed'}",
        input_ch_rtc
        )

    run_IndelRealigner_GATK(
        params.reference_fasta,
        params.reference_fasta_fai,
        params.reference_fasta_dict,
        params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz,
        params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz_tbi,
        params.bundle_known_indels_vcf_gz,
        params.bundle_known_indels_vcf_gz_tbi,
        run_RealignerTargetCreator_GATK.out.ir_targets
        )

    def basename_map = [:]
    params.samples_to_process.each { sample_data ->
        basename_map["${file(sample_data.path).baseName}"] = sample_data.id
    }

    run_IndelRealigner_GATK.out.output_ch_indel_realignment
        .flatMap{
            def flattened_out = []
            ((it[0] in List) ? it[0] : [it[0]]).eachWithIndex{ elem, index ->
                flattened_out.add(
                    [
                        'bam': elem,
                        'bam_index': (it[1] in List) ? it[1][index] : it[1],
                        'interval_id': it[2],
                        'interval_path': it[3],
                        'has_unmapped': it[4],
                        'id': basename_map[file(elem).baseName.replace(it[5], "")]
                    ]
                )
            }

            return flattened_out;
        }
        .set{ output_ch_realign_indels }

    emit:
    output_ch_realign_indels = output_ch_realign_indels
}
