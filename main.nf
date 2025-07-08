#!/usr/bin/env nextflow

nextflow.enable.dsl=2


log.info """\
==============================================
R E C A L I B R A T E - B A M  P I P E L I N E
==============================================
Boutros Lab

Current Configuration:

    - pipeline:
        name: ${workflow.manifest.name}
        version: ${workflow.manifest.version}

    - processing:
        run_indel_realignment: ${params.run_indel_realignment}
        run_base_recalibration: ${params.run_base_recalibration}

    - input:
        samples: ${params.samples_to_process}
        aligner: ${params.aligner}
        bundle_mills_and_1000g_gold_standard_indels_vcf_gz: ${params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz}
        bundle_known_indels_vcf_gz: ${params.bundle_known_indels_vcf_gz}
        bundle_v0_dbsnp138_vcf_gz: ${params.bundle_v0_dbsnp138_vcf_gz}
        bundle_contest_hapmap_3p3_vcf_gz: ${params.bundle_contest_hapmap_3p3_vcf_gz}
        intervals: ${(params.is_targeted) ?: 'WGS'}
        Recalibration tables: ${params.input.recalibration_table}

    - output:
        output: ${params.output_dir}
        output_dir_base: ${params.output_dir_base}
        log_output_dir: ${params.log_output_dir}

    Tools Used:
        tool GATK: ${params.docker_image_gatk}
        tool PipeVal: ${params.docker_image_pipeval}
        tool GATK3: ${params.docker_image_gatk3}

    All parameters:
        ${params}

------------------------------------
Starting workflow...
------------------------------------
        """

include { run_validate_PipeVal_with_metadata } from './external/pipeline-Nextflow-module/modules/PipeVal/validate/main.nf' addParams(
    options: [
        docker_image_version: params.pipeval_version,
        main_process: "./" //Save logs in <log_dir>/process-log/run_validate_PipeVal
        ]
    )
include { run_SplitIntervals_GATK } from './module/split-intervals.nf'
include { extract_GenomeIntervals } from './external/pipeline-Nextflow-module/modules/common/extract_genome_intervals/main.nf' addParams(
    options: [
        save_intermediate_files: params.save_intermediate_files,
        output_dir: params.output_dir_base
        ]
    )
include { realign_indels } from './module/indel-realignment.nf'
include { recalibrate_base } from './module/base-recalibration.nf'
include { merge_bams } from './module/merge-bam.nf'
include {
    run_GetPileupSummaries_GATK
    concatenate_PileupSummaries
    run_CalculateContamination_GATK
    run_DepthOfCoverage_GATK
} from './module/summary-qc.nf'
include { delete_input } from './module/delete-input.nf'
include { sanitize_string } from './external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
include { remove_intermediate_files } from './external/pipeline-Nextflow-module/modules/common/intermediate_file_removal/main.nf' addParams(
    options: [
        save_intermediate_files: params.save_intermediate_files,
        output_dir: params.output_dir_base,
        log_output_dir: "${params.log_output_dir}/process-log"
        ]
    )

// Returns the index file for the given bam or vcf
def indexFile(bam_or_vcf) {
    if (bam_or_vcf.endsWith('.bam')) {
        return "${bam_or_vcf}.bai"
    } else if (bam_or_vcf.endsWith('vcf.gz')) {
        return "${bam_or_vcf}.tbi"
    } else {
        throw new Exception("Index file for ${bam_or_vcf} file type not supported. Use .bam or .vcf.gz files.")
    }
}

workflow {
    /**
    *   Parameter validation
    */
    if (!params.run_indel_realignment && !params.run_base_recalibration) {
        error "Error: At least one of run_indel_realignment or run_base_recalibration must be true"
    }

    /**
    *   Input channel processing
    */

    /**
    *   Input validation
    */
    Channel.from(params.samples_to_process)
        .flatMap { sample ->
            def all_metadata = sample.findAll { it.key != "path" }
            return [
                [sample.path, [all_metadata, "path"]],
                [indexFile(sample.path), [[id: sample.id], "index"]]
            ]
        } | run_validate_PipeVal_with_metadata

    run_validate_PipeVal_with_metadata.out.validation_result
        .collectFile(
            name: 'input_validation.txt',
            storeDir: "${params.output_dir_base}/validation"
        )

    run_validate_PipeVal_with_metadata.out.validated_file
        .map { filename, metadata -> [metadata[0].id, metadata[0] + [(metadata[1]): filename]] }
        .groupTuple()
        .map { it[1].inject([:]) { result, i -> result + i } }
        .set { validated_samples_with_index }

    // The elements of validated_samples_with_index are the same as
    // params.samples_to_process, with the following changes:
    // * sample.path is the validated BAM file
    // * sample.index is the validated BAI file (new key)

    /**
    *   Interval extraction and splitting
    */
    extract_GenomeIntervals(params.reference_fasta_dict)

    run_SplitIntervals_GATK(
        extract_GenomeIntervals.out.genomic_intervals,
        params.reference_fasta,
        params.reference_fasta_fai,
        params.reference_fasta_dict
    )

    run_SplitIntervals_GATK.out.interval_list
        .flatten()
        .map{ interval_path ->
            [
                'interval_id': file(interval_path).getName().replace('-contig.interval_list', ''),
                'interval_path': interval_path
            ]
        }
        .filter{ it.interval_id == 'noncanonical' }
        .concat(
            run_SplitIntervals_GATK.out.interval_list
                .flatten()
                .map{ interval_path ->
                    [
                        'interval_id': file(interval_path).getName().replace('-contig.interval_list', ''),
                        'interval_path': interval_path
                    ]
                }
                .filter{ it.interval_id != 'noncanonical' }
        )
        .set{ input_ch_intervals }


    /**
    *   Conditional workflow processing
    */

    // Prepare input for processing with individual sample information preserved
    validated_samples_with_index
        .combine(input_ch_intervals)
        .map{ sample, interval ->
            [
                'bams': [sample.path],          // Keep as list for compatibility but preserve individual samples
                'indices': [sample.index],      // Keep as list for compatibility but preserve individual samples
                'interval_id': interval.interval_id,
                'interval_path': interval.interval_path,
                'sample_id': sample.id          // Preserve individual sample ID
            ]
        }
        .set{ input_ch_processing }

    /**
    *   Conditional workflow execution using channel manipulation
    */

    // Prepare channels for input deletion (common across all modes)
    validated_samples_with_index
        .filter{ params.metapipeline_states_to_delete.contains(it.sample_type) }
        .map{ sample -> sample.path }
        .flatten()
        .unique()
        .set{ input_bams_for_deletion }

    // Step 1: Conditionally run indel realignment
    if (params.run_indel_realignment) {
        realign_indels(input_ch_processing)

        processed_samples = realign_indels.out.output_ch_realign_indels
        completion_signal = realign_indels.out.output_ch_realign_indels
            .map{ it.is_noncanonical_contig }
            .collect()
    } else {
        // For base recalibration only mode, prepare raw BAM input with sample IDs
        validated_samples_with_index
            .combine(input_ch_intervals)
            .map{ sample, interval ->
                [
                    'bam': sample.path,
                    'bam_index': sample.index,
                    'interval_id': interval.interval_id,
                    'interval': interval.interval_path,
                    'is_noncanonical_contig': (interval.interval_id == 'noncanonical' || interval.interval_id == '0000'),
                    'sample_id': sample.id  // Include the actual sample ID
                ]
            }
            .set{ processed_samples }

        // No completion signal needed for raw BAMs
        Channel.empty().set{ completion_signal }
    }

    // Step 2: Conditionally run base recalibration
    if (params.run_base_recalibration) {
        // Prepare sample ID mapping for base recalibration only mode
        if (params.run_indel_realignment) {
            // For indel realignment + base recalibration, no mapping needed
            recalibrate_base(processed_samples)
        } else {
            // For base recalibration only mode, sample IDs are already included
            recalibrate_base(processed_samples)
        }

        samples_for_merge = recalibrate_base.out.recalibrated_samples
        // recalibrated_samples now contains interval_id
        processing_completion_signal = recalibrate_base.out.recalibrated_samples
            .map{ it.sample }
            .collect()
    } else {
        // Transform indel realignment output to merge format
        realign_indels.out.output_ch_realign_indels
            .map{
                [
                    'sample': it.sample_id,
                    'bam': it.bam,
                    'interval_id': it.interval_id
                ]
            }
            .set{ samples_for_merge }

        processing_completion_signal = completion_signal
    }

    // Step 3: Input deletion after processing completion
    input_bams_for_deletion
        .combine(processing_completion_signal)
        .map{ it[0] }
        .set{ input_ch_bams_to_delete }

    delete_input(input_ch_bams_to_delete)

    /**
    *   Merge BAMs
    */
    merge_bams(samples_for_merge)


    /**
    *   Summary and QC processes - Parallelized Pileup Summaries
    */
    // Always run GetPileupSummaries on interval BAMs (from either base recalibration or indel realignment)
    if (params.run_base_recalibration) {
        // Use interval BAMs from base recalibration workflow
        recalibrate_base.out.interval_bams
            .map{ interval_data ->
                [interval_data.sample, interval_data.bam, interval_data.bam_index, interval_data.interval_id, interval_data.interval_path]
            }
            .set{ input_ch_interval_pileups }
    } else {
        // Use interval BAMs from indel realignment workflow (since at least one must be run)
        realign_indels.out.interval_bams
            .map{ interval_data ->
                [interval_data.sample, interval_data.bam, interval_data.bam_index, interval_data.interval_id, interval_data.interval_path]
            }
            .set{ input_ch_interval_pileups }
    }

    run_GetPileupSummaries_GATK(
        params.reference_fasta,
        params.reference_fasta_fai,
        params.reference_fasta_dict,
        params.bundle_contest_hapmap_3p3_vcf_gz,
        params.bundle_contest_hapmap_3p3_vcf_gz_tbi,
        input_ch_interval_pileups
    )

    // Group interval results by sample and concatenate
    run_GetPileupSummaries_GATK.out.interval_pileupsummaries
        .map{ sample_id, interval_id, pileup_file -> [sample_id, pileup_file] }
        .groupTuple()
        .set{ grouped_interval_pileups }

    concatenate_PileupSummaries(grouped_interval_pileups)

    pileup_summaries_output = concatenate_PileupSummaries.out.concatenated_pileupsummaries

    validated_samples_with_index
        .filter{ it.sample_type == 'normal' }
        .map{ it -> [sanitize_string(it.id)] }
        .join(pileup_summaries_output)
        .set{ normal_pileupsummaries }

    validated_samples_with_index
        .filter{ it.sample_type == 'tumor' }
        .map{ it -> [sanitize_string(it.id)] }
        .join(pileup_summaries_output)
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

    // Setup DepthOfCoverage - always use merged BAMs regardless of pileup summaries parallelization
    merge_bams.out.output_ch_merge_bams
        .map{ [it.sample, it.bam, it.bam_index] }
        .set{ input_ch_merged_bams }

    summary_intervals = (params.is_targeted) ?
        Channel.from(params.intervals).collect() :
        extract_GenomeIntervals.out.genomic_intervals

    summary_intervals.combine(input_ch_merged_bams)
        .map{ it[0] }
        .set{ input_ch_summary_intervals }

    run_DepthOfCoverage_GATK(
        params.reference_fasta,
        params.reference_fasta_fai,
        params.reference_fasta_dict,
        input_ch_summary_intervals,
        input_ch_merged_bams
    )

    /**
    *   Per-sample coordinated deletion of interval BAMs
    *   Delete each sample's interval BAMs as soon as that sample completes both workflows
    */
    if (params.run_base_recalibration || params.run_indel_realignment) {
        // Get per-sample merge completion signals
        merge_bams.out.output_ch_merge_bams
            .map{ merge_data -> [merge_data.sample, "merge_complete"] }
            .set{ per_sample_merge_completion }
        
        // Get per-sample pileup completion signals
        concatenate_PileupSummaries.out.concatenated_pileupsummaries
            .map{ sample_id, pileup_file -> [sample_id, "pileup_complete"] }
            .set{ per_sample_pileup_completion }
        
        // Organize interval BAMs by sample for deletion
        if (params.run_base_recalibration) {
            // Use recalibrated interval BAMs
            recalibrate_base.out.interval_bams
                .map{ interval_data -> [interval_data.sample, interval_data.bam] }
                .groupTuple()
                .set{ interval_bams_by_sample }
        } else {
            // Use indel realigned interval BAMs (indel realignment only mode)
            realign_indels.out.interval_bams
                .map{ interval_data -> [interval_data.sample, interval_data.bam] }
                .groupTuple()
                .set{ interval_bams_by_sample }
        }
        
        // Per-sample coordinated deletion: delete as soon as THIS sample completes both workflows
        per_sample_merge_completion
            .join(per_sample_pileup_completion, by: 0)  // Join by sample_id, wait for both workflows
            .join(interval_bams_by_sample, by: 0)       // Get this sample's interval BAMs
            .map{ sample_id, merge_signal, pileup_signal, interval_bams -> 
                [sample_id, interval_bams]
            }
            .set{ ready_for_deletion_by_sample }
        
        // Flatten the BAMs for deletion
        ready_for_deletion_by_sample
            .map{ sample_id, interval_bams -> interval_bams }
            .flatten()
            .set{ coordinated_deletion_files }

        remove_intermediate_files(
            coordinated_deletion_files,
            "sample_workflows_complete"
        )
    }
}

