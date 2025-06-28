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
    run_CalculateContamination_GATK
    run_DepthOfCoverage_GATK
} from './module/summary-qc.nf'
include { delete_input } from './module/delete-input.nf'
include { sanitize_string } from './external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

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
        .set{ input_ch_intervals }


    /**
    *   Conditional workflow processing
    */
    
    // Prepare input for processing
    validated_samples_with_index
        .reduce( ['bams': [], 'indices': []] ){ a, b ->
            a.bams.add(b.path);
            a.indices.add(b.index);
            return a
        }
        .combine(input_ch_intervals)
        .map{ it -> it[0] + it[1] }
        .set{ input_ch_processing }

    // Execute processing based on parameters
    if (params.run_indel_realignment && params.run_base_recalibration) {
        // Both processes (original behavior)
        realign_indels(input_ch_processing)
        
        // Input file deletion after indel realignment
        validated_samples_with_index
            .filter{ params.metapipeline_states_to_delete.contains(it.sample_type) }
            .map{ sample -> sample.path }
            .flatten()
            .unique()
            .set{ input_bams }

        realign_indels.out.output_ch_realign_indels
            .map{ it.has_unmapped }
            .collect()
            .set{ ir_complete_signal }

        input_bams
            .combine(ir_complete_signal)
            .map{ it[0] }
            .set{ input_ch_bams_to_delete }

        delete_input(input_ch_bams_to_delete)
        
        recalibrate_base(
            realign_indels.out.output_ch_realign_indels,
            validated_samples_with_index.map{ sample -> sample.id }.flatten()
        )
        
        samples_for_merge = recalibrate_base.out.recalibrated_samples
        
    } else if (params.run_indel_realignment && !params.run_base_recalibration) {
        // Indel realignment only
        realign_indels(input_ch_processing)
        
        // Input file deletion after indel realignment
        validated_samples_with_index
            .filter{ params.metapipeline_states_to_delete.contains(it.sample_type) }
            .map{ sample -> sample.path }
            .flatten()
            .unique()
            .set{ input_bams }

        realign_indels.out.output_ch_realign_indels
            .map{ it.has_unmapped }
            .collect()
            .set{ ir_complete_signal }

        input_bams
            .combine(ir_complete_signal)
            .map{ it[0] }
            .set{ input_ch_bams_to_delete }

        delete_input(input_ch_bams_to_delete)
        
        // Transform indel realignment output to merge format
        realign_indels.out.output_ch_realign_indels
            .map{
                [
                    'sample': it.bam.baseName.split("_")[3], // extract sample ID
                    'bam': it.bam
                ]
            }
            .set{ samples_for_merge }
            
    } else if (!params.run_indel_realignment && params.run_base_recalibration) {
        // Base recalibration only on raw BAMs
        // Create interval-based entries to enable ApplyBQSR parallelization
        // BaseRecalibrator will deduplicate identical BAMs automatically
        validated_samples_with_index
            .combine(input_ch_intervals)
            .map{ sample, interval ->
                [
                    'bam': sample.path,
                    'bam_index': sample.index,
                    'interval_id': interval.interval_id,
                    'interval': interval.interval_path,
                    'has_unmapped': (interval.interval_id == 'nonassembled' || interval.interval_id == '0000'),
                    'sample_id': sample.id  // Preserve sample ID with each BAM
                ]
            }
            .set{ raw_bam_input }
        
        recalibrate_base(
            raw_bam_input,
            validated_samples_with_index.map{ sample -> sample.id }.flatten()
        )
        
        // Input file deletion (consistent with other modes)
        validated_samples_with_index
            .filter{ params.metapipeline_states_to_delete.contains(it.sample_type) }
            .map{ sample -> sample.path }
            .flatten()
            .unique()
            .set{ input_bams }

        recalibrate_base.out.recalibrated_samples
            .map{ it.sample } // Use any signal from base recalibration completion
            .collect()
            .set{ bqsr_complete_signal }

        input_bams
            .combine(bqsr_complete_signal)
            .map{ it[0] }
            .set{ input_ch_bams_to_delete }

        delete_input(input_ch_bams_to_delete)
        
        samples_for_merge = recalibrate_base.out.recalibrated_samples
    }

    /**
    *   Merge BAMs
    */
    merge_bams(samples_for_merge)


    /**
    *   Summary and QC processes
    */
    merge_bams.out.output_ch_merge_bams
        .map{ [it.sample, it.bam, it.bam_index] }
        .set{ input_ch_merged_bams }

    summary_intervals = (params.is_targeted) ?
        Channel.from(params.intervals).collect() :
        extract_GenomeIntervals.out.genomic_intervals

    summary_intervals.combine(input_ch_merged_bams)
        .map{ it[0] }
        .set{ input_ch_summary_intervals }

    run_GetPileupSummaries_GATK(
        params.reference_fasta,
        params.reference_fasta_fai,
        params.reference_fasta_dict,
        params.bundle_contest_hapmap_3p3_vcf_gz,
        params.bundle_contest_hapmap_3p3_vcf_gz_tbi,
        input_ch_summary_intervals,
        input_ch_merged_bams
    )

    validated_samples_with_index
        .filter{ it.sample_type == 'normal' }
        .map{ it -> [sanitize_string(it.id)] }
        .join(run_GetPileupSummaries_GATK.out.pileupsummaries)
        .set{ normal_pileupsummaries }

    validated_samples_with_index
        .filter{ it.sample_type == 'tumor' }
        .map{ it -> [sanitize_string(it.id)] }
        .join(run_GetPileupSummaries_GATK.out.pileupsummaries)
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

    run_DepthOfCoverage_GATK(
        params.reference_fasta,
        params.reference_fasta_fai,
        params.reference_fasta_dict,
        input_ch_summary_intervals,
        input_ch_merged_bams
    )
}
