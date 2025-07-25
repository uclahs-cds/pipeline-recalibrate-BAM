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

    Channel.from(params.samples_to_process)
        .map { sample ->
            sample["index"] = indexFile(sample.path)
            return sample
        }
        .set{ samples_with_index }

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
        .map{ interval_data ->
            interval_data['has_unmapped'] = (interval_data.interval_id == 'nonassembled' || interval_data.interval_id == '0000');
            return interval_data;
        }
        .set{ input_ch_intervals }


    /**
    *   Indel realignment
    */
    if (params.run_indelrealignment) {
        samples_with_index
            .reduce( ['bams': [], 'indices': []] ){ a, b ->
                a.bams.add(b.path);
                a.indices.add(b.index);
                return a
            }
            .combine(input_ch_intervals)
            .map{ it -> it[0] + it[1] }
            .set{ input_ch_indel_realignment }

        realign_indels(input_ch_indel_realignment)

        output_ch_ir = realign_indels.out.output_ch_realign_indels

        realign_indels.out.output_ch_realign_indels
            .set{ output_ch_ir }

        /**
        *   Input file deletion
        */
        samples_with_index
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
    } else {
        // Reformat channels to match expected output of IndelRealignment
        samples_with_index
            .map{ raw_samples ->
                [
                    'bam': [raw_samples.path],
                    'bam_index': [raw_samples.index],
                    'id': raw_samples.id
                ]
            }
            .combine(input_ch_intervals)
            .map{ it -> it[0] + it[1] }
            .set{ output_ch_ir }
    }

    /**
    *   Base recalibration
    */
    if (params.run_bqsr) {
        recalibrate_base(
            output_ch_ir
        )

        output_ch_bqsr = recalibrate_base.out.recalibrated_samples
    } else {
        // Pass through IR output channel since formats now match
        output_ch_bqsr = output_ch_ir
    }


    /**
    *   Merge BAMs
    */
    merge_bams(output_ch_bqsr)


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

    samples_with_index
        .filter{ it.sample_type == 'normal' }
        .map{ it -> [sanitize_string(it.id)] }
        .join(run_GetPileupSummaries_GATK.out.pileupsummaries)
        .set{ normal_pileupsummaries }

    samples_with_index
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
