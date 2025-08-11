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
include {
    remove_intermediate_files as remove_interval_BAMs
    } from './external/pipeline-Nextflow-module/modules/common/intermediate_file_removal/main.nf' addParams(
        options: [
            save_intermediate_files: params.save_intermediate_files,
            output_dir: params.output_dir_base,
            log_output_dir: "${params.log_output_dir}/process-log"
            ]
        )
include { indexFile } from './external/pipeline-Nextflow-module/modules/common/indexFile/main.nf'
include { realign_indels } from './module/indel-realignment.nf'
include { recalibrate_base } from './module/base-recalibration.nf'
include { merge_bams } from './module/merge-bam.nf'
include {
    contamination_qc
    run_DepthOfCoverage_GATK
} from './module/summary-qc.nf'
include { delete_input } from './module/delete-input.nf'
include { sanitize_string } from './external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

workflow {
    /**
    *   Input channel processing
    */

    /**
    *   Input validation
    */
    Channel.from(params.samples_to_process)
        .map { sample ->
            sample["index"] = indexFile(sample.bam)
            return sample
        }
        .set{ samples_with_index }

    samples_with_index
        .flatMap { full_sample ->
            def all_metadata = full_sample.findAll { it.key != "bam" }
            return [
                [full_sample.bam, [all_metadata, "bam"]],
                [full_sample.index, [[id: full_sample.id], "index"]]
            ]
        } | run_validate_PipeVal_with_metadata

    run_validate_PipeVal_with_metadata.out.validation_result
        .collectFile(
            name: 'input_validation.txt',
            storeDir: "${params.output_dir_base}/validation"
        )


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
                a.bams.add(b.bam);
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
            .map{ sample -> sample.bam }
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
                    'bam': [raw_samples.bam],
                    'bam_index': [raw_samples.index],
                    'id': raw_samples.sample_id
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
    *   Contamination calculation
    */
    contamination_qc(
        output_ch_bqsr,
        samples_with_index
    )

    /**
    *   Remove interval BAMs
    */
    output_ch_bqsr
        .map{ bqsred_bams -> [bqsred_bams.id, bqsred_bams.bam] }
        .groupTuple()
        .map{ grouped_bams -> [grouped_bams[0], ['bams': grouped_bams[1]]] }
        .set{ interval_bams_to_delete }

    merge_bams.out.output_ch_merge_bams
        .map{ merged_bam -> merged_bam.sample }
        .set{ samples_merged }

    interval_bams_to_delete
        .join(samples_merged)
        .map{ joined_on_merged -> [joined_on_merged[0], joined_on_merged[1]] }
        .join(contamination_qc.out.pileupsgathered)
        .map{ joined_on_contamination -> joined_on_contamination[1]['bams'] }
        .flatten()
        .set{ input_ch_delete_interval_bams }

    remove_interval_BAMs(
        input_ch_delete_interval_bams,
        "ready_to_delete"
    )

    /**
    *   Depth of Coverage
    */
    merge_bams.out.output_ch_merge_bams
        .map{ merged_bam -> [merged_bam.sample, merged_bam.bam, merged_bam.bam_index] }
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
}
