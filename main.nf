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

include { run_validate_PipeVal } from './external/pipeline-Nextflow-module/modules/PipeVal/validate/main.nf' addParams(
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
    Channel.from(params.samples_to_process)
        .map{ sample -> ['index': indexFile(sample.path)] + sample }
        .set{ input_ch_samples_with_index }

    input_ch_samples_with_index
        .map{ sample -> [sample.path, sample.index] }
        .flatten()
        .set{ input_ch_validate }

    input_ch_sample_with_index
        .map{ sample -> sample.id }
        .flatten()
        .set{ input_ch_sample_ids }

    input_ch_samples_with_index
        .map{ sample -> [sample.path, sample.index] }
        .flatten()
        .collect()
        .set{ input_ch_all_bams_and_indices }

    run_validate_PipeVal(input_ch_validate)

    run_validate_PipeVal.out.validation_result
        .collectFile(
            name: 'input_validation.txt',
            storeDir: "${params.output_dir_base}/validation"
        )

    extract_GenomeIntervals("${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict")

    run_SplitIntervals_GATK(
        extract_GenomeIntervals.out.genomic_intervals,
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict"
    )

    run_SplitIntervals_GATK.out.interval_list
        .flatten()
        .map{ interval_path ->
            [
                'interval_id': file(interval_path).getName().replace('-scattered.interval_list', ''),
                'interval_path': interval_path
            ]
        }
        .set{ input_ch_intervals }


}
