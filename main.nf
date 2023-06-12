#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include processes and workflows here
include { run_validate_PipeVal } from '.external/pipeline-Nextflow-module/modules/PipeVal/validate/main.nf' addParams(
    options: [
        docker_image_version: params.pipeval_version,
        main_process: "./" //Save logs in <log_dir>/process-log/run_validate_PipeVal
    ]
)
include { tool_name_command_name } from './module/module-name'

// Log info here
log.info """\    
        ======================================
        T E M P L A T E - N F  P I P E L I N E
        ======================================
        Boutros Lab

        Current Configuration:
        - pipeline:
            name: ${workflow.manifest.name}
            version: ${workflow.manifest.version}

        - input:
            input a: ${params.variable_name}
            ...

        - output: 
            output a: ${params.output_path}
            ...

        - options:
            option a: ${params.option_name}
            ...

        Tools Used:
            tool a: ${params.docker_image_name}

        ------------------------------------
        Starting workflow...
        ------------------------------------
        """
        .stripIndent()

// Channels here
// Decription of input channel
Channel
    .from( params.input['tumor'] )
    .multiMap{ it ->
        tumor_bam: it['BAM']
        tumor_index: indexFile(it['BAM'])
    }
    .set { tumor_input }

Channel
    .from( params.input['normal'] )
    .multiMap{ it ->
        normal_bam: it['BAM']
        normal_index: indexFile(it['BAM'])
    }
    .set { normal_input }

// Decription of input channel
Channel
    .fromPath(params.variable_name)
    .ifEmpty { error "Cannot find: ${params.variable_name}" }
    .into { input_ch_variable_name } // copy into two channels, one is for validation

// Pre-validation steps
tumor_input
    .mix ( normal_input )
    .set { ich_validation }

// Main workflow here
workflow {
    // Validation process
    run_validate_PipeVal(
        ich_validation
        )

    run_validate_PipeVal.out.val_file.collectFile(
        name: 'input_validation.txt', newLine: true,
        storeDir: "${params.output_dir_base}/validation"
        )

    // Workflow or process
    tool_name_command_name(
        tumor_input.tumor_bam,
        tumor_input.tumor_index,
        normal_input.normal_bam,
        normal_input.normal_index,
        input_ch_variable_name
        )
}
