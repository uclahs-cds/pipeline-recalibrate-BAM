#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow {
    System.out.println(params)
    System.out.println(params.samples_to_process)
}
