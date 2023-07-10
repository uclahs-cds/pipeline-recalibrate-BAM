# recalibrate-BAM

1. [Overview](#Overview)
2. [How To Run](#How-To-Run)
3. [Flow Diagram](#Flow-Diagram)
4. [Pipeline Steps](#Pipeline-Steps)
5. [Inputs](#Inputs)
5. [Outputs](#Outputs)
6. [Benchmarking](#Benchmarking)
7. [References](#References)

## Overview

This pipeline takes BAMs and corresponding indices from [pipeline-align-DNA](https://github.com/uclahs-cds/pipeline-align-DNA) and performs indel realignment and BQSR. It can be run on a single normal sample, a single tumor sample, a normal-tumor paired sample, and a 1 normal-n tumor set of samples.

---

## How To Run

**The pipeline is currently configured to run on a SINGLE NODE mode with normal only, tumour only, normal-tumour paired, or single normal-multiple tumour samples.**

1. Update the params section of the .config file ([Example config](config/template.config)).

2. Update the YAML.

3. Download the submission script (submit_nextflow_pipeline.py) from [here](https://github.com/uclahs-cds/tool-submit-nf), and submit your pipeline below.

> **Note**: Because this pipeline uses an image stored in the GitHub Container Registry, you must follow the steps listed in the [Docker Introduction](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3223396/Container+Registry+-+GitHub+Packages) on Confluence to set up a PAT for your GitHub account and log into the registry on the cluster before running this pipeline.

- YAML input
```
python submit_nextflow_pipeline.py \
       --nextflow_script /path/to/main.nf \
       --nextflow_config /path/to/call-gSNP.config \
       --nextflow_yaml /path/to/sample.yaml \
       --pipeline_run_name job_name \
       --partition_type <type> \
       --email email_address
```
---

## Flow Diagram

![call-gSNP flow diagram](call-gSNP-DSL2.png)

---

## Pipeline Steps

### 1. Split genome into intervals for parallelization
Use the input target intervals and split them into intervals for parallel processing.

### 2. Realign Indels
Generate indel realignment targets and realign indels.

### 3. Generate BQSR (Base Quality Score Recalibration)
Assess how sequencing errors correlate with four covariates (assigned quality score, read group the read belongs, machine cycle producing this base, and current and immediately upstream base), and output base quality score recalibration table.

### 4. Apply BQSR per split interval in parallel
Print out interval-level recalibrated BAM.

### 5. Reheader interval-level BAMs
In paired mode, reheader the interval-level BAMs.

### 6. Index reheadered BAMs
Index each reheadered interval-level BAM. After this step, the workflow splits into two: one path (7-10) merges the BAMs for Depth of Coverage and contamination calculations while the other path proceeds with the HaplotypeCaller (11-17).

### 7. Merge interval-level BAMs
Merge BAMs from each interval to generate whole sample BAM.

### 8. Get pileup summaries
Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results will be used in the next Calculate Contamination step.

### 9. Calculate contamination
Calculates the fraction of reads coming from cross-sample contamination, given results from Step 8. Generates a tumor segmentation file.

### 10.	DepthOfCoverage
Calculate depth of coverage using the whole sample BAM from step 7.

### 11.	HC – call raw VCF on each interval in parallel
Generate raw VCF for each split interval using HaplotypeCaller. Generate GVCF for SNPs and INDELs.

### 12. Merge raw VCFs
Merge raw variants from each interval to generate whole sample raw VCF.

### 13. VQSR (Variant Quality Score Recalibration): Generate VQSR model for SNPs.

### 14. VQSR: Generate VQSR model for INDELs.

### 15. VQSR: Take the whole sample raw VCF from Step 11 as input, and apply the model in Step 13 to generate variants in which only SNPs are recalibrated.

### 16. VQSR: Take the output from Step 15 as input, and apply the model in Step 14 to recalibrate only INDELs.

#### Steps 13 through 16 model the technical profile of variants in a training set and uses that to filter out probable artifacts from the raw VCF. After these four steps, a recalibrated VCF is generated.

### 17. Filter gSNP – Filter out ambiguous variants
Use customized Perl script to filter out ambiguous variants.

### 18. Generate sha512 checksum
Generate sha512 checksum for final BAM, filtered VCF, and GVCFs for SNPs and INDELs.

---

## Inputs

### Input CSV

| Field | Type | Description |
|:------|:-----|:------------|
| patient_id | string | Patient ID (will be standardized according to data storage structure in the near future) |
| sample_id | string | Sample ID |
| normal_id | string | Must be strictly set to the sample tag (SM:) in the BAM header @RG line (should be also in the pipeline-align-DNA input .csv file) |
| normal_BAM | path | Set to absolute path to normal BAM |
| tumour_id | string | Must be strictly set to the sample tag (SM:) in the BAM header @RG line (should be also in the pipeline-align-DNA input .csv file) |
| tumour_BAM | path | Set to absolute path to tumour BAM |

For normal-only or tumour-only samples, exclude the fields for the other state.

For inputs with one normal sample and multiple tumour samples, add rows. Keep the non-tumour related fields identical for each row and update the tumour fields.

### Config

| Input Parameter | Required | Type | Description |
|:----------------|:---------|:-----|:------------|
| `dataset_id` | Yes | string | Dataset ID |
| `avere_prefix` | Yes | string | Prefix for location of avere cache |
| `blcds_registered_dataset` | Yes | boolean | Set to true when using BLCDS folder structure; use false for now |
| `output_dir` | Yes | string | Need to set if `blcds_registered_dataset = false` |
| `input_csv` | Yes | path | Absolute path to input CSV file |
| `save_intermediate_files` | Yes | boolean | Set to false to disable publishing of intermediate files; true otherwise; disabling option will delete intermediate files to allow for processing of large BAMs |
| `aligner` | Yes | string | Original aligner used to align input BAMs; formatted as \<aligner\>-\<aligner-version\> |
| `cache_intermediate_pipeline_steps` | No | boolean | Set to true to enable process caching from Nextflow; defaults to false |
| `is_emit_original_quals` | Yes | boolean | Set to true to emit original quality scores; false to omit |
| `is_NT_paired` | Yes | boolean | Set to true for normal-tumour paired mode, and to false for normal only mode |
| `is_DOC_run` | Yes | boolean | Set to true to run GATK DepthOfCoverage (very time-consuming for large BAMs); false otherwise |
| `scatter_count` | Yes | integer | Number of intervals to divide into for parallelization |
| `intervals` | Yes | path | Use all .list in inputs for WGS; Set to absolute path to targeted exome interval file (with .interval_list, .list, .intervals, or .bed suffix) |
| `gatk_ir_compression` | No | integer | Compression level for BAMs output by IndelRealigner. Default: 0. Range: 0-9 |
| `reference_fasta` | Yes | path | Absolute path to reference genome fasta file, e.g., `/hot/ref/reference/GRCh38-BI-20160721/Homo_sapiens_assembly38.fasta` |
| `bundle_mills_and_1000g_gold_standard_indels_vcf_gz` | Yes | path | Absolute path to Mills & 1000G Gold Standard Indels file, e.g., `/hot/ref/tool-specific-input/GATK/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` |
| `bundle_known_indels_vcf_gz` | Yes | path | Absolute path to known indels file, e.g., `/hot/ref/tool-specific-input/GATK/GRCh38/Homo_sapiens_assembly38.known_indels.vcf.gz` |
| `bundle_v0_dbsnp138_vcf_gz` | Yes | path | Absolute path to dbsnp file, e.g., `/hot/ref/tool-specific-input/GATK/GRCh38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz` |
| `bundle_hapmap_3p3_vcf_gz` | Yes | path | Absolute path to HapMap 3.3 file, e.g., `/hot/ref/tool-specific-input/GATK/GRCh38/hapmap_3.3.hg38.vcf.gz` |
| `bundle_omni_1000g_2p5_vcf_gz` | Yes | path | Absolute path to 1000 genomes OMNI 2.5 file, e.g., `/hot/ref/tool-specific-input/GATK/GRCh38/1000G_omni2.5.hg38.vcf.gz` |
| `bundle_phase1_1000g_snps_high_conf_vcf_gz` | Yes | path | Absolute path to 1000 genomes phase 1 high-confidence file, e.g., `/hot/ref/tool-specific-input/GATK/GRCh38/1000G_phase1.snps.high_confidence.hg38.vcf.gz` |
| `bundle_contest_hapmap_3p3_vcf_gz` | Yes | path | Absolute path to HapMap 3.3 biallelic sites file, e.g., `/hot/ref/tool-specific-input/GATK/GRCh38/Biallelic/hapmap_3.3.hg38.BIALLELIC.PASS.2021-09-01.vcf.gz` |
| `work_dir` | optional | path | Path of working directory for Nextflow. When included in the sample config file, Nextflow intermediate files and logs will be saved to this directory. With ucla_cds, the default is `/scratch` and should only be changed for testing/development. Changing this directory to `/hot` or `/tmp` can lead to high server latency and potential disk space limitations, respectively. |
| `docker_container_registry` | optional | string | Registry containing tool Docker images. Default: `ghcr.io/uclahs-cds` |
| `metapipeline_delete_input_bams` | optional | boolean | Set to true to delete the input BAM files once the initial processing step is complete. **WARNING**: This option should NOT be used for individual runs of call-gSNP; it's intended for metapipeline-DNA to optimize disk space usage by removing files that are no longer needed from the `workDir`. |
| `metapipeline_final_output_dir` | optional | string | Absolute path for the final output directory of metapipeline-DNA that's expected to contain the output BAM from align-DNA. **WARNING**: This option should not be used for individual runs of call-gSNP; it's intended for metapipeline-DNA to optimize disk space usage. |

---

## Outputs

| Output | Description |
|:-------|:------------|
| `<aligner>_<GATK>_<dataset_id>_<normal_id>.bam` | Post-processed normal BAM |
| `<aligner>_<GATK>_<dataset_id>_<normal_id>.bam.bai` | Post-processed normal BAM index |
| `<aligner>_<GATK>_<dataset_id>_<normal_id>.bam.sha512` | Post-processed normal BAM sha512 checksum |
| `<aligner>_<GATK>_<dataset_id>_<tumor_id>.bam` | Post-processed tumour BAM if in normal-tumour paired mode |
| `<aligner>_<GATK>_<dataset_id>_<tumor_id>.bam.bai` | Post-processed tumour BAM index if in normal-tumour paired mode |
| `<aligner>_<GATK>_<dataset_id>_<tumor_id>.bam.sha512` | Post-processed tumour BAM sha512 checksum if in normal-tumour paired mode |
| `<GATK>_<dataset_id>_<normal_id\|tumor_id>.g.vcf.gz` | Per-sample GVCF |
| `<GATK>_<dataset_id>_<normal_id\|tumor_id>.g.vcf.gz.sha512` | Per-sample GVCF checksum |
| `<GATK>_<dataset_id>_<normal_id\|tumor_id>.g.vcf.gz.tbi` | Per-sample GVCF index |
| `<GATK>_<dataset_id>_<normal_id\|tumor_id>.g.vcf.gz.tbi.sha512` | Per-sample GVCF index checksum |
| `<GATK>_<dataset_id>_<sample_id>.vcf` | Raw variant calls |
| `<GATK>_<dataset_id>_<sample_id>.vcf.idx` | Raw variant calls index |
| `<GATK>_<dataset_id>_<sample_id>_VQSR-SNP-AND_INDEL.vcf.gz` | SNP and INDEL recalibrated variants |
| `<GATK>_<dataset_id>_<sample_id>_VQSR-SNP-AND_INDEL.vcf.gz.sha512` | SNP and INDEL recalibrated variants checksum |
| `<GATK>_<dataset_id>_<sample_id>_VQSR-SNP-AND_INDEL.vcf.gz.tbi` | SNP and INDEL recalibrated variants index |
| `<GATK>_<dataset_id>_<sample_id>_VQSR-SNP-AND_INDEL.vcf.gz.tbi.sha512` | SNP and INDEL recalibrated variants index checksum |
| `<GATK>_<dataset_id>_<sample_id>_snv.vcf.gz` | Filtered SNVs with non-germline and ambiguous variants removed |
| `<GATK>_<dataset_id>_<sample_id>_snv.vcf.gz.tbi` | Filtered germline SNVs index |
| `<GATK>_<dataset_id>_<sample_id>_snv.vcf.gz.sha512` | Filtered germline SNVs sha512 checksum |
| `<GATK>_<dataset_id>_<sample_id>_indel.vcf.gz` | Filtered INDELs with non-germline and ambiguous variants removed |
| `<GATK>_<dataset_id>_<sample_id>_indel.vcf.gz.tbi` | Filtered germline INDELs index |
| `<GATK>_<dataset_id>_<sample_id>_indel.vcf.gz.sha512` | Filtered germline INDELs sha512 checksum |
| `report.html`, `timeline.html` and `trace.txt` | Nextflow report, timeline and trace files |
| `*.command.*` | Process specific logging files created by nextflow |

---

## Benchmarking

### Test Data Set

1. A-mini: A subset dataset consisting of read alignment from chr8, chr21, and chrX.

### Results

---

## References

--

## License

Authors: Yash Patel (YashPatel@mednet.ucla.edu), Shu Tao (shutao@mednet.ucla.edu), Stefan Eng (stefaneng@mednet.ucla.edu)

Call-gSNP is licensed under the GNU General Public License version 2. See the file LICENSE for the terms of the GNU GPL license.

Call-gSNP takes BAM files and utilizes GATK to call short germline variants (SNP and INDEL).

Copyright (C) 2021-2023 University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
