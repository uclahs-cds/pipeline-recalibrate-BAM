# Changelog
All notable changes to the pipeline-recalibrate-BAM pipeline.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### [Changed]
- Use `run_validate_PipeVal_with_metadata` to gate on validation
- Move index/dictionary file discovery to configuration stage
- Include all parameters from `default.config` in the README

---

## [1.0.1] - 2024-05-29

### [Changed]
- Update Nextflow configuration test workflows
- Use `methods.setup_process_afterscript()` for process logs

### [Fixed]
- Resolve interval path for original given intervals to ensure the index file for intervals is found when present

---

## [1.0.0] - 2024-03-13
### [Added]
- Sort BAMs before merging for consistent order of output BAM @PG header lines
- Add NFTest case
- Add new flow diagram to README
- Add additional details to Pipeline Steps section of README
- Option to provide base recalibration tables for any subset of samples to skip `BaseRecalibrator`
- Output pipeline parameters to log directory using `store_object_as_json`
- Add Action to deploy documentation to GitHub Pages
- Add Nextflow configuration test action and two regression tests
- Add requirements to README

### [Changed]
- Use modularized `set_env` function
- Use modularized `check_limits` function

---

## [1.0.0-rc.4] - 2023-09-29
### [Fixed]
- Resource updater to allow update for all processes
- Channel handling to save contamintion estimate files

---

## [1.0.0-rc.3] - 2023-08-09
### [Added]
- Customization options for which input BAMs to delete with metapipeline

---

## [1.0.0-rc.2] - 2023-08-02
### [Added]
- Custom resource allocation updates through configuration parameters

---

## [1.0.0-rc.1] - 2023-07-27
### [Added]
- Basic configuration setup
- Input YAML templates and handling
- Per-patient input YAML handling
- Interval splitting process, either per-chromosome by default or scattered
- Indel Realignment workflow
- BQSR workflow
- Merging workflow
- Contamination and depth of coverage processes
- Input deletion
- Intermediate file deletion
- Schema validation
- Documentation
