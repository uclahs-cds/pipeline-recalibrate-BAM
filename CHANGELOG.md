# Changelog
All notable changes to the pipeline-recalibrate-BAM pipeline.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]
### [Added]
- Sort BAMs before merging for consistent order of output BAM @PG header lines
- Add NFTest case

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
