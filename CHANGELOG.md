# Changelog
All notable changes to the pipeline-recalibrate-BAM pipeline.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

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
