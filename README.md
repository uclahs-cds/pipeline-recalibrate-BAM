# Pipeline Name

- [Pipeline Name](#pipeline-name)
  - [Overview](#overview)
  - [How To Run](#how-to-run)
  - [Flow Diagram](#flow-diagram)
  - [Pipeline Steps](#pipeline-steps)
    - [1. Step/Proccess 1](#1-stepproccess-1)
    - [2. Step/Proccess 2](#2-stepproccess-2)
    - [3. Step/Proccess n](#3-stepproccess-n)
  - [Inputs](#inputs)
  - [Outputs](#outputs)
  - [Testing and Validation](#testing-and-validation)
    - [Test Data Set](#test-data-set)
    - [Validation <version number\>](#validation-version-number)
    - [Validation Tool](#validation-tool)
  - [References](#references)
  - [Discussions](#discussions)
  - [Contributors](#contributors)
  - [License](#license) 
## Overview

A 3-4 sentence summary of the pipeline, including the pipeline's purpose, the type of expected scientific inputs/outputs of the pipeline (e.g: FASTQs and BAMs), and a list of tools/steps in the pipeline.

---

## How To Run

1. Update the params section of the .config file

2. Update the input yaml

3. See the submission script, [here](https://github.com/uclahs-cds/tool-submit-nf), to submit your pipeline

---

## Flow Diagram

A directed acyclic graph of your pipeline.

![alt text](pipeline-name-DAG.png?raw=true)

---

## Pipeline Steps

### 1. Step/Proccess 1

> A 2-3 sentence description of each step/proccess in your pipeline that includes the purpose of the step/process, the tool(s) being used and their version, and the expected scientific inputs/outputs (e.g: FASTQs and BAMs) of the pipeline.

### 2. Step/Proccess 2

> A 2-3 sentence description of each step/proccess in your pipeline that includes the purpose of the step/process, the tool(s) being used and their version, and the expected scientific inputs/outputs (e.g: FASTQs and BAMs) of the pipeline.

### 3. Step/Proccess n

> A 2-3 sentence description of each step/proccess in your pipeline that includes the purpose of the step/process, the tool(s) being used and their version, and the expected scientific inputs/outputs (e.g: FASTQs and BAMs) of the pipeline.

---

## Inputs

### Input YAML

> include an example of the organization structure within the YAML. Example:
```yaml
input 1: 'patient_id'
input:
    normal:
      - id: <normal id>
        BAM: </path/to/normal.bam>
    tumor:
      - id: <tumor id>
        BAM: </path/to/tumor.bam>
```

### Config

| Field | Type | Required | Description |
| ----- | ---- | ------------ | ------------------------ |
| param 1 | _type_ | yes/no | 1-2 sentence description of the parameter, including any defaults if any. |
| param 2 | _type_ | yes/no | 1-2 sentence description of the parameter, including any defaults if any. |
| param n | _type_ | yes/no | 1-2 sentence description of the parameter, including any defaults if any. |
| `work_dir` | path | no | Path of working directory for Nextflow. When included in the sample config file, Nextflow intermediate files and logs will be saved to this directory. With ucla_cds, the default is `/scratch` and should only be changed for testing/development. Changing this directory to `/hot` or `/tmp` can lead to high server latency and potential disk space limitations, respectively. |

> Include the optional param `work_dir` in the inputs accompanied by a warning of the potentials dangers of using the param. Update the warning if necessary.

---

## Outputs

<!-- List and describe the final output(s) of the pipeline. -->

| Output | Description |
| ------------ | ------------------------ |
| ouput 1 | 1 - 2 sentence description of the output. |
| ouput 2 | 1 - 2 sentence description of the output. |
| ouput n | 1 - 2 sentence description of the output. |

---

## Testing and Validation

### Test Data Set

A 2-3 sentence description of the test data set(s) used to validate and test this pipeline. If possible, include references and links for how to access and use the test dataset

### Validation <version number\>

 Input/Output | Description | Result  
 | ------------ | ------------------------ | ------------------------ |
| metric 1 | 1 - 2 sentence description of the metric | quantifiable result |
| metric 2 | 1 - 2 sentence description of the metric | quantifiable result |
| metric n | 1 - 2 sentence description of the metric | quantifiable result |

- [Reference/External Link/Path 1 to any files/plots or other validation results](<link>)
- [Reference/External Link/Path 2 to any files/plots or other validation results](<link>)
- [Reference/External Link/Path n to any files/plots or other validation results](<link>)

### Validation Tool

Included is a template for validating your input files. For more information on the tool check out: https://github.com/uclahs-cds/tool-validate-nf

---

## References

1. [Reference 1](<links-to-papers/external-code/documentation/metadata/other-repos/or-anything-else>)
2. [Reference 2](<links-to-papers/external-code/documentation/metadata/other-repos/or-anything-else>)
3. [Reference n](<links-to-papers/external-code/documentation/metadata/other-repos/or-anything-else>)

---

## Discussions

- [Issue tracker](<link-to-repo-issues-page>) to report errors and enhancement ideas.
- Discussions can take place in [<pipeline> Discussions](<link-to-repo-discussions-page>)
- [<pipeline> pull requests](<link-to-repo-pull-requests>) are also open for discussion

---

## Contributors

> Update link to repo-specific URL for GitHub Insights Contributors page.

Please see list of [Contributors](https://github.com/uclahs-cds/template-NextflowPipeline/graphs/contributors) at GitHub.

---

## License

[pipeline name] is licensed under the GNU General Public License version 2. See the file LICENSE for the terms of the GNU GPL license.

<one line to give the program's name and a brief idea of what it does.>

Copyright (C) 2023 University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
