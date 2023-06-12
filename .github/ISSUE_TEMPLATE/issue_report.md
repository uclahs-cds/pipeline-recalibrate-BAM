---
name: Issue report
about: Create a report to help us improve our pipelines
title: "[ISSUE]: "
labels: ''
assignees: ''

---

**Describe the issue**<br>
A clear and concise description of what the issue is. Please include the following in your issue report along with any explicit errors observed

**Error:**
```
paste relevant error message here, if applicable
```

**Screenshots**<br>
If applicable, add screenshots to help explain your problem.

---

**Pipeline release version:**<br>
`version`

**Cluster you are using:**
- [ ] `Slurm-Dev`
- [ ] `Slurm-Test`
- [ ] `Regeneron`

**Node type:**
<!-- specify another scratch size here if it's not the default. eg.) 3TB for F72+3TB  --->
- [ ] `F2`   &numsp;&emsp;  (scratch: `32GB`)
- [ ] `F16`         &emsp;  (scratch: `500GB`)
- [ ] `F32`         &emsp;  (scratch: `500GB`)
- [ ] `F72`         &emsp;  (scratch: `2TB`)
- [ ] `M64`         &emsp;  (scratch: `2TB`)

**Submission method**
- [ ] interactive
- [ ] submission script

**Actual submission command**<br>
Paste the command used inside your script if you used one.
<!-- 
EXAMPLE:
    - `python3 /path/to/submit_nextflow_pipeline.py --nextflow_script ...`
    - `nextflow run ...`
    - etc...
--->

```
submission command
```

**Paths**
<!-- Please enter the paths under the table below --->

| Description | Path |
|-------------|------|
| sbatch/qsub command and logs <br> (if applicable) | ```/path/to/logs```|
| Config File                                       | ```/path/to/config/file``` |
| Working Directory                                 | ```/path/to/work/dir``` |
| Logs Produceded by the Pipeline                   | ```/path/to/logs``` |


---

**To Reproduce**<br>
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**<br>
A clear and concise description of what you expected to happen.

**Additional context**<br>
Add any other context about the problem here.
