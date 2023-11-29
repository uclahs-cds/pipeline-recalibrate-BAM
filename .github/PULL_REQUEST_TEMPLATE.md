# Description
<!-- Briefly describe the changes included in this pull request and the paths
     to the test cases below starting with 'Closes #...' if appropriate -->

### Closes #...

## Testing

### NFTest Results

<!-- Please include the last 10 lines of your NFTest logfile in the code block
     below - you can collect those by modifying the `tail` command with the name
     of your branch and the relevant logfile.

     Include the full path to the NFTest log file.

     If you did not run NFTest, please justify _why_ you feel your changes
     don't need to be tested.
     -->

```console
$ tail /hot/software/pipeline/pipeline-recalibrate-BAM/Nextflow/development/unreleased/<branchname>/log-nftest-<datestamp>.log

... last 10 lines of NFTest output here ...
```

### Other Test Results

<!-- Include any non-NFTest test details here. If you have modified anything
     that is is not tested by NFTest, please detail:
       * How you tested it
       * How others can test it (including paths to config and YAML files)
       * Why you are confident in your test results
-->

# Checklist
<!-- Please read each of the following items and confirm by replacing the [ ] with a [X] -->

- [ ] I have read the [code review guidelines](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3187646/Code+Review+Guidelines) and the [code review best practice on GitHub check-list](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3189956/Code+Review+Best+Practice+on+GitHub+-+Check+List).

- [ ] I have reviewed the [Nextflow pipeline standards](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3193890/Nextflow+pipeline+standardization).

- [ ] The name of the branch is meaningful and well formatted following the [standards](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3189956/Code+Review+Best+Practice+on+GitHub+-+Check+List), using \[AD\_username (or 5 letters of AD if AD is too long)]-\[brief\_description\_of\_branch\].

- [ ] I have set up or verified the branch protection rule following the [github standards](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3190380/GitHub+Standards#GitHubStandards-Branchprotectionrule) before opening this pull request.

- [ ] I have added my name to the contributors listings in the ``manifest`` block in the `nextflow.config` as part of this pull request, am listed
already, or do not wish to be listed. (*This acknowledgement is optional.*)

- [ ] I have added the changes included in this pull request to the `CHANGELOG.md` under the next release version or unreleased, and updated the date.

- [ ] I have updated the version number in the `metadata.yaml` and `manifest` block of the `nextflow.config` file following [semver](https://semver.org/), or the version number has already been updated. (*Leave it unchecked if you are unsure about new version number and discuss it with the infrastructure team in this PR.*)

- [ ] I have tested the pipeline using NFTest, _or_ I have justified why I did not need to run NFTest above.
