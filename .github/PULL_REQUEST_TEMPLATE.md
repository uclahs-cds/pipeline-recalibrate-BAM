# Description
<!-- Briefly describe the changes included in this pull request and the paths
     to the test cases below starting with 'Closes #...' if appropriate -->

### Closes #...

## Testing Results

<!-- If you did not run NFTest, please justify _why_ you feel your changes
     don't need to be tested. -->

- NFTest
  - log:      /hot/software/pipeline/pipeline-recalibrate-BAM/Nextflow/development/unreleased/\<branchname\>/log-nftest-\<datestamp\>.log
  - cases:    default set <!-- update this if you made any changes to nftest.yml or ran default disabled test cases explicitly -->

<!--
    Include any non-NFTest test details here, using the "Additional Case"
    template below as appropriate. Please make sure that a reviewer can
    understand:
    * How you tested the pipeline
    * How others can test it (including paths to config and YAML files)
    * How the results demonstrate the changes in this PR. Examples of this
      might include:
        * The run completed successfully, so nothing was broken
        * An output file changed (this might require two pipeline runs, before
          and after your changes, to demonstrate the difference)
    -->
- Additional Case 1
  - sample:    <!-- e.g. A-mini S2.T-1, A-mini S2.T-n1 -->
  - input csv: <!-- path/to/input.csv -->
  - config:    <!-- path/to/xxx.config -->
  - output:    <!-- path/to/output -->

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
