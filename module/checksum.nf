/*
  Nextflow module for calculating SHA512 checksum

  input:
      file_for_calc: path to file for which to calculate checksum

  params:
      params.output_dir_base: string(path)
      params.log_output_dir: string(path)
      params.docker_image_pipeval: string
*/
process calculate_sha512 {
    container params.docker_image_pipeval
    publishDir path: "${params.output_dir_base}/output",
      mode: "copy",
      pattern: "*.sha512",
      saveAs: { filename -> (filename.endsWith(".bai.sha512") && !filename.endsWith(".bam.bai.sha512")) ? "${file(file(filename).baseName).baseName}.bam.bai.sha512" : "${filename}"}

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}/${task.process.split(':')[-1]}-${task.index}/log${file(it).getName()}" }

    input:
    path(file_for_calc)

    output:
    path(".command.*")
    path("*.sha512"), emit: sha512_sum

    script:
    """
    set -euo pipefail
    sha512sum ${file_for_calc} > ${file_for_calc}.sha512
    """
}
