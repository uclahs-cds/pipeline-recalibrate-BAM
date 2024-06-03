nextflow.enable.dsl=2

include { remove_intermediate_files } from '../external/pipeline-Nextflow-module/modules/common/intermediate_file_removal/main.nf' addParams(
    options: [
        save_intermediate_files: !params.metapipeline_delete_input_bams || params.save_intermediate_files,
        output_dir: params.output_dir_base,
        log_output_dir: "${params.log_output_dir}/process-log"
        ]
    )

String get_root_directory(String directory) {
    Path root_directory = new File(directory).toPath()

    // Keep traversing parent directories until only the root directory is left
    while (root_directory.getParent().getParent()) {
        root_directory = root_directory.getParent()
    }

    return root_directory
}

/*
    Nextflow process to check if given input can be deleted.

    input:
        file_to_check: path to file to be checked for deletion
*/
process check_deletion_status {
    container params.docker_image_pipeval
    containerOptions "--volume ${get_root_directory(params.metapipeline_final_output_dir)}:${get_root_directory(params.metapipeline_final_output_dir)}"

    ext log_dir_suffix: { "-${task.index}" }

    input:
    path(file_to_check)

    output:
    tuple path(file_to_check), env(DELETE_INPUT), emit: file_to_delete

    when:
    params.metapipeline_delete_input_bams

    script:
    """
    set -euo pipefail

    FILE_NAME=`basename ${file_to_check}`
    until [ -f ${params.metapipeline_final_output_dir}/\$FILE_NAME ]
    do
        sleep 30
    done

    EXPECTED_PATH_TO_CHECK="\$(readlink ${file_to_check})"
    FINAL_PATH_TO_CHECK="${params.metapipeline_final_output_dir}/\$FILE_NAME"

    if [[ "\$EXPECTED_PATH_TO_CHECK" == "\$FINAL_PATH_TO_CHECK" ]]
    then
        DELETE_INPUT='false'
    else
        EXPECTED_SIZE=`stat --printf="%s" \$EXPECTED_PATH_TO_CHECK`
        COPIED_SIZE=`stat --printf="%s" \$FINAL_PATH_TO_CHECK`

        while [ \$EXPECTED_SIZE != \$COPIED_SIZE ]
        do
            sleep 30
            COPIED_SIZE=`stat --printf="%s" \$FINAL_PATH_TO_CHECK`
        done
        DELETE_INPUT='true'
    fi
    """
}

/*
    Delete pipeline input files.
    Intended to delete pipeline input files when recalibrate-BAM is run within a metapipeline.
    Requires metapipeline_delete_input_bams and metapipeline_final_output_dir params to be set.
    Deletes the input files (that reside in /scratch) once they're been copied to the final metapipeline output destination.
*/
workflow delete_input {
    take:
    files_to_delete

    main:
    check_deletion_status(files_to_delete)

    check_deletion_status.out.file_to_delete
        .filter{ it[1] == 'true' }
        .map{ it -> it[0] }
        .set{ files_to_delete }

    remove_intermediate_files(
        files_to_delete,
        "ready_to_delete"
        )
}
