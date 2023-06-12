/*
*   Module/process description here
*
*   @input  <name>  <type>  <description>
*   @params <name>  <type>  <description>
*   @output <name>  <type>  <description>
*/
process tool_name_command_name {
    container params.docker_image_name

    label "resource_allocation_tool_name_command_name"

    publishDir path: "${params.output_dir}/${task.process.replace(':','/')}-${task.index}",
        pattern: "<file name pattern>",
        mode: "copy",
        enabled: true

    publishDir path: "${params.log_output_dir}",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    // Additional directives here
    
    input: 
        tuple val(row_1_name), path(row_2_name_file_extension)
        val(variable_name)

    output:
        path("${variable_name}.command_name.file_extension"), emit: output_ch_tool_name_command_name

    script:
    """
    # make sure to specify pipefail to make sure process correctly fails on error
    set -euo pipefail

    # the script should ideally only have call to a tool
    # to make the command more human readable:
    #  - seperate components of the call out on different lines
    #  - when possible by explict with command options, spelling out their long names
    tool_name \
        command_name \
        --option_1_long_name ${row_1_name} \
        --input ${row_2_name_file_extension} \
        --output ${variable_name}.command_name.file_extension
    """
}
