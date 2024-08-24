process PIGZ {
    tag "$meta.id"
    conda "${moduleDir}/environment.yml"
    
    input:
    tuple val(meta), val(reads), val(white_list)

    output:
    tuple val(meta), val(reads), path("10x_V*_barcode_whitelist.txt"), emit: white_list

    script:
    """
    pigz \\
        -d -k -c \\
        -p $task.cpus \\
        ${workflow.projectDir}/assets/10x_${white_list}_barcode_whitelist.txt.gz \\
        > 10x_${white_list}_barcode_whitelist.txt
    """

    stub:
    """
    touch 10x_${white_list}_barcode_whitelist.txt
    """
}
