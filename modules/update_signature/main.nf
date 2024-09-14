process UPDATE_SIGNATURE {
    tag "$meta.id"

    input:
    tuple val(meta), path(json_signature)

    output:
    tuple val(meta), path("updated_signature.json"), emit: updated_signature
    
    script:
    """
    update_signature.py \\
        -r ${workflow.projectDir}/assets/retinoblastoma_signatures.json \\
        -i ${json_signature}        
    """

    stub:
    """
    touch updated_signature.json
    """
}
