process MTX_TO_SEURAT {
    tag "$meta.id"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/${meta.id}/mtx_seurat", failOnError: false

    input:
    tuple val(meta), path(filtered_dir)

    output:
    tuple val(meta), path("*_filtered_matrix.rds"), emit: seurat_object

    script:
    """
    mtx_to_seurat.R \\
        "${filtered_dir}/matrix.mtx" \\
        "${filtered_dir}/barcodes.tsv" \\
        "${filtered_dir}/features.tsv" \\
        ${meta.id}_filtered_matrix.rds \\
        ${meta.id}
    """

    stub:
    """
    touch ${meta.id}_filtered_matrix.rds
    """
}
