process SINGLE_SEURAT {
    tag "$meta.id"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/${meta.id}/single_seurat_analysis", failOnError: false

    input:
    tuple val(meta), path(filtered_clean_matrix)

    output:
    tuple val(meta), path("*.png"), emit: single_seurat_plots
    tuple val(meta), path("*.rds"), emit: single_seurat_object

    script:
    """
    cluster.R \\
        ${filtered_clean_matrix}
    """

    stub:
    """
    touch heatmap.png
    touch single_final_cluster.rds
    """
}
