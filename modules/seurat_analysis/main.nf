process SEURAT {
    tag "$meta.id"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/${meta.id}/seurat_plots", failOnError: false

    input:
    tuple val(meta), path(filtered_matrix)

    output:
    tuple val(meta), path("plots_*/*.png"), emit: seurat_plots
    tuple val(meta), path("plots_*/*.rds"), emit: cluster

    script:
    """
    mkdir plots_${meta.id}
    cluster.R \\
        ${filtered_matrix} \\
        plots_${meta.id}
    """

    stub:
    """
    mkdir plots_${meta.id}
    touch plots_${meta.id}/heatmap.png
    touch plots_${meta.id}/qc_before.png
    touch plots_${meta.id}/cluster_final.rds
    """
}
