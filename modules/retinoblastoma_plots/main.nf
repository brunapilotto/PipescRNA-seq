process RETINOBLASTOMA_PLOTS {
    tag "$meta.id"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/${meta.id}/retinoblastoma_plots", failOnError: false

    input:
    tuple val(meta), path(seurat_object), path(signature)

    output:
    tuple val(meta), path("*.png"), emit: retinoblastoma_plots
    tuple val(meta), path("*.txt"), emit: retinoblastoma_metrics
    tuple val(meta), path("*.pdf"), emit: auc_plots

    script:
    """
    retinoblastoma_plots.R \\
        ${seurat_object} \\
        ${meta.id} \\
        ${signature}
    """

    stub:
    """
    touch cones_${meta.id}_vln_plot.png
    touch cell_type_counts.txt
    touch Rplots.pdf
    """
}
