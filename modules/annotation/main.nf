process ANNOTATION_PLOTS {
    tag "$meta.id"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/${meta.id}/annotation_plots", failOnError: false

    input:
    tuple val(meta), path(seurat_object)

    output:
    tuple val(meta), path("*.png"), emit: annotation_plots
    tuple val(meta), path("*.txt"), emit: annotation_metrics

    script:
    """
    annotation_plots.R \\
        ${seurat_object} \\
        ${meta.id}
    """

    stub:
    """
    touch predictions_plot.png
    touch cell_type_counts.txt
    """
}
