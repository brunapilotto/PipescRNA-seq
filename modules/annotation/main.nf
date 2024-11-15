process ANNOTATION_PLOTS {
    tag "$meta.id"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/${meta.id}/annotation_plots", failOnError: false

    input:
    tuple val(meta), path(seurat_object)

    output:
    tuple val(meta), path("*.png") , emit: annotation_plots
    tuple val(meta), path("*.txt") , emit: annotation_metrics
    tuple val(meta), path("*.rds") , emit: annotation_object
    tuple val(meta), path("*.json"), emit: immuno_signature

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
    touch ${meta.id}_annotation.rds
    touch immuno_signature.json
    """
}
