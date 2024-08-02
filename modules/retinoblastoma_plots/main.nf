process RETINOBLASTOMA_PLOTS {
    tag "$meta.id"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/${meta.id}/retinoblastoma_plots", failOnError: false

    input:
    tuple val(meta), path(seurat_object)

    output:
    tuple val(meta), path("*.png"), emit: retinoblastoma_plots

    script:
    """
    retinoblastoma_plots.R \\
        ${seurat_object} \\
        ${meta.id}
    """

    stub:
    """
    touch cones_${meta.id}_vln_plot.png
    """
}
