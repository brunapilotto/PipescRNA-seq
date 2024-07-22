process QC_SEURAT {
    tag "$meta.id"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/${meta.id}/qc", failOnError: false

    input:
    tuple val(meta), path(filtered_matrix)

    output:
    tuple val(meta), path("*.png"), emit: qc_plots
    tuple val(meta), path("*.rds"), emit: clean_object
    tuple val(meta), path("*.txt"), emit: qc_metrics

    script:
    """
    QC.R \\
        ${filtered_matrix} \\
        ${meta.id}
    """

    stub:
    """
    touch qc_before.png
    touch qc_after.png
    touch feature-feature.png
    touch clean_object.rds
    touch qc_metrics_before.txt
    touch qc_metrics_after.txt
    """
}
