process INTEGRATION {
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/Integrated", failOnError: false

    input:
    path(objects)

    output:
    path("*.rds"), emit: merged_object
    path("*.png"), emit: merged_plots
    path("*.txt"), emit: merged_tables

    script:
    """
    integration.R \\
        $objects
    """

    stub:
    """
    touch merged_obj.rds
    touch merged_plots.png
    touch merged_tables.txt
    """
}
