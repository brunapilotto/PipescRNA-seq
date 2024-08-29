process MULTIQC {
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/multiqc", failOnError: false
    
    input:
    path(files)

    output:
    path("multiqc_data")       , emit: multiqc_dir
    path("multiqc_report.html"), emit: multiqc_report

    script:
    """
    multiqc --clean-up .
    """

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html
    """
}
