SCRIPT_PATH = "/home/mbruna/PipescRNA-seq/bin/emptydrops_cell_calling.R"

process EMPTYDROPS_CELL_CALLING {
    label 'process_medium'

    conda "bioconda::bioconductor-dropletutils"

    input:
    path matrix
    path barcodes
    path features
    path outdir

    output:
    path("${outdir}/emptydrops_filtered"), emit: filtered_matrices

    script:
    """
    if [ -z "\$(find "${outdir}/emptydrops_filtered" -mindepth 1 -print -quit)" ]; then
        mkdir -p ${outdir}/emptydrops_filtered
        ${SCRIPT_PATH} \\
            ${matrix} \\
            ${barcodes} \\
            ${features} \\
            ${outdir}/emptydrops_filtered \\
            star \\
            0
    fi
    """
}
