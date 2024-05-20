SCRIPT_PATH = "/home/mbruna/PipescRNA-seq/bin/cluster.R"

process SEURAT {
    label 'process_low'

    conda "/home/mbruna/anaconda3/envs/seurat"

    input:
    path seurat_obj
    path outdir

    output:
    path "${outdir}/seurat/*.png", emit: seuratPlots

    script:

    """
    ${SCRIPT_PATH} \\
        ${seurat_obj} \\
        ${outdir}/seurat \\
    """
}
