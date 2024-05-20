SCRIPT_PATH = "/home/mbruna/PipescRNA-seq/bin/mtx_to_seurat.R"

process MTX_TO_SEURAT {
    label 'process_single'

    conda "/home/mbruna/anaconda3/envs/seurat"

    input:
    path input
    path outdir
    val sample_name

    output:
    path "${outdir}/seurat/*.rds", emit: seuratObject

    script:

    matrix   = "${input}/matrix.mtx"
    barcodes = "${input}/barcodes.tsv"
    features = "${input}/features.tsv"


    """
    mkdir -p ${outdir}/seurat
    ${SCRIPT_PATH} \\
        $matrix \\
        $barcodes \\
        $features \\
        ${outdir}/seurat/${sample_name}_filtered_matrix.rds \\
    """
}
