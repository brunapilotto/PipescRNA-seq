SCRIPT_PATH = "/home/mbruna/PipescRNA-seq/bin/mtx_to_h5ad.py"
INDEX_PATH = "/data1/projects/LGMB-018/reference_cell_ranger/star_index/"

process MTX_TO_H5AD {
    label 'process_low'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    
    input:
    path inputs
    path txp2gene
    val sample_name
    path outputdir

    output:
    tuple val(input_type), path("${outputdir}/h5ad/*h5ad") , emit: h5ad

    script:
    // Get a file to check input type. Some aligners bring arrays instead of a single file.
    def input_to_check = (inputs instanceof String) ? inputs : inputs[0]

    // check input type of inputs
    input_type = (input_to_check.toUriString().contains('unfiltered') || input_to_check.toUriString().contains('raw')) ? 'raw' : 'filtered'
    if (input_to_check.toUriString().contains('emptydrops')) { input_type = 'custom_emptydrops_filter' }


    mtx_dir      = (input_type == 'custom_emptydrops_filter') ? 'emptydrops_filtered' : "${input_type}"
    suffix       = (input_type == 'custom_emptydrops_filter') || (input_type == 'filtered') ? '' : '.gz'
    mtx_matrix   = "${mtx_dir}/matrix.mtx${suffix}"
    barcodes_tsv = "${mtx_dir}/barcodes.tsv${suffix}"
    features_tsv = "${mtx_dir}/features.tsv${suffix}"


    """
    mkdir -p ${outputdir}/h5ad
    ${SCRIPT_PATH} \\
    --task_process ${task.process} \\
    --aligner star \\
    --sample ${sample_name} \\
    --input $mtx_matrix \\
    --barcode $barcodes_tsv \\
    --feature $features_tsv \\
    --txp2gene ${txp2gene} \\
    --star_index ${INDEX_PATH} \\
    --out ${outputdir}/h5ad/${sample_name}_${input_type}_matrix.h5ad
    --verbose True
    """
}
