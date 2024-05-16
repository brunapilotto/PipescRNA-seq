process PIGZ {
    label 'process_low'

    conda "/home/mbruna/PipescRNA-seq/envs/pigz.yml"

    input:
    path raw_files

    output:
    path("${raw_files}/matrix.mtx.gz") , emit: matrix_zip
    path("${raw_files}/barcodes.tsv.gz") , emit: barcodes_zip
    path("${raw_files}/features.tsv.gz") , emit: features_zip

    script:
    """
    pigz -f -p $task.cpus -k ${raw_files}/matrix.mtx > ${raw_files}/matrix.mtx.gz
    pigz -f -p $task.cpus -k ${raw_files}/barcodes.tsv > ${raw_files}/barcodes.tsv.gz
    pigz -f -p $task.cpus -k ${raw_files}/features.tsv > ${raw_files}/features.tsv.gz
    """
}
