STAR_PATH = "/home/mbruna/STAR-2.7.1a/bin/Linux_x86_64_static/STAR"
INDEX_PATH = "/data1/projects/LGMB-018/reference_cell_ranger/refdata-gex-GRCh38-2024-A/star/"
WHITE_LIST_PATH = "/home/mbruna/cellranger-8.0.0/lib/python/cellranger/barcodes/737K-august-2016.txt"

process STAR {
    
    input:
    path read_1
    path read_2
    val sample_name
    path outdir

    output:
    path "${outdir}/star/${sample_name}Log.final.out"
    
    shell:
    """
    if [ -z "\$(find "${outdir}/star" -mindepth 1 -print -quit)" ]; then
        mkdir -p ${outdir}/star
        ${STAR_PATH} --genomeDir ${INDEX_PATH} \\
        --readFilesIn ${read_2} ${read_1} \\
        --runThreadN 20 \\
        --soloType Droplet \\
        --outFileNamePrefix ${outdir}/star/${sample_name} \\
        --soloCBwhitelist ${WHITE_LIST_PATH} \\
        --soloBarcodeReadLength 150
    fi
    """
}