STAR_PATH = "/home/mbruna/STAR-2.7.11b/bin/Linux_x86_64_static/STAR"
INDEX_PATH = "/data1/projects/LGMB-018/reference_cell_ranger/refdata-gex-GRCh38-2024-A/star/"
GTF_PATH = "/data1/projects/LGMB-018/reference_cell_ranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
WHITE_LIST_PATH = "/home/mbruna/PipescRNA-seq/assets/10x_V2_barcode_whitelist.txt"

process STAR {
    label 'process_high'
    
    input:
    path read_1
    path read_2
    val sample_name
    path outdir

    output:
    path("${outdir}/star/*Log.final.out")                , emit: log_final
    path("${outdir}/star/*Solo.out")                     , emit: counts
    path("${outdir}/star/*Log.out")                      , emit: log_out
    path("${outdir}/star/*Log.progress.out")             , emit: log_progress
    path("${outdir}/star/*_STARgenome")                  , emit: star_genome
    path("${outdir}/star/*Solo.out/Gene/raw")            , emit: raw_counts
    path("${outdir}/star/*Solo.out/Gene/filtered")       , emit: filtered_counts
    path("${outdir}/star/*Solo.out/Gene/Summary.csv")    , emit: summary_csv

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p ${outdir}/star
    ${STAR_PATH} --genomeDir ${INDEX_PATH} \\
    --readFilesIn ${read_2} ${read_1} \\
    --runThreadN $task.cpus \\
    --outFileNamePrefix ${outdir}/star/${sample_name} \\
    --soloCBwhitelist ${WHITE_LIST_PATH} \\
    --sjdbGTFfile ${GTF_PATH} \\
    $args
    """
}