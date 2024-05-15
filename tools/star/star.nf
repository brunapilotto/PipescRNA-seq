STAR_PATH = "/home/mbruna/STAR-2.7.11b/bin/Linux_x86_64_static/STAR"
INDEX_PATH = "/data1/projects/LGMB-018/reference_cell_ranger/star_index"
GTF_PATH = "/data1/projects/LGMB-018/reference_cell_ranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
WHITE_LIST_PATH = "/home/mbruna/cellranger-8.0.0/lib/python/cellranger/barcodes/737K-august-2016.txt"

process STAR {
    
    input:
    path read_1
    path read_2
    val sample_name
    path outdir

    output:
    path("${outdir}/star/*Log.final.out")                         , emit: log_final
    path("${outdir}/star/*Aligned.out.bam")                       , emit: bam
    path("${outdir}/star/*Solo.out")                              , emit: counts
    path("${outdir}/star/*Log.out")                               , emit: log_out
    path("${outdir}/star/*Log.progress.out")                      , emit: log_progress
    path("${outdir}/star/*_STARgenome")                           , emit: star_genome
    path("${outdir}/star/*Solo.out/Gene/raw")                     , emit: raw_counts
    path("${outdir}/star/*Solo.out/Gene/filtered")                , emit: filtered_counts
    path("${outdir}/star/*Aligned.sortedByCoord.out.bam")         , emit: bam_sorted
    path("${outdir}/star/*.tab")                                  , emit: tab
    path("${outdir}/star/*Solo.out/Gene/Summary.csv")             , emit: summary_csv

    
    shell:
    """
    if [ -z "\$(find "${outdir}/star" -mindepth 1 -print -quit)" ]; then
        mkdir -p ${outdir}/star
        ${STAR_PATH} --genomeDir ${INDEX_PATH} \\
        --readFilesIn ${read_2} ${read_1} \\
        --runThreadN 12 \\
        --soloType CB_UMI_Simple \\
        --outFileNamePrefix ${outdir}/star/${sample_name} \\
        --soloCBwhitelist ${WHITE_LIST_PATH} \\
        --soloFeatures Gene \\
        --soloUMIlen 10 \\
        --sjdbGTFfile ${GTF_PATH} \\
        --outWigType bedGraph \\
        --twopassMode Basic \\
        --outSAMtype BAM Unsorted SortedByCoordinate \\
        --runDirPerm All_RWX \\
        --soloBarcodeReadLength 150 \\
        --soloCellFilter EmptyDrops_CR
    fi
    """
}