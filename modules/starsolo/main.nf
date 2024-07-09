process STARSOLO {
    tag "$meta.id"
    maxForks 1
    debug true
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/${meta.id}/star", failOnError: false
    
    input:
    tuple val(meta), val(reads)

    output:
    tuple val(meta), path("*Log.final.out")         , emit: log_final
    tuple val(meta), path("*Log.out")               , emit: log_out
    tuple val(meta), path("*Log.progress.out")      , emit: log_progress
    tuple val(meta), path("*Solo.out/Gene/filtered"), emit: filtered_counts

    script:
    """
    STAR --genomeDir ${params.index_dir} \\
    --readFilesIn ${reads.fastq2} ${reads.fastq1} \\
    --runThreadN $task.cpus \\
    --outFileNamePrefix ${meta.id} \\
    --soloCBwhitelist ${params.white_list_path} \\
    --sjdbGTFfile ${params.gtf_path} \\
    --readFilesCommand zcat \\
    --soloType CB_UMI_Simple \\
    --soloFeatures Gene \\
    --soloUMIlen 10 \\
    --twopassMode Basic \\
    --outSAMtype None \\
    --runDirPerm All_RWX \\
    --soloBarcodeReadLength 150 \\
    --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \\
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
    --soloUMIfiltering MultiGeneUMI_CR \\
    --soloUMIdedup 1MM_CR
    """

    stub:
    """
    mkdir -p ${meta.id}Solo.out/Gene/filtered
    touch ${meta.id}Log.final.out
    touch ${meta.id}Log.out
    touch ${meta.id}Log.progress.out
    """
}