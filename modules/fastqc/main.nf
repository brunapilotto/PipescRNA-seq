process FASTQC {
    tag "$meta.id"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/${meta.id}/fastqc", failOnError: false
    
    input:
    tuple val(meta), val(reads), val(white_list)

    output:
    tuple val(meta), path("*.zip") , emit: fastqc_zip
    tuple val(meta), path("*.html"), emit: fastqc_html
    
    script:
    """
    fastqc -o ./ -f fastq ${reads.fastq1} ${reads.fastq2} -t $task.cpus
    """

    stub:
    """
    touch ${meta.id}.fastqc.zip
    touch ${meta.id}.fastqc.html
    """
}
