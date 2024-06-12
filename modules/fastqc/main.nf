process FASTQC {
    tag "$meta.id"
    debug true
    publishDir "${params.outdir}/${meta.id}/fastqc", failOnError: false
    
    input:
    tuple val(meta), path(read_1), path(read_2)

    output:
    tuple val(meta), path("*.zip") , emit: fastqc_zip
    tuple val(meta), path("*.html"), emit: fastqc_html
    
    script:
    """
    ${params.fastqc_exec} -o ./ -f fastq ${read_1} ${read_2} -t $task.cpus
    """

    stub:
    """
    touch ${meta.id}.fastqc.zip
    touch ${meta.id}.fastqc.html
    """
}
