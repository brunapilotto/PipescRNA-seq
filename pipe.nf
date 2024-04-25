// params.fastq_1 = null
// params.fastq_2 = null

process fastp {
    conda 'bioconda::fastp=0.23.2'
    
    // input:
    // path fastq_1
    // path fastq_2

    // output:
    // path 'out.R1.fq.gz' into clean1
    // path 'out.R2.fq.gz' into clean2
    
    shell:
    '''
    fastp -?
    '''
}

workflow {
    fastp()
}