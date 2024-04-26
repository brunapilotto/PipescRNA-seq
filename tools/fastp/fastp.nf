process FASTP {
    conda 'bioconda::fastp=0.23.2'
    
    // input:
    // path fastq_1
    // path fastq_2

    // output:
    // path "${lastFolderName}_1.fq.gz" into clean1
    // path "${lastFolderName}_2.fq.gz" into clean2
    
    shell:
    '''
    fastp -?
    '''
}