FASTQC_PATH = "/data1/biotools/FastQC/fastqc"

process FASTQC {
    label 'process_single'
    
    input:
    path read_1
    path read_2
    path outdir

    output:
    path "${outdir}/fastqc/*.zip", emit: fastqc_zip
    path "${outdir}/fastqc/*.html", emit: fastqc_html
    
    script:
    """
    mkdir -p ${outdir}/fastqc
    ${FASTQC_PATH} -o ${outdir}/fastqc -f fastq ${read_1} ${read_2} -t $task.cpus
    """
}