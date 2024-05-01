FASTQC_PATH = "/data1/biotools/FastQC/fastqc"

process FASTQC {
    
    input:
    path clean_1
    path clean_2
    path outdir

    output:
    stdout
    
    shell:
    """
    if [ -z "\$(find "${outdir}/fastqc" -mindepth 1 -print -quit)" ]; then
        mkdir -p ${outdir}/fastqc
        ${FASTQC_PATH} -o ${outdir}/fastqc -f fastq ${clean_1} ${clean_2} -t 20
        rm ${outdir}/fastqc/*.zip
    fi
    """
}