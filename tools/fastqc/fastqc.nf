FASTQC_PATH = "/data1/biotools/FastQC/fastqc"

process FASTQC {
    
    input:
    path read_1
    path read_2
    path outdir

    output:
    stdout
    
    shell:
    """
    if [ -z "\$(find "${outdir}/fastqc" -mindepth 1 -print -quit)" ]; then
        mkdir -p ${outdir}/fastqc
        ${FASTQC_PATH} -o ${outdir}/fastqc -f fastq ${read_1} ${read_2} -t 20
        rm ${outdir}/fastqc/*.zip
    fi
    """
}