process FASTP {
    conda 'bioconda::fastp=0.23.2'
    cpus 16
    
    input:
    path fastq_1
    path fastq_2
    val sample_name 
    path outdir

    output:
    path "${outdir}/fastp/${sample_name}_1_clean.fastq" 
    path "${outdir}/fastp/${sample_name}_2_clean.fastq"
    
    shell:
    """
    if [ ! -d "${outdir}/fastp" ]; then
        mkdir -p ${outdir}/fastp
        fastp fastp -i ${fastq_1} -I ${fastq_2} -o ${outdir}/fastp/${sample_name}_1_clean.fastq -O ${outdir}/fastp/${sample_name}_2_clean.fastq
        mv fastp.html ${outdir}/fastp/${sample_name}.html
        mv fastp.json ${outdir}/fastp/${sample_name}.json
    fi
    """
}