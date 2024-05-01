process FASTP {
    conda '/home/mbruna/PipescRNA-seq/envs/fastp.yml'
    
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
    if [ -z "\$(find "${outdir}/fastp" -mindepth 1 -print -quit)" ]; then
        mkdir -p ${outdir}/fastp
        fastp -i ${fastq_1} -I ${fastq_2} -o ${outdir}/fastp/${sample_name}_1_clean.fastq -O ${outdir}/fastp/${sample_name}_2_clean.fastq -q 20 --thread 16
        mv fastp.html ${outdir}/fastp/${sample_name}.html
        mv fastp.json ${outdir}/fastp/${sample_name}.json
    fi
    """
}