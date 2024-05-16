process MULTIQC {
    conda "/home/mbruna/PipescRNA-seq/envs/multiqc.yml"
    
    input:
    path fastqc_zip
    path star_log
    path outdir

    output:
    stdout
    
    shell:
    """
    if [ -z "\$(find "${outdir}/multiqc" -mindepth 1 -print -quit)" ]; then
        mkdir -p ${outdir}/multiQC
        multiqc . -o ${outdir}/multiqc
    fi
    """
}