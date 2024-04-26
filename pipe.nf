params.outdir = null
params.fastq_path = null

include { FASTP } from './tools/fastp/fastp'

log.info """\
    Pipeline lncsc-RNAseq
    ======================
    fastq_path: ${params.fastq_path}
    output: ${params.outdir}
""".stripIndent()

workflow {
    def lastFolderName = file(params.fastq_path).getName()
    fastq_files = channel.fromPath( "${params.fastq_path}/${lastFolderName}_{1,2}.fastq", checkIfExists: true )
    .toSortedList()
    
    // Run Fastp
    FASTP()

}