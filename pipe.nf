params.outdir = null
params.fastq_path = null

include { FASTP } from './tools/fastp/fastp'
include { FASTQC } from './tools/fastqc/fastqc'

log.info """\
    Pipeline lncsc-RNAseq
    ======================
    fastq_path: ${params.fastq_path}
    output: ${params.outdir}
""".stripIndent()

def getFastqFiles(directory) {
    def dir = file(directory)
    def files = dir.listFiles().findAll { element -> element.getName().endsWith(".fastq") }.sort()
    return files
}

workflow {
    def lastFolderName = file(params.fastq_path).getName()
    def fastq_files = getFastqFiles(params.fastq_path)

    // Criação dos canais
    fastq1 = Channel.fromPath(fastq_files[0])
    fastq2 = Channel.fromPath(fastq_files[1])
    
    // Run Fastp
    ( clean_files_1, clean_files_2 ) = FASTP(fastq1, fastq2, lastFolderName, params.outdir)

    //Run FastQC
    FASTQC(clean_files_1, clean_files_2, params.outdir)
}