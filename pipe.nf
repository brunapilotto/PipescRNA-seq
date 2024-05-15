params.outdir = null
params.fastq_path = null

include { FASTQC } from './tools/fastqc/fastqc'
include { MULTIQC } from './tools/multiqc/multiqc'
include { STAR } from './tools/star/star'

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

    // Run FastQC
    FASTQC(fastq1, fastq2, params.outdir)

    // Run STAR
    STAR(fastq1, fastq2, lastFolderName, params.outdir)

    // Run MultiQC
    MULTIQC(FASTQC.out.fastqc_zip, STAR.out.log_final, params.outdir)
}