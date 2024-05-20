output_dir = params.outdir ? file( params.outdir, checkIfExists: true ) : null
fastqs_path = params.fastq_path ? file( params.fastq_path, checkIfExists: true ) : null

include { FASTQC } from './modules/fastqc/fastqc'
include { STAR } from './modules/star/star'
include { MTX_TO_SEURAT } from './modules/mtx_conversion/mtx_to_seurat'
include { SEURAT } from './modules/seurat_analysis/seurat'

log.info """\
    Pipeline lncsc-RNAseq
    ======================
    fastq_path: ${params.fastq_path}
    output: ${params.outdir}
""".stripIndent()

def getFastqFiles(directory) {
    def dir = file(directory)
    def files = dir.listFiles().findAll { element -> element.getName().endsWith(".fastq.gz") }.sort()
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

    // Convert matrix to seurat obj
    MTX_TO_SEURAT(STAR.out.filtered_counts, params.outdir, lastFolderName)

    // Make clusters
    SEURAT(MTX_TO_SEURAT.out.seuratObject, params.outdir)
}