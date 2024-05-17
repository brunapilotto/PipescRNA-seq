output_dir = params.outdir ? file( params.outdir, checkIfExists: true ) : null
fastqs_path = params.fastq_path ? file( params.fastq_path, checkIfExists: true ) : null

include { FASTQC } from './modules/fastqc/fastqc'
include { MULTIQC } from './modules/multiqc/multiqc'
include { STAR } from './modules/star/star'
include { PIGZ } from './utils/pigz'
include { EMPTYDROPS_CELL_CALLING } from './modules/bioconductor/emptydrops'
include { MTX_TO_H5AD } from './modules/mtx_conversion/mtx_to_h5ad'

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
    ch_mtx_matrices = Channel.empty()
    ch_txp2gene = []

    // Criação dos canais
    fastq1 = Channel.fromPath(fastq_files[0])
    fastq2 = Channel.fromPath(fastq_files[1])

    // Run FastQC
    FASTQC(fastq1, fastq2, params.outdir)

    // Run STAR
    STAR(fastq1, fastq2, lastFolderName, params.outdir)
    ch_mtx_matrices = ch_mtx_matrices.mix(STAR.out.raw_counts, STAR.out.filtered_counts)

    // Run PIGZ
    PIGZ(STAR.out.raw_counts)

    // Run EmptyDrops
    EMPTYDROPS_CELL_CALLING(PIGZ.out.matrix_zip, PIGZ.out.barcodes_zip, PIGZ.out.features_zip, params.outdir)
    ch_mtx_matrices = ch_mtx_matrices.mix( EMPTYDROPS_CELL_CALLING.out.filtered_matrices )

    // Convert matrix to h5ad
    MTX_TO_H5AD(ch_mtx_matrices, ch_txp2gene, lastFolderName, params.outdir)

    // Run MultiQC
    // MULTIQC(FASTQC.out.fastqc_zip, STAR.out.log_final, params.outdir)
}