include { FASTQC        } from './modules/fastqc'
include { STAR          } from './modules/star'
include { MTX_TO_SEURAT } from './modules/mtx_conversion'
include { SEURAT        } from './modules/seurat_analysis'

log.info """\
    Pipeline lncsc-RNAseq
    ======================
    output: ${params.outdir}
""".stripIndent()

workflow {
    ch_metadado = Channel.fromPath(params.metadados).splitCsv(sep: ',', header: true)
                    .map { row -> tuple([id:row['ID']], row['fastq1'], row['fastq2'])}

    FASTQC( ch_metadado )
    STAR( ch_metadado )
    MTX_TO_SEURAT( STAR.out.filtered_counts )
    SEURAT( MTX_TO_SEURAT.out.seurat_object )
}
