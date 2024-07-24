// Include dotenv module to load environment variables
include { dotenv } from "plugin/nf-dotenv"

// Load environment variables
params.index_dir = dotenv("INDEX_DIR")
params.gtf_path  = dotenv("GTF_PATH")
params.star_exec = dotenv("STAR_EXEC")

// Define other parameters
params.white_list_path = "${workflow.projectDir}/assets/10x_V2_barcode_whitelist.txt"

// Include modules
include { FASTQC        } from './modules/fastqc'
include { STARSOLO      } from './modules/starsolo'
include { MTX_TO_SEURAT } from './modules/mtx_conversion'
include { QC_SEURAT     } from './modules/qc'
include { SINGLE_SEURAT } from './modules/single_seurat_analysis'

log.info """\
    Pipeline lncsc-RNAseq
    ======================
    output: ${params.outdir}
""".stripIndent()

workflow {
    ch_metadata = Channel.fromPath(params.metadata)
                    .splitCsv(sep: ',', header: true)
                    .map { row -> 
                            tuple([id: row['ID']], [fastq1: row['fastq1'], fastq2: row['fastq2']])
                        }
    ch_metadata.count().view()
    FASTQC ( ch_metadata )
    STARSOLO ( ch_metadata )
    MTX_TO_SEURAT ( STARSOLO.out.filtered_counts )
    QC_SEURAT ( MTX_TO_SEURAT.out.seurat_object )
    SINGLE_SEURAT ( QC_SEURAT.out.clean_object )
}
