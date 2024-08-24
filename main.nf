// Include dotenv module to load environment variables
include { dotenv } from "plugin/nf-dotenv"

// Load environment variables
params.index_dir = dotenv("INDEX_DIR")
params.gtf_path  = dotenv("GTF_PATH")
params.star_exec = dotenv("STAR_EXEC")

// Include nf-schema module to validate the pipeline
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
validateParameters()

// Print summary of supplied parameters
log.info """\
Pipeline lncsc-RNAseq
======================
""".stripIndent()
log.info paramsSummaryLog(workflow)

// Include modules
include { FASTQC                                                   } from './modules/fastqc'
include { STARSOLO                                                 } from './modules/starsolo'
include { MTX_TO_SEURAT                                            } from './modules/mtx_conversion'
include { QC_SEURAT                                                } from './modules/qc'
include { SINGLE_SEURAT                                            } from './modules/single_seurat_analysis'
include { INTEGRATION                                              } from './modules/seurat_integration'
include { RETINOBLASTOMA_PLOTS as SINGLE_RETINOBLASTOMA_PLOTS      } from './modules/retinoblastoma_plots'
include { RETINOBLASTOMA_PLOTS as INTEGRATION_RETINOBLASTOMA_PLOTS } from './modules/retinoblastoma_plots'
include { ANNOTATION_PLOTS as SINGLE_ANNOTATION_PLOTS              } from './modules/annotation'
include { ANNOTATION_PLOTS as INTEGRATION_ANNOTATION_PLOTS         } from './modules/annotation'

workflow {
    ch_metadata = Channel.fromList(samplesheetToList(params.metadata, "assets/schema_input.json"))
                .map { meta, fastq1, fastq2, white_list -> 
                            tuple(meta, [fastq1: fastq1, fastq2: fastq2], white_list)
                        }

    FASTQC ( ch_metadata )
    STARSOLO ( ch_metadata )
    MTX_TO_SEURAT ( STARSOLO.out.filtered_counts )
    QC_SEURAT ( MTX_TO_SEURAT.out.seurat_object )
    SINGLE_SEURAT ( QC_SEURAT.out.clean_object )
    SINGLE_RETINOBLASTOMA_PLOTS ( SINGLE_SEURAT.out.single_seurat_object )
    SINGLE_ANNOTATION_PLOTS ( SINGLE_SEURAT.out.single_seurat_object )

    if ( params.sample_count >= 2 ) {
        single_objects = SINGLE_SEURAT.out.single_seurat_object
                    .map { meta, path -> path }
                    .collect()
        INTEGRATION ( single_objects )
        integrated_object = INTEGRATION.out.merged_object
                            .map { path -> tuple([id: 'Integrated'], path) }
        INTEGRATION_RETINOBLASTOMA_PLOTS ( integrated_object )
        INTEGRATION_ANNOTATION_PLOTS ( integrated_object )
    }
}
