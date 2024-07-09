# PipescRNA-seq

This pipeline is a bioinformatics analysis for processing 10x Genomics single-cell RNA-seq data.

## STARsolo

[STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) is a tool for analyzing droplet single cell RNA sequencing data and it is used in this pipeline.

To make STARsolo raw gene counts (almost) identical to CellRanger's, first, you have to use [CellRanger's human reference](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads) genome. Annotations affect the counts, and, to match CellRanger counts, CellRanger annotations have to be used.

### Genome index generation

The FASTA and GTF files, for the latest release have to be used in STAR genome index generation step before mapping. To make the agreement between STARsolo and CellRanger even more perfect, you can add `--genomeSAsparseD 3` to the genome generation options, which is used by CellRanger to generate STAR genomes:

```bash
STAR --runMode genomeGenerate --runThreadN ... --genomeDir ./ --genomeFastaFiles refdata-gex-GRCh38-2020-A/fasta/genome.fa  --sjdbGTFfile refdata-gex-GRCh38-2020-A/genes/genes.gtf --genomeSAsparseD 3
```

### Cell filtering (calling)

In addition to unfiltered output of gene/cell counts, STARsolo performs cell filtering (a.k.a. cell calling), which aims to select a subset of cells that are likely to be "real" cells as opposed to empty droplets (containing ambient RNA). This option is used in this pipeline:

```bash
STAR --genomeDir /path/to/index/ \
--readFilesIn /path/to/fastq2 /path/to/fastq1 \
--runThreadN <CPUS> \
--outFileNamePrefix <SAMPLE_NAME> \
--soloCBwhitelist /path/to/white_list.txt \
--sjdbGTFfile /path/to/genes.gtf \
--readFilesCommand zcat \
--soloType CB_UMI_Simple \
--outSAMtype None \
--soloBarcodeReadLength 0 \
--runDirPerm All_RWX \
--soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
--soloUMIfiltering MultiGeneUMI_CR \
--soloUMIdedup 1MM_CR
```

The 10X Chromium whitelist files can be found inside the [CellRanger distribution](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist). A copy of the whitelist for 5' v2 chemistry is available in the pipeline `assets/` folder. Make sure that the whitelist is compatible with the specific version of the 10X chemistry.

## Pipeline Usage

> [!NOTE]
> Before running the pipeline, make sure you have all dependencies in your machine (check [Dependencies](#dependencies)).

### Samplesheet

First, prepare a sample sheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
ID,fastq_1,fastq_2
SRR14800535,SRR14800535_1.fastq.gz,SRR14800535_2.fastq.gz
SRR14800536,SRR14800536_1.fastq.gz,SRR14800536_2.fastq.gz
```

Each row represents a pair of fastq files (paired end).

### Command

Now, you can run the pipeline using:

```bash
nextflow run PipescRNA-seq/main.nf --outdir <OUTDIR> --metadata samplesheet.csv
```

where:

- **--metadata**: path to comma-separated file containing information about the samples in the experiment.

- **--outdir**: the output directory where the results will be saved.

### After finishing

Once the pipeline is finished, you can check the results in the [Nextflow work directory](#nextflow-work-directory) and/or your output directory.

Once you finish the analysis, you can remove the last run files from the Nextflow work directory.

Ex: if your Nextflow work directory is `/tmp/nextflow/`, you can use the following command:

```bash
rm -rf /tmp/nextflow/*
```

> [!WARNING]
> Be careful! If the pipeline was stopped by an error, do not delete the information on this folder, because the `-resume` option needs this information to continue executions that were stopped by an error. Only delete it if the pipeline was completed with no error, and you are not running other analysis.

## Dependencies

### Nextflow

It is required to have Nextflow installed on your machine. Please visit the [Nextflow Documentation](https://www.nextflow.io/docs/latest/install.html) for more information.

### Conda

It is required to have Conda installed on your machine. Please visit the [Conda Documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) for more information.

Once installed, please change the `params.config` file with the path of the installation.

### Environment variables (.env)

In the project folder, it is necessary to create a `.env` file with the following information:

- `INDEX_DIR`: path to the precomputed STAR index
- `GTF_PATH`: reference GTF annotation file

Example .env file:

```bash
INDEX_DIR="/some/path/reference_cell_ranger/star_index"
GTF_PATH="/some/path/reference_cell_ranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
```

### Nextflow work directory

Please add this line at your `.bashrc` file to set nextflow work directory:

```bash
export NXF_WORK=/tmp/nextflow
```

### Tower config

You can monitor your pipeline execution in [Seqera Cloud](https://cloud.seqera.io/) page. To do this, you must have an account and access token, which need to be specified at the `tower.config` file in the `conf/` folder.

```groovy
tower {
    enabled     = true
    accessToken = "token"
}
```

### Pipeline parameters

Another pipeline parameters must be specified at the `params.config` file in the `conf/` folder. The following information are required:

- `ANACONDA3`: path to Conda folder

File example:

```groovy
params {
    ANACONDA3 = "/opt/anaconda3"
}
```
