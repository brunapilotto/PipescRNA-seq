#!/usr/bin/env Rscript
library(Seurat)

args <- commandArgs(trailingOnly=TRUE)

mtx_file      <- args[1]
barcode_file  <- args[2]
feature_file  <- args[3]
out.file      <- args[4]
sample_name   <- args[5]

expression.matrix <- ReadMtx(mtx = mtx_file, features = feature_file, cells = barcode_file)
seurat.object <- CreateSeuratObject(counts = expression.matrix, project = sample_name)

saveRDS(seurat.object, file = out.file)
