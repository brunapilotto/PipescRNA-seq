#!/usr/bin/env Rscript
BiocManager::install(version = "3.18")
bioc_packages <- c("mixtools")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, force = TRUE, quietly = TRUE)
  }
}

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(jsonlite)
library(AUCell)

args <- commandArgs(trailingOnly=TRUE)
seurat_object <- args[1]
sample_name <- args[2]
retinoblastoma_json <- args[3]

seurat_object <- readRDS(seurat_object)
retinoblastoma_cell_types <- fromJSON(retinoblastoma_json)

########## RB cell types ##########

genes_list <- list()

cell_types <- c(
  "Glial",
  "Cancer-associated_Fibroblasts",
  "Rod-like",
  "Cone-like",
  "Cone_precursor-like",
  "MKI67+Cone_precursor",
  "Neural_cells"
)

for (category in cell_types) {
  genes_list[[category]] <- retinoblastoma_cell_types[[category]]$geneSymbols
}

genes_in_seurat <- rownames(seurat_object)
filtered_genes_list <- list()
for (category in names(genes_list)) {
  genes_present <- intersect(genes_list[[category]], genes_in_seurat)
  percentage_present <- length(genes_present) / length(genes_list[[category]]) * 100
  
  if (percentage_present > 20) {
    filtered_genes_list[[category]] <- genes_list[[category]]
  }
}
rm(genes_list)

if (sample_name == "Integrated") {
  counts <- GetAssayData(seurat_object, layer = "data")
} else {
  counts <- GetAssayData(seurat_object, layer = "counts")
}

ranking <- AUCell_buildRankings(counts)

if (!"RBcells" %in% colnames(seurat_object@meta.data)) {
  seurat_object$RBcells <- NA
}

for (category in names(filtered_genes_list)) {
  cell_AUC <- AUCell_calcAUC(filtered_genes_list[[category]], ranking)
  cell_assigment <- AUCell_exploreThresholds(cell_AUC, plotHist = TRUE, assign=TRUE)
  seurat_object$RBcells <- ifelse(colnames(seurat_object) %in% cell_assigment$geneSet$assignment, 
                              gsub("_", " ", category), 
                              seurat_object$RBcells)
}

rb_cells_plot <- DimPlot(seurat_object, label = TRUE, repel = TRUE, 
                            label.size = 3, group.by = "RBcells") +
                            ggtitle("Cell types in human RB")
ggsave(filename = "rb_cells_plot.png", plot = rb_cells_plot, dpi = 300, height=7, width=10, units = "in")

########## RB cells ##########

genes_present <- intersect(retinoblastoma_cell_types$Retinoblastoma$geneSymbols, genes_in_seurat)
percentage_present <- length(genes_present) / length(retinoblastoma_cell_types$Retinoblastoma$geneSymbols) * 100
if (percentage_present > 20) {
  cell_AUC <- AUCell_calcAUC(retinoblastoma_cell_types$Retinoblastoma$geneSymbols, ranking)
  cell_assigment <- AUCell_exploreThresholds(cell_AUC, plotHist = TRUE, assign=TRUE)
  seurat_object$RB <- ifelse(colnames(seurat_object) %in% cell_assigment$geneSet$assignment, 
                              "RB", 
                              NA)
  RB_plot <- DimPlot(seurat_object, repel = TRUE, 
                            label.size = 3, group.by = "RB") +
                            ggtitle("Retinoblastoma cells")
  ggsave(filename = "RB_plot.png", plot = RB_plot, dpi = 300, height=5, units = "in")
} else {
  print("Retinoblastoma cells not present in the dataset")
}
