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

color_list <- c(
    "#264653",  # Dark Teal
    "#2a9d8f",  # Emerald Green
    "#e9c46a",  # Mustard Yellow
    "#f4a261",  # Pastel Orange
    "#e76f51",  # Burnt Sienna
    "#6b705c",  # CadetDark Olive Green
    "#457b9d",  # Steel Blue
    "#a8dadc",  # Light Blue
    "#f7a399",  # Peach Pink
    "#ffe8d6",  # Light Beige
    "#a5a58d",  # Sage Gray
    "#474747",  # Dark Gray
    "#9a8c98",  # Grayish Lilac
    "#c3aed6",  # Light Lavender
    "#5e6472",  # Slate Blue
    "#e63946",  # Crimson Red
    "#adc178",  # Light Olive Green
    "#93e1d8",  # Pastel Aqua
    "#606c38",  # Dark Olive Green
    "#1d3557"   # Deep Navy Blue
)

args <- commandArgs(trailingOnly=TRUE)
seurat_object <- args[1]
sample_name <- args[2]
signatures_json <- args[3]

seurat_object <- readRDS(seurat_object)
signatures <- fromJSON(signatures_json)

if (sample_name == "Integrated") {
  counts <- GetAssayData(seurat_object, layer = "data")
} else {
  counts <- GetAssayData(seurat_object, layer = "counts")
}
ranking <- AUCell_buildRankings(counts)
genes_in_seurat <- rownames(seurat_object)

########## RB cells ##########

genes_list <- list()
cell_types <- c(
  "Retinoblastoma",
  "B cells",
  "T cells",
  "NK cells",
  "CAF",
  "Progenitors"
)

for (category in cell_types) {
  genes_list[[category]] <- signatures[[category]]$geneSymbols
}

filtered_genes_list <- list()
for (category in names(genes_list)) {
  genes_present <- intersect(genes_list[[category]], genes_in_seurat)
  percentage_present <- length(genes_present) / length(genes_list[[category]]) * 100
  
  if (percentage_present > 20) {
    filtered_genes_list[[category]] <- genes_list[[category]]
  } else {
    print(paste(category, "Cell type not present in the dataset"))
  }
}
rm(genes_list)

if (!"CellType" %in% colnames(seurat_object@meta.data)) {
  seurat_object$CellType <- NA
}

for (category in names(filtered_genes_list)) {
  cell_AUC <- AUCell_calcAUC(filtered_genes_list[[category]], ranking)
  cell_assigment <- AUCell_exploreThresholds(cell_AUC, plotHist = TRUE, assign=TRUE)
  seurat_object$CellType <- ifelse(colnames(seurat_object) %in% cell_assigment$geneSet$assignment, 
                              gsub("_", " ", category), 
                              seurat_object$CellType)
}

cells_type_plot_umap <- DimPlot(seurat_object, label = TRUE, repel = TRUE, reduction = "umap",
                            label.size = 3, group.by = "CellType", cols=color_list) +
                            ggtitle("Cell types")
ggsave(filename = "cells_type_plot_umap.png", plot = cells_type_plot_umap, dpi = 300, height=7, width=10, units = "in")
cells_type_plot_tsne <- DimPlot(seurat_object, label = TRUE, repel = TRUE, reduction = "tsne",
                            label.size = 3, group.by = "CellType", cols=color_list) +
                            ggtitle("Cell types")
ggsave(filename = "cells_type_plot_tsne.png", plot = cells_type_plot_tsne, dpi = 300, height=7, width=10, units = "in")

cell_type_counts_df <- as.data.frame(table(seurat_object$CellType))
colnames(cell_type_counts_df) <- c("CellType", "Count")
write.table(cell_type_counts_df, "cell_type_counts.txt", sep="\t", quote = F, row.names = F)

########## RB cell types ##########

genes_list <- list()
cell_in_retinoblastoma <- c(
  "Glial",
  "CAF",
  "Rod-like",
  "Cone-like",
  "Cone_precursor-like",
  "MKI67+Cone_precursor",
  "Neural_cells"
)

for (category in cell_in_retinoblastoma) {
  genes_list[[category]] <- signatures[[category]]$geneSymbols
}

filtered_genes_list <- list()
for (category in names(genes_list)) {
  genes_present <- intersect(genes_list[[category]], genes_in_seurat)
  percentage_present <- length(genes_present) / length(genes_list[[category]]) * 100
  
  if (percentage_present > 20) {
    filtered_genes_list[[category]] <- genes_list[[category]]
  } else {
    print(paste(category, "Cell type not present in the dataset"))
  }
}
rm(genes_list)

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

rb_cells_plot_umap <- DimPlot(seurat_object, label = TRUE, repel = TRUE, reduction = "umap",
                            label.size = 3, group.by = "RBcells", cols=color_list) +
                            ggtitle("Cell types in human RB")
ggsave(filename = "rb_cells_plot_umap.png", plot = rb_cells_plot_umap, dpi = 300, height=7, width=10, units = "in")
rb_cells_plot_tsne <- DimPlot(seurat_object, label = TRUE, repel = TRUE, reduction = "tsne",
                            label.size = 3, group.by = "RBcells", cols=color_list) +
                            ggtitle("Cell types in human RB")
ggsave(filename = "rb_cells_plot_tsne.png", plot = rb_cells_plot_tsne, dpi = 300, height=7, width=10, units = "in")

cell_type_in_rb_counts <- table(seurat_object$RBcells)
cell_type_counts_df <- as.data.frame(cell_type_in_rb_counts)
colnames(cell_type_counts_df) <- c("CellType", "Count")
write.table(cell_type_counts_df, "cell_type_in_rb_counts.txt", sep="\t", quote = F, row.names = F)

######### EMTome #########

genes_present <- intersect(signatures$EMT$geneSymbols, genes_in_seurat)
percentage_present <- length(genes_present) / length(signatures$EMT$geneSymbols) * 100
if (percentage_present > 20) {
  cell_AUC <- AUCell_calcAUC(signatures$EMT$geneSymbols, ranking)
  cell_assigment <- AUCell_exploreThresholds(cell_AUC, plotHist = TRUE, assign=TRUE)
  seurat_object$EMT <- ifelse(colnames(seurat_object) %in% cell_assigment$geneSet$assignment, 
                              "EMT", 
                              NA)
  emt_plot_umap <- DimPlot(seurat_object, repel = TRUE, reduction = "umap",
                            label.size = 3, group.by = "EMT") +
                            ggtitle("EMT")
  ggsave(filename = "emt_plot_umap.png", plot = emt_plot_umap, dpi = 300, height=5, units = "in")
  emt_plot_tsne <- DimPlot(seurat_object, repel = TRUE, reduction = "tsne",
                            label.size = 3, group.by = "EMT") +
                            ggtitle("EMT")
  ggsave(filename = "emt_plot_tsne.png", plot = emt_plot_tsne, dpi = 300, height=5, units = "in")

  EMT_counts <- table(seurat_object$EMT)
  EMT_counts_df <- as.data.frame(EMT_counts)
  colnames(EMT_counts_df) <- c("CellType", "Count")
  write.table(EMT_counts_df, "EMT_counts.txt", sep="\t", quote = F, row.names = F)
} else {
  print("EMT hallmark not present in the dataset")
}

saveRDS(seurat_object, file = paste0(sample_name, "_signatures.rds"))
