#!/usr/bin/env Rscript

BiocManager::install(version = "3.18")
bioc_packages <- c("IRanges", "SingleCellExperiment")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, force = TRUE, quietly = TRUE)
  }
}

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(Cairo)
library(SingleCellExperiment)
library(scater)
library(scuttle)

args <- commandArgs(trailingOnly=TRUE)
seurat_mtx <- args[1]
sample_name <- args[2]

seurat_object <- readRDS(seurat_mtx)
rm(seurat_mtx)

output_string <- sprintf("num_cells_before = %d\nnum_features_before = %d", ncol(seurat_object), nrow(seurat_object))
writeLines(output_string, "qc_metrics_before.txt")

seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
seurat_object[["percent.rb"]] <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL]")

# Visualize QC metrics as a violin plot
qc_before_plot <- VlnPlot(seurat_object, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"), ncol = 4,
                          cols = c("#e6867a"), pt.size = 0)
qc_before_plot[[1]] <- qc_before_plot[[1]] + ggtitle("nGenes") + theme(axis.title.x = element_blank())
qc_before_plot[[2]] <- qc_before_plot[[2]] + ggtitle("nUMIs") + theme(axis.title.x = element_blank())
qc_before_plot[[3]] <- qc_before_plot[[3]] + ggtitle("%Mit") + theme(axis.title.x = element_blank())
qc_before_plot[[4]] <- qc_before_plot[[4]] + ggtitle("%Rb") + theme(axis.title.x = element_blank())

qc_before_plot <- wrap_plots(qc_before_plot)
ggsave(filename = "qc_before.png", plot = qc_before_plot, dpi = 300, width = 12,
       height = 7, units = "in")

qc.count <- isOutlier(seurat_object@meta.data[["nCount_RNA"]], log=TRUE, type="lower")
qc.feat <- isOutlier(seurat_object@meta.data[["nFeature_RNA"]], log=TRUE, type="lower")
qc.mito <- isOutlier(seurat_object@meta.data[["percent.mt"]], type="higher")
qc.rib <- isOutlier(seurat_object@meta.data[["percent.rb"]], type="higher")

seurat_object[["discard"]] <- qc.feat | qc.count | qc.rib | qc.mito

plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", 
               group.by = "discard", cols = c("#c7fcd7", "#ed4a6a")) +
  ggtitle("UMIs vs %Mito") +
  xlab("nUMIs") +  # Altera o nome do eixo X
  ylab("%Mito") +  # Altera o nome do eixo Y
  labs(color = "Discard") +
  NoLegend()

plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
               group.by = "discard", cols = c("#c7fcd7", "#ed4a6a")) +
  ggtitle("UMIs vs Genes") +
  xlab("nUMIs") +  # Altera o nome do eixo X
  ylab("nGenes") +  # Altera o nome do eixo Y
  labs(color = "Discard")

feature_feature_relationships_plot <- plot1 + plot2

ggsave(filename = "feature-feature.png", plot = feature_feature_relationships_plot, dpi = 300, width = 11,
       height = 7, units = "in")

seurat_object <- subset(seurat_object, subset = discard == FALSE)

output_string <- sprintf("num_cells_after = %d\nnum_features_after = %d", ncol(seurat_object), nrow(seurat_object))
writeLines(output_string, "qc_metrics_after.txt")
saveRDS(seurat_object, file = paste0(sample_name, "_clean_object.rds"))

# Visualize QC metrics as a violin plot
qc_after_plot <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4,
                         cols = c("#e6867a"), pt.size = 0)
qc_after_plot[[1]] <- qc_after_plot[[1]] + ggtitle("nGenes") + theme(axis.title.x = element_blank())
qc_after_plot[[2]] <- qc_after_plot[[2]] + ggtitle("nUMIs") + theme(axis.title.x = element_blank())
qc_after_plot[[3]] <- qc_after_plot[[3]] + ggtitle("%Mit") + theme(axis.title.x = element_blank())
qc_after_plot[[4]] <- qc_after_plot[[4]] + ggtitle("%Rb") + theme(axis.title.x = element_blank())

qc_after_plot <- wrap_plots(qc_after_plot)
ggsave(filename = "qc_after.png", plot = qc_after_plot, dpi = 300, width = 12,
       height = 7, units = "in")
