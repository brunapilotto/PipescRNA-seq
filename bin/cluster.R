#!/usr/bin/env Rscript

BiocManager::install(version = "3.18")
BiocManager::install("SingleCellExperiment", force = TRUE)
devtools::install_github('immunogenomics/presto')

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(devtools)
library(celldex)
library(SingleR)
library(Cairo)
library(SingleCellExperiment)

args <- commandArgs(trailingOnly=TRUE)
cone <- c("GNGT2", "RXRG", "CRX", "PDE6H")
retinoblastoma <- c("MKI67", "CDK1", "TOP2A", "KIF14", "CDC25C")

clean_seurat_object <- args[1]

seurat_object <- readRDS(clean_seurat_object)

seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object), 10)

# plot variable features with labels
variable_features_plot <- LabelPoints(plot = VariableFeaturePlot(seurat_object), points = top10, 
                                          repel = TRUE, xnudge = 0, ynudge = 0)
ggsave(filename = "variable_features.png", plot = variable_features_plot, dpi = 300, width = 10, 
        height = 7, units = "in")


all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
seurat_object <- RunTSNE(seurat_object, dims = 1:10, verbose = FALSE)

# Plot
tsne_plot <- DimPlot(seurat_object, reduction = "tsne", label.size = 4, repel = TRUE, label = TRUE,
                     group.by = "seurat_clusters") + ggtitle("Seurat Clustering")
ggsave(filename = "tsne_plot.png", plot = tsne_plot, dpi = 300, height = 5, units = "in")

# Cone expression plot
cone_expression_feature_plot <- FeaturePlot(seurat_object, features = cone)
ggsave(filename = "cone_expression_feature.png", plot = cone_expression_feature_plot, dpi = 300, 
       width = 10, units = "in")

cone_expression_violin_plot <- VlnPlot(seurat_object, features = cone)
ggsave(filename = "cone_expression.png", plot = cone_expression_violin_plot, dpi = 300, width = 11,
       height = 7, units = "in")

# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
cone_expression_ridge_plot <- RidgePlot(seurat_object, features = cone, ncol = 2)
ggsave(filename = "cone_expression_ridge.png", plot = cone_expression_ridge_plot, dpi = 300, width = 7,
       height = 7, units = "in")

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
cone_expression_dot_plot <- DotPlot(seurat_object, features = cone) + RotatedAxis()
ggsave(filename = "cone_expression_dot.png", plot = cone_expression_dot_plot, dpi = 300, width = 8,
       height = 8, units = "in")

# Retinoblastoma expression plot
retinoblastoma_expression_violin_plot <- VlnPlot(seurat_object, features = retinoblastoma)
ggsave(filename = "retinoblastoma_expression.png", plot = retinoblastoma_expression_violin_plot, dpi = 300,
       width = 8, units = "in")

# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
retinoblastoma_expression_ridge_plot <- RidgePlot(seurat_object, features = retinoblastoma, ncol = 2)
ggsave(filename = "retinoblastoma_expression_ridge.png", plot = retinoblastoma_expression_ridge_plot, dpi = 300, width = 8,
       height = 10, units = "in")

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
retinoblastoma_expression_dot_plot <- DotPlot(seurat_object, features = retinoblastoma) + RotatedAxis()
ggsave(filename = "retinoblastoma_expression_dot.png", plot = retinoblastoma_expression_dot_plot, dpi = 300, width = 8,
       height = 8, units = "in")

retinoblastoma_expression_feature_plot <- FeaturePlot(seurat_object, features = retinoblastoma)
ggsave(filename = "retinoblastoma_expression_feature.png", plot = retinoblastoma_expression_feature_plot, dpi = 300)

# find markers for every cluster compared to all remaining cells, report only the positive ones
cluster_markers <- FindAllMarkers(seurat_object, only.pos = TRUE)
# Mostrar os principais marcadores para cada cluster
cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
heatmap <- DoHeatmap(seurat_object, features = top5$gene)
ggsave(filename = "heatmap.png", plot = heatmap, dpi = 300, width = 15, height = 8, units = "in")

saveRDS(seurat_object, file = "single_final_cluster.rds")

monaco.ref <- celldex::MonacoImmuneData()
sce <- as.SingleCellExperiment(DietSeurat(seurat_object))
predictions <- SingleR(test=sce, assay.type.test=1, 
                       ref=monaco.ref, labels=monaco.ref$label.main)
table(predictions$pruned.labels)

seurat_object@meta.data$predictions <- predictions$pruned.labels
seurat_object <- SetIdent(seurat_object, value = "monaco.main")
predictions_plot <- DimPlot(seurat_object, label = TRUE, repel = TRUE, 
                            label.size = 3, group.by = "predictions") +
                            ggtitle("SingleR Predictions")
ggsave(filename = "predictions_plot.png", plot = predictions_plot, dpi = 300, height=5, units = "in")
