#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(devtools)
library(PCAtools)

devtools::install_github('immunogenomics/presto')


load_and_prepare <- function(paths) {
  objects <- lapply(paths, readRDS)
  return(objects)
}

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

args <- commandArgs(trailingOnly = TRUE)
seurat_objects <- load_and_prepare(args)


######## Perform integration ########

anchors <- FindIntegrationAnchors(object.list = seurat_objects, dims = 1:30)
integrated_seurat <- IntegrateData(anchorset = anchors, dims = 1:30)

rm(seurat_objects)
rm(anchors)

DefaultAssay(integrated_seurat) <- "RNA"

integrated_seurat <- NormalizeData(integrated_seurat)
integrated_seurat <- FindVariableFeatures(integrated_seurat, selection.method = "vst", nfeatures = 2000)
integrated_seurat <- ScaleData(integrated_seurat)
integrated_seurat <- RunPCA(integrated_seurat)
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "pca", dims = 1:30)

before_integration_plot <- DimPlot(integrated_seurat, reduction = "umap", cols=color_list) + ggtitle("Before Integration")
ggsave(filename = "before_integration_plot.png", plot = before_integration_plot, dpi = 300, width = 14, 
        height = 7, units = "in")

DefaultAssay(integrated_seurat) <- "integrated"
integrated_seurat <- ScaleData(integrated_seurat)
integrated_seurat <- RunPCA(integrated_seurat)

pca_stdev <- Stdev(integrated_seurat, reduction = "pca")
pca_variance <- pca_stdev^2
pca_variance_explained_percent <- (pca_variance / sum(pca_variance)) * 100
chosen.elbow <- findElbowPoint(pca_variance_explained_percent)
output_string_1 <- sprintf("Elbow point in the percentage of variance explained by successive PCs: %d", chosen.elbow)
writeLines(output_string_1, "elbow_point.txt")
output_string_2 <- sprintf("Percentage of variance explained by the first %d PCs: %f", chosen.elbow, sum(pca_variance_explained_percent[1:chosen.elbow]))
writeLines(output_string_2, "percentage_of_variance.txt")
pdf(file = "elbow_plot.pdf", width = 8, height = 6)
par(mar = c(5, 4, 4, 2) + 0.1)
plot(pca_variance_explained_percent, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
title(main = "Variance Explained by Principal Components")
dev.off()

integrated_seurat <- RunTSNE(integrated_seurat, reduction = "pca", dims = 1:chosen.elbow)
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "pca", dims = 1:chosen.elbow)

after_integration_plot <- DimPlot(integrated_seurat, reduction = "umap", cols=color_list) + ggtitle("After Integration")
ggsave(filename = "after_integration_plot.png", plot = after_integration_plot, dpi = 300, width = 8, 
        height = 6, units = "in")

integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:chosen.elbow)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.5)
cluster_dim_plot_umap <- DimPlot(integrated_seurat, reduction = "umap", label.size = 4, repel = TRUE, label = TRUE,
        group.by = "seurat_clusters", cols=color_list) + ggtitle("Seurat Clustering")
ggsave(filename = "cluster_dim_plot_umap.png", plot = cluster_dim_plot_umap, dpi = 300, width = 8, 
        height = 6, units = "in")
cluster_dim_plot_tsne <- DimPlot(integrated_seurat, reduction = "tsne", label.size = 4, repel = TRUE, label = TRUE,
        group.by = "seurat_clusters", cols=color_list) + ggtitle("Seurat Clustering")
ggsave(filename = "cluster_dim_plot_tsne.png", plot = cluster_dim_plot_tsne, dpi = 300, width = 8, 
        height = 6, units = "in")

saveRDS(integrated_seurat, file = "integrated_obj.rds")

count_table <- table(integrated_seurat@meta.data$seurat_clusters, integrated_seurat@meta.data$orig.ident)
prop_table <- prop.table(count_table, margin = 1)
percent_table <- prop_table * 100
percent_df <- as.data.frame(as.table(percent_table))
colnames(percent_df) <- c("Cluster", "Sample", "Percentage")
cells_in_cluster <- ggplot(percent_df, aes(fill=Sample, y=Cluster, x=Percentage)) + 
  geom_bar(position="fill", stat="identity", colour = "black") +
  scale_fill_brewer(palette = "Pastel1") +
  scale_x_continuous(labels = scales::percent_format(scale = 100)) +
  xlab("Proportion") +
  ylab("Clusters") +
  ggtitle("Fraction of cells in each cell cluster") +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"))
ggsave(filename = "cells_in_cluster.png", plot = cells_in_cluster, dpi = 300, width = 8, 
       height = 6, units = "in")

output_string_3 <- sprintf("num_cells_final = %d\nnum_features_final = %d", ncol(integrated_seurat), nrow(integrated_seurat))
writeLines(output_string_3, "qc_metrics_final.txt")

# find markers for every cluster compared to all remaining cells, report only the positive ones
cluster_markers <- FindAllMarkers(integrated_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(cluster_markers, "total_cluster_marker_genes.txt", sep="\t", quote = F, row.names = F)
cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
write.table(top5, "top5_cluster_marker_genes.txt", sep="\t", quote = F, row.names = F)
heatmap <- DoHeatmap(integrated_seurat, features = top5$gene, size = 4) + NoLegend()
ggsave(filename = "heatmap.png", plot = heatmap, dpi = 300, width = 20, height = 10, units = "in")
