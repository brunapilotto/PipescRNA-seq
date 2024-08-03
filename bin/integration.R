#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(devtools)

devtools::install_github('immunogenomics/presto')


load_and_prepare <- function(paths) {
  objects <- lapply(paths, readRDS)
  return(objects)
}


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

before_integration_plot <- DimPlot(integrated_seurat, reduction = "umap") + ggtitle("Before Integration")
ggsave(filename = "before_integration_plot.png", plot = before_integration_plot, dpi = 300, width = 14, 
        height = 7, units = "in")

DefaultAssay(integrated_seurat) <- "integrated"
integrated_seurat <- ScaleData(integrated_seurat)
integrated_seurat <- RunPCA(integrated_seurat)
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "pca", dims = 1:30)

after_integration_plot <- DimPlot(integrated_seurat, reduction = "umap") + ggtitle("After Integration")
ggsave(filename = "after_integration_plot.png", plot = after_integration_plot, dpi = 300, width = 8, 
        height = 6, units = "in")

integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:30)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.5)
cluster_dim_plot <- DimPlot(integrated_seurat, reduction = "umap", label.size = 4, repel = TRUE, label = TRUE,
        group.by = "seurat_clusters") + ggtitle("Seurat Clustering")
ggsave(filename = "cluster_dim_plot.png", plot = cluster_dim_plot, dpi = 300, width = 8, 
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
