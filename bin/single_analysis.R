#!/usr/bin/env Rscript

BiocManager::install(version = "3.18")
BiocManager::install()
BiocManager::install("IRanges", force = TRUE)
BiocManager::install("SingleCellExperiment", force = TRUE)
BiocManager::install("PCAtools", force = TRUE)
devtools::install_github('immunogenomics/presto')
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(devtools)
library(Cairo)
library(SingleCellExperiment)
library(DoubletFinder)
library(PCAtools)

args <- commandArgs(trailingOnly=TRUE)
clean_seurat_object <- args[1]
sample_name <- args[2]

seurat_object <- readRDS(clean_seurat_object)

seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object), 10)

# plot variable features with labels
variable_features_plot <- LabelPoints(plot = VariableFeaturePlot(seurat_object), points = top10, 
                                          repel = TRUE, xnudge = 0, ynudge = 0)
ggsave(filename = "variable_features.png", plot = variable_features_plot, dpi = 300, width = 10, 
        height = 7, units = "in")


all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)
rm(all.genes)
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

pca_stdev <- Stdev(seurat_object, reduction = "pca")
pca_variance <- pca_stdev^2
pca_variance_explained_percent <- (pca_variance / sum(pca_variance)) * 100
chosen.elbow <- findElbowPoint(pca_variance_explained_percent)
if (chosen.elbow < 10) {
  chosen.elbow <- 10
}
output_string <- sprintf("Elbow point in the percentage of variance explained by successive PCs: %d", chosen.elbow)
writeLines(output_string, "dimensionality.txt")
pdf(file = "elbow_plot.pdf", width = 8, height = 6)
par(mar = c(5, 4, 4, 2) + 0.1)
plot(pca_variance_explained_percent, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
title(main = "Variance Explained by Principal Components")
dev.off()

seurat_object <- FindNeighbors(seurat_object, dims = 1:chosen.elbow)
seurat_object <- FindClusters(seurat_object, resolution = 0.55)
seurat_object <- RunUMAP(seurat_object, dims = 1:chosen.elbow)


######## Remove Doublets ########

## pK Identification
x.res.list <- paramSweep(seurat_object, PCs = 1:chosen.elbow, sct = FALSE)
x.stats <- summarizeSweep(x.res.list, GT = FALSE)
bcmvn <- find.pK(x.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

# Homotypic Doublet Proportion Estimate
DoubletRate <- nrow(seurat_object@meta.data)*8*1e-6
homotypic.prop <- modelHomotypic(seurat_object$seurat_clusters) 
nExp_poi <- round(DoubletRate*nrow(seurat_object@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies 
seurat_object <- doubletFinder(seurat_object, PCs = 1:chosen.elbow, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

DF.name = colnames(seurat_object@meta.data)[grepl("DF.classification", colnames(seurat_object@meta.data))]
doublets_plot <- DimPlot(seurat_object, reduction = "umap", group.by = DF.name, 
                            label.size = 4, repel = TRUE) + ggtitle("DF Classification")
ggsave(filename = "doublets_plot.png", plot = doublets_plot, dpi = 300, height = 5, units = "in")
seurat_object <- seurat_object[, seurat_object@meta.data[, DF.name] == "Singlet"]

output_string <- sprintf("num_cells_final = %d\nnum_features_final = %d", ncol(seurat_object), nrow(seurat_object))
writeLines(output_string, "qc_metrics_final.txt")

umap_plot <- DimPlot(seurat_object, reduction = "umap", label.size = 4, repel = TRUE, label = TRUE,
                     group.by = "seurat_clusters") + ggtitle("Seurat Clustering")
ggsave(filename = "umap_plot.png", plot = umap_plot, dpi = 300, height = 5, width = 7, units = "in")


# find markers for every cluster compared to all remaining cells, report only the positive ones
cluster_markers <- FindAllMarkers(seurat_object, only.pos = TRUE)
write.table(cluster_markers, "total_cluster_marker_genes.txt", sep="\t", quote = F, row.names = F)
cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
write.table(top5, "top5_cluster_marker_genes.txt", sep="\t", quote = F, row.names = F)
heatmap <- DoHeatmap(seurat_object, features = top5$gene)
ggsave(filename = "heatmap.png", plot = heatmap, dpi = 300, width = 15, height = 8, units = "in")

saveRDS(seurat_object, file = paste0(sample_name, "_single_final_cluster.rds"))
