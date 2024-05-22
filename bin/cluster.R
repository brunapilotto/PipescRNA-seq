#!/usr/bin/env Rscript
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

args <- commandArgs(trailingOnly=TRUE)

seurat_mtx  <- args[1]
out_dir     <- args[2]

seurat_object <- readRDS(seurat_mtx)

#QC and selecting cells for further analysis
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

# Visualize QC metrics as a violin plot
qc_before_plot <- vln_plot <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Adicionar títulos personalizados aos gráficos
qc_before_plot[[1]] <- qc_before_plot[[1]] + ggtitle("nGenes")
qc_before_plot[[2]] <- qc_before_plot[[2]] + ggtitle("nUMIs")
qc_before_plot[[3]] <- qc_before_plot[[3]] + ggtitle("%Mit")

# Usar wrap_plots para juntar os gráficos com títulos personalizados
qc_before_plot <- wrap_plots(qc_before_plot)

qc_before_file <- paste0(out_dir, "/qc_before.png")
ggsave(filename = qc_before_file, plot = qc_before_plot, dpi = 300, width = 11,
       height = 7, units = "in")

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
feature_feature_relationships_plot <- plot1 + plot2

feature_feature_relationships_file <- paste0(out_dir, "/feature-feature.png")
ggsave(filename = feature_feature_relationships_file, plot = feature_feature_relationships_plot, dpi = 300, width = 11,
       height = 7, units = "in")

seurat_object <- subset(seurat_object, subset = percent.mt < 5)

# 1. Normalização dos dados
seurat_object <- NormalizeData(seurat_object)

# 2. Identificação das variáveis características
seurat_object <- FindVariableFeatures(seurat_object)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object), 10)

# plot variable features with and without labels
variable_features_plot <- LabelPoints(plot = VariableFeaturePlot(seurat_object), points = top10, repel = TRUE)
variable_features_file <- paste0(out_dir, "/variable_features.png")
ggsave(filename = variable_features_file, plot = variable_features_plot, dpi = 300, width = 10, 
        height = 7, units = "in")

# 3. Escalonamento dos dados
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)

# 4. Análise de componentes principais (PCA)
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

# 5. Construção do gráfico de vizinhos mais próximos
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)

# 6. Agrupamento das células
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

# 7. Redução de dimensionalidade com t-SNE
seurat_object <- RunTSNE(seurat_object, dims = 1:10)

# Plotando o t-SNE com os clusters
tsne_plot <- DimPlot(seurat_object, reduction = "tsne", label = TRUE)
tsne_plot_file <- paste0(out_dir, "/tsne_plot.png")
ggsave(filename = tsne_plot_file, plot = tsne_plot, dpi = 300)

# gráfico com os genes específicos
cone_expression_violin_plot <- VlnPlot(seurat_object, features = c("GNGT2", "RXRG", "CRX", "PDE6H"))
cone_expression_file <- paste0(out_dir, "/cone_expression.png")
ggsave(filename = cone_expression_file, plot = cone_expression_violin_plot, dpi = 300, width = 11,
       height = 7, units = "in")

# outra visualização
# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
cone_expression_ridge_plot <- RidgePlot(seurat_object, features = c("GNGT2", "RXRG", "CRX", "PDE6H"), ncol = 2)
cone_expression_ridge_file <- paste0(out_dir, "/cone_expression_ridge.png")
ggsave(filename = cone_expression_ridge_file, plot = cone_expression_ridge_plot, dpi = 300, width = 10,
       height = 7, units = "in")

# outra visualização
# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
cone_expression_dot_plot <- DotPlot(seurat_object, features = c("GNGT2", "RXRG", "CRX", "PDE6H")) + RotatedAxis()
cone_expression_dot_file <- paste0(out_dir, "/cone_expression_dot.png")
ggsave(filename = cone_expression_dot_file, plot = cone_expression_dot_plot, dpi = 300, width = 8,
       height = 8, units = "in")

retinoblastoma_expression_violin_plot <- VlnPlot(seurat_object, features = c("MKI67", "CDK1", "TOP2A", "KIF14", "CDC25C"))
retinoblastoma_expression_file <- paste0(out_dir, "/retinoblastoma_expression.png")
ggsave(filename = retinoblastoma_expression_file, plot = retinoblastoma_expression_violin_plot, dpi = 300)
# outra visualização
# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
retinoblastoma_expression_ridge_plot <- RidgePlot(seurat_object, features = c("MKI67", "CDK1", "TOP2A", "KIF14", "CDC25C"), ncol = 2)
retinoblastoma_expression_ridge_file <- paste0(out_dir, "/retinoblastoma_expression_ridge.png")
ggsave(filename = retinoblastoma_expression_ridge_file, plot = retinoblastoma_expression_ridge_plot, dpi = 300, width = 8,
       height = 10, units = "in")

# outra visualização
# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
retinoblastoma_expression_dot_plot <- DotPlot(seurat_object, features = c("MKI67", "CDK1", "TOP2A", "KIF14", "CDC25C")) + RotatedAxis()
retinoblastoma_expression_dot_file <- paste0(out_dir, "/retinoblastoma_expression_dot.png")
ggsave(filename = retinoblastoma_expression_dot_file, plot = retinoblastoma_expression_dot_plot, dpi = 300, width = 8,
       height = 8, units = "in")

cone_expression_feature_plot <- FeaturePlot(seurat_object, features = c("GNGT2", "RXRG", "CRX", "PDE6H"))
cone_expression_feature_file <- paste0(out_dir, "/cone_expression_feature.png")
ggsave(filename = cone_expression_feature_file, plot = cone_expression_feature_plot, dpi = 300)

retinoblastoma_expression_feature_plot <- FeaturePlot(seurat_object, features = c("MKI67", "CDK1", "TOP2A", "KIF14", "CDC25C"))
retinoblastoma_expression_feature_file <- paste0(out_dir, "/retinoblastoma_expression_feature.png")
ggsave(filename = retinoblastoma_expression_feature_file, plot = retinoblastoma_expression_feature_plot, dpi = 300)

# 8. Identificação de biomarcadores de cluster
# find markers for every cluster compared to all remaining cells, report only the positive ones
cluster_markers <- FindAllMarkers(seurat_object, only.pos = TRUE)
# Mostrar os principais marcadores para cada cluster
cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
heatmap <- DoHeatmap(seurat_object, features = top5$gene)
heatmap_file <- paste0(out_dir, "/heatmap.png")
ggsave(filename = heatmap_file, plot = heatmap, dpi = 300, width = 15, height = 8, units = "in")

cluster_seurat_file <- paste0(out_dir, "/cluster_final.rds")
saveRDS(seurat_object, file = cluster_seurat_file)
