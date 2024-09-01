#!/usr/bin/env Rscript

BiocManager::install(version = "3.18")
bioc_packages <- c("IRanges", "SingleCellExperiment", "AUCell", "mixtools")
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
library(celldex)
library(SingleR)
library(SingleCellExperiment)
library(jsonlite)
library(AUCell)

args <- commandArgs(trailingOnly=TRUE)
seurat_object <- args[1]
sample_name <- args[2]
GSEA_hallmarks_json_path <- args[3]

seurat_object <- readRDS(seurat_object)


######## SingleR ########

monaco.ref <- celldex::MonacoImmuneData()
if (sample_name == "Integrated") {
  sce <- as.SingleCellExperiment(DietSeurat(seurat_object, assays = "integrated"))
} else {
  sce <- as.SingleCellExperiment(DietSeurat(seurat_object))
}
predictions <- SingleR(test=sce, assay.type.test=1, 
                       ref=monaco.ref, labels=monaco.ref$label.main)
rm(sce)
rm(monaco.ref)
seurat_object@meta.data$predictions <- predictions$pruned.labels
seurat_object <- SetIdent(seurat_object, value = "monaco.main")
predictions_plot <- DimPlot(seurat_object, label = TRUE, repel = TRUE, 
                            label.size = 3, group.by = "predictions") +
                            ggtitle("SingleR Predictions")
ggsave(filename = "predictions_plot.png", plot = predictions_plot, dpi = 300, height=5, units = "in")

count_table <- table(predictions$pruned.labels)
cell_type_counts <- as.data.frame(count_table)
colnames(cell_type_counts) <- c("Cell_Type", "Count")
total_cells <- sum(cell_type_counts$Count)
cell_type_counts$Proportion <- cell_type_counts$Count / total_cells * 100
cell_type_counts$Group <- sample_name
write.table(cell_type_counts, "cell_type_counts.txt", sep="\t", quote = F, row.names = F)
rm(predictions)
annotations <- ggplot(cell_type_counts, aes(fill=Cell_Type, y=Group, x=Proportion)) + 
  geom_bar(position="fill", stat="identity", colour = "black", width = 0.1) +
  scale_fill_brewer(palette = "Pastel1") +
  scale_x_continuous(labels = scales::percent_format(scale = 100)) +
  ylab(NULL) +
  xlab("Proportion (%)") +
  ggtitle("Fraction of cells types") +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black")) +
  labs(fill = "Cell Type") 
ggsave(filename = "annotations.png", plot = annotations, dpi = 300, width = 6, 
       height = 4, units = "in")

######### Hallmaks of cancer #########

hallmarks <- c('HALLMARK_TNFA_SIGNALING_VIA_NFKB',
 'HALLMARK_HYPOXIA',
 'HALLMARK_CHOLESTEROL_HOMEOSTASIS',
 'HALLMARK_MITOTIC_SPINDLE',
 'HALLMARK_WNT_BETA_CATENIN_SIGNALING',
 'HALLMARK_TGF_BETA_SIGNALING',
 'HALLMARK_IL6_JAK_STAT3_SIGNALING',
 'HALLMARK_DNA_REPAIR',
 'HALLMARK_G2M_CHECKPOINT',
 'HALLMARK_APOPTOSIS',
 'HALLMARK_NOTCH_SIGNALING',
 'HALLMARK_ADIPOGENESIS',
 'HALLMARK_ANDROGEN_RESPONSE',
 'HALLMARK_MYOGENESIS',
 'HALLMARK_PROTEIN_SECRETION',
 'HALLMARK_INTERFERON_ALPHA_RESPONSE',
 'HALLMARK_INTERFERON_GAMMA_RESPONSE',
 'HALLMARK_APICAL_JUNCTION',
 'HALLMARK_APICAL_SURFACE',
 'HALLMARK_HEDGEHOG_SIGNALING',
 'HALLMARK_COMPLEMENT',
 'HALLMARK_UNFOLDED_PROTEIN_RESPONSE',
 'HALLMARK_PI3K_AKT_MTOR_SIGNALING',
 'HALLMARK_E2F_TARGETS',
 'HALLMARK_MYC_TARGETS_V1',
 'HALLMARK_MYC_TARGETS_V2',
 'HALLMARK_INFLAMMATORY_RESPONSE',
 'HALLMARK_XENOBIOTIC_METABOLISM',
 'HALLMARK_FATTY_ACID_METABOLISM',
 'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
 'HALLMARK_GLYCOLYSIS',
 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',
 'HALLMARK_P53_PATHWAY',
 'HALLMARK_UV_RESPONSE_UP',
 'HALLMARK_UV_RESPONSE_DN',
 'HALLMARK_ANGIOGENESIS',
 'HALLMARK_HEME_METABOLISM',
 'HALLMARK_COAGULATION',
 'HALLMARK_IL2_STAT5_SIGNALING',
 'HALLMARK_BILE_ACID_METABOLISM',
 'HALLMARK_PEROXISOME',
 'HALLMARK_ALLOGRAFT_REJECTION',
 'HALLMARK_SPERMATOGENESIS',
 'HALLMARK_KRAS_SIGNALING_UP',
 'HALLMARK_KRAS_SIGNALING_DN',
 'HALLMARK_PANCREAS_BETA_CELLS',
 'HALLMARK_ESTROGEN_RESPONSE_EARLY',
 'HALLMARK_ESTROGEN_RESPONSE_LATE',
 'HALLMARK_MTORC1_SIGNALING')
gsea_hallmarks <- fromJSON(GSEA_hallmarks_json_path)
genes_list <- list()

for (category in hallmarks) {
  genes_list[[category]] <- gsea_hallmarks[[category]]$geneSymbols
}

genes_in_seurat <- rownames(seurat_object)
filtered_genes_list <- list()
for (category in names(genes_list)) {
  genes_present <- intersect(genes_list[[category]], genes_in_seurat)
  percentage_present <- length(genes_present) / length(genes_list[[category]]) * 100
  
  if (percentage_present >= 20) {
    filtered_genes_list[[category]] <- genes_list[[category]]
  }
}
rm(genes_list)

counts <- GetAssayData(seurat_object, layer = "data")
ranking <- AUCell_buildRankings(counts)

if (!"hmc" %in% colnames(seurat_object@meta.data)) {
  seurat_object$hmc <- NA
}

for (category in names(filtered_genes_list)) {
  cell_AUC <- AUCell_calcAUC(filtered_genes_list[[category]], ranking)
  cell_assigment <- AUCell_exploreThresholds(cell_AUC, plotHist = TRUE, assign=TRUE)
  seurat_object$hmc <- ifelse(colnames(seurat_object) %in% cell_assigment$geneSet$assignment, 
                              gsub("_", " ", sub("^HALLMARK_", "", category)), 
                              seurat_object$hmc)
}

hmc_plot <- DimPlot(seurat_object, label = TRUE, repel = TRUE, 
                            label.size = 3, group.by = "hmc") +
                            ggtitle("Hallmarks")
ggsave(filename = "hmc_plot.png", plot = hmc_plot, dpi = 300, height=5, units = "in")


######### EMT #########

genes_present <- intersect(gsea_hallmarks$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION$geneSymbols, genes_in_seurat)
percentage_present <- length(genes_present) / length(gsea_hallmarks$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION$geneSymbols) * 100
if (percentage_present >= 20) {
  cell_AUC <- AUCell_calcAUC(gsea_hallmarks$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION$geneSymbols, ranking)
  cell_assigment <- AUCell_exploreThresholds(cell_AUC, plotHist = TRUE, assign=TRUE)
  seurat_object$EMT <- ifelse(colnames(seurat_object) %in% cell_assigment$geneSet$assignment, 
                              "EMT", 
                              NA)
  emt_plot <- DimPlot(seurat_object, repel = TRUE, 
                            label.size = 3, group.by = "EMT") +
                            ggtitle("EMT")
  ggsave(filename = "emt_plot.png", plot = emt_plot, dpi = 300, height=5, units = "in")
} else {
  print("EMT hallmark not present in the dataset")
}
