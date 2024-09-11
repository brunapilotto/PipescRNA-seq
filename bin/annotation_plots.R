#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(Cairo)
library(celldex)
library(SingleR)
library(SingleCellExperiment)
library(jsonlite)

args <- commandArgs(trailingOnly=TRUE)
seurat_object <- args[1]
sample_name <- args[2]

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
ggsave(filename = "predictions_plot.png", plot = predictions_plot, dpi = 300, height=6, width=9, units = "in")

count_table <- table(predictions$pruned.labels)
cell_type_counts <- as.data.frame(count_table)
colnames(cell_type_counts) <- c("Cell_Type", "Count")
total_cells <- sum(cell_type_counts$Count)
cell_type_counts$Proportion <- cell_type_counts$Count / total_cells * 100
cell_type_counts$Group <- sample_name
write.table(cell_type_counts, "cell_type_counts.txt", sep="\t", quote = F, row.names = F)

all.markers <- metadata(predictions)$de.genes
rm(predictions)
top_genes_list <- list()
for (cell_type in cell_type_counts$Cell_Type) {
  unique_genes <- unique(unlist(all.markers[[cell_type]]))
  top_genes <- head(unique_genes, 200)
  top_genes_list[[cell_type]] <- top_genes
}
write_json(top_genes_list, "top_genes_by_cell_type.json")
rm(all.markers)

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
