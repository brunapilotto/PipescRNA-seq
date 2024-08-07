#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

args <- commandArgs(trailingOnly=TRUE)
seurat_object <- args[1]
sample_name <- args[2]

retinoblastoma <- c("MKI67", "CDK1", "TOP2A", "KIF14", "UBE2C", "CDC25C")
rod_precursor <- c("NRL", "RCVRN")
retinoma <- c("CDCA7", "HELLS", "MCM3", "PCNA")
rods <- c("CNGA1", "GNAT1", "NR2E3", "PDE6A", "PDE6G", "RHO")
bipolar_cells <- c("CA10", "LRTM1", "PCP2", "PRKCA", "TRPM1", "VSX1", "VSX2")
muller_glia <- c("MIR9-1HG", "CLU", "GLUL", "RLBP1", "SPP1", "VIM")
microglia <- c("AIF1", "C1QA", "HLA-DPA1", "HLA-DRA", "PTPRC")
cones <- c("RXRG", "CRX", "THRB", "ARR3", "GNB3", "OPN1LW", "PDE6H", "GNGT2")
lnc_1 <- c("LINC00152", "LINC00115", "LINC00858", "LINC00202", "LINC00324",
         "LINC00488", "PVT1", "CASC9", "HEIH", "HCP5", 
         "XIST", "FTX", "HOTTIP", "ROR", "DANCR", 
         "THOR", "PANDAR", "PlncRNA-1", "BANCR", "CCAT1", 
         "MIAT", "HIF1A‐AS1", "ANRIL", "TP73‐AS1", "UCA1")
lnc_2 <- c("TUG1", "MALAT1", "FOXD2‐AS1", "PROX1‐AS1", "ELFN1‐AS1", 
         "ADPGK‐AS1", "ZNRD1‐AS1", "ILF3‐AS1", "AFAP1-AS1", "LEF1‐AS1", 
         "TMPO‐AS1", "TRPM2‐AS", "HOXA11‐AS", "SNHG14", "SNHG16", 
         "SNHG20", "ZFPM2‐AS1", "FEZF1‐AS1", "MIR17HG", "MIR7‐3HG", 
         "HOTAIR", "KCNQ1OT1", "SND1-IT1", "NEAT1")
Cell_types <- c("cones", "lnc_1", "lnc_2", "retinoblastoma")

seurat_object <- readRDS(seurat_object)

genes_to_remove <- c()
for (Cell in Cell_types) {
  for (gene in get(Cell)) {
    if (!(gene %in% rownames(seurat_object))) {
      message(paste("Gene not found: ", gene))
      genes_to_remove <- c(genes_to_remove, gene)
    }
  }
}

for(Cell in Cell_types){
  genes <- setdiff(get(Cell), genes_to_remove)
  plot <- VlnPlot(seurat_object, features = genes, pt.size = 0)
  ggsave(filename = paste0(Cell, "_", sample_name, "_vln_plot.png"), plot = plot, dpi = 300, width = 15, height = 8, units = "in")
}
