#conda activate seurat_v4
setwd("/projects/chicken_hip/expression_matrix_new2/deg_and_enrichment")
rm(list=ls())
gc()

library(Seurat)
library(ggplot2)
library(dplyr)
library(MAST)
library(future)
options(future.globals.maxSize = 1e9)

anno <- read.table("/projects/hulr/chicken_hip_yunmei/GRCg7b.bioMart.v111.csv", sep=",", header = TRUE, stringsAsFactors = FALSE)
anno <- anno[, c(1:2,3,6)]
colnames(anno) <- c("EnsID", "GeneName", "Chr", "GeneType")
mtGene <- unique(anno[anno$Chr=="MT",]$EnsID)

scRNA <- readRDS("chicken_hip_scRNA_harmony_annotation_8_samples_final.rds")
names(scRNA@meta.data)
table(scRNA@meta.data$Sample)

dim(scRNA)
scRNA <- subset(scRNA, features = setdiff(rownames(scRNA), mtGene))
dim(scRNA)

celltype_all <- sort(unique(scRNA$Celltype))
deg_results_celltype <- list()

for (cl in celltype_all) {
  cells_cl <- WhichCells(scRNA, idents = cl)
  if (length(cells_cl) < 3) next
  Idents(scRNA) <- "Celltype"
  obj_cl <- subset(scRNA, cells = cells_cl)
  
  markers <- FindMarkers(obj_cl,
                         ident.1 = "WLH",
                         ident.2 = "RJF",
                         group.by = "Group",
                         min.pct = 0.1,test.use = "MAST",latent.vars = "Batch",
                         only.pos = FALSE,
                         p.adjust.method = "BH")
  
  markers$celltype_all <- cl
  markers$gene <- rownames(markers)
  markers$higher_in <- ifelse(markers$avg_log2FC > 0, "WLH", "RJF")
  
  deg_results_celltype[[as.character(cl)]] <- markers
}

deg_combined_celltype <- bind_rows(deg_results_celltype, .id = "celltype")
head(deg_combined_celltype)
colnames(deg_combined_celltype)[colnames(deg_combined_celltype) == "gene"] <- "GeneID"
head(deg_combined_celltype)

anno <- read.table("/projects/chicken_hip/GRCg7b.bioMart.v111.csv", sep=",", header = TRUE, stringsAsFactors = FALSE)
anno <- anno[, c(1,3,6,7)]
colnames(anno) <- c("GeneID","Chr", "GeneType","GeneName")
anno$GeneName[is.na(anno$GeneName) | anno$GeneName == ""] <- anno$GeneID[is.na(anno$GeneName) | anno$GeneName == ""]

deg_combined_celltype <- deg_combined_celltype %>%
  left_join(anno %>% select(GeneID, Chr, GeneType, GeneName), by = "GeneID")

write.csv(deg_combined_celltype, "DEG_all_celltype_WLH_vs_RJF_MAST.csv", row.names = FALSE)
