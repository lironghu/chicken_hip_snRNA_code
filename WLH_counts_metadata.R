setwd("/projects/chicken_hip/expression_matrix_new2/pyscenic")
rm(list=ls())
gc()

library(dplyr)
library(Seurat)

#-----------------------------------------------------------one_to_one_orthologs
orth <- read.csv("/projects/hulr/0.snRNA/05.one_to_one_orthologs/AAA.One_to_One_Orthologs.RefMouse_modified2.csv")

#-----------------------------------------------------------matrix
data <- readRDS("/projects/chicken_hip/expression_matrix_new2/chicken_hip_scRNA_harmony_annotation_8samples_final.rds")

Idents(data) <- "Group"
data<- subset(data,idents="WLH")
names(data@meta.data)

counts <- GetAssayData(object = data, 
                             assay = "RNA", 
                             layer= "counts")
head(rownames(counts))
head(orth)

counts<- as.data.frame(counts)
counts$gene_id_chicken <- rownames(counts)
chicken_mouse_gene <- orth[,c("gene_id_chicken","gene_name_mouse")]
counts <- inner_join(counts,chicken_mouse_gene,by="gene_id_chicken")
dim(counts)

counts <- counts[!is.na(counts$gene_id_chicken) & 
                              counts$gene_id_chicken != "#N/A", ]
dim(counts)

rownames(counts) <- counts$gene_name_mouse
counts$gene_id_chicken <- NULL
counts$gene_name_mouse <- NULL
head(rownames(counts))

counts<- t(as.matrix(counts))
counts[1:4,1:4]
write.csv(counts,file = "WLH_counts.csv")

metadata <- data@meta.data
write.csv(metadata, file = "WLH_metadata.csv", row.names = TRUE)