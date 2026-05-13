#conda activate seurat_v4
setwd("/projects/chicken_hip_yunmei/expression_matrix_new2")
rm(list=ls())
gc()

library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(future)
options(future.globals.maxSize = 1e9)

anno <- read.table("/projects/chicken_hip/GRCg7b.bioMart.v111.csv", sep=",", header = TRUE, stringsAsFactors = FALSE)
anno <- anno[, c(1:2,3,6)]
colnames(anno) <- c("EnsID", "GeneName", "Chr", "GeneType")
mtGene <- unique(anno[anno$Chr=="MT",]$EnsID)

root_path <- '/projects/chicken_hip/expression_matrix_new2'
subdirectories <- list.dirs(path = root_path, recursive = FALSE)
file_names <- "filtered_cell_gene_matrix"
file_paths <- file.path(subdirectories, "outs", file_names)

#generate Seurat objects
sceList = lapply(file_paths,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder,gene.column = 1), 
                     project = folder,min.cells = 50, min.features = 200)
})

head(rownames(sceList[[1]][["RNA"]]@counts))

le_names <- character(length(sceList))
for(i in seq_along(sceList)) {
  sample_name <- basename(dirname(dirname(file_paths[i])))
  sample_names[i] <- sample_name
  sceList[[i]]$orig.ident <- sample_name
  Idents(sceList[[i]]) <- sample_name
}

#Calculate nCount, nFeature, and the proportion of mitochondrial gene<
s- lapply(sceList, function(obj) {
  mt.genes <- intersect(rownames(obj), mtGene)
  obj[["percent_mt"]] <- PercentageFeatureSet(obj, features = mt.genes)
  obj$nCount_nFeature <- obj$nCount_RNA / obj$nFeature_RNA
  obj
})

before_ncells <- sapply(sceList, ncol)
names(before_ncells) <- sample_names
cat("before_ncells")

sceList <- lapply(X = sceList, FUN = function(x){
  x <- subset(x,
              subset = nFeature_RNA >= 200 & nFeature_RNA <= 7500 & 
                percent_mt <= 5 & nCount_nFeature >= 1.5 & nCount_nFeature <= 10 & 
                nCount_RNA >= 1000 & nCount_RNA <= 60000)})

after_ncells <- sapply(sceList, ncol)
names(after_ncells) <- sample_names
cat("after_ncells)

head(sceList[[1]]@meta.data)

#Using DoubletFinder to remove doublets
library(DoubletFinder)
for(i in 1:length(sceList)){
  sc <- sceList[[i]]
  sc <- NormalizeData(sc, verbose = FALSE)
  sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sc <- ScaleData(sc, verbose = FALSE)
  sc <- RunPCA(sc)
  sc <- RunUMAP(sc, reduction = "pca", dims = 1:30, verbose = FALSE)

  ##pK Identification (no ground-truth)
  sweep_res.sc <- paramSweep(sc, PCs = 1:30, sct = FALSE)
  sweep_stats.sc <- summarizeSweep(sweep_res.sc, GT = FALSE)
  bcmvn.sc <- find.pK(sweep_stats.sc)
  mpK <- as.numeric(as.vector(bcmvn.sc$pK[which.max(bcmvn.sc$BCmetric)]))
  print("mpK")

  ##Homotypic Doublet Proportion Estimate
  annotation.sc <- sc@meta.data$seurat_clusters
  nExp_poi.sc <- round(0.06*nrow(sc@meta.data))  ## Assuming 6% doublet formation rate - tailor for your dataset
  print("nExp_poi.sc")

  ##Run DoubletFinder with varying classification stringencies
  sc <- doubletFinder(sc, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.sc, sct = FALSE)
  sc@meta.data$dbnm <- sc@meta.data[, grep("DF.classifications_", colnames(sc@meta.data))]
  sc@meta.data$pANN <- sc@meta.data[, grep("pANN_", colnames(sc@meta.data))]
  print("table dbnm")
  table(sc@meta.data[, "dbnm"])
  sc_singlet <- subset(sc, subset = dbnm == "Singlet")
  sceList[[i]] <- sc_singlet
}

after_ncells <- sapply(sceList, ncol)
names(after_ncells) <- sample_names
cat("after_ncells")

#merge
all_list <- list()
meta_list <- list()
sample_names <- character(length(sceList))

for(i in seq_along(sceList)){
    sc <- sceList[[i]]
    DefaultAssay(sc) <- "RNA"
    meta <- sc@meta.data[, "orig.ident", drop = FALSE]
    meta_list[[i]] <- meta
    sample_names[i] <- unique(sc$orig.ident)[1]
    count_tmp <- GetAssayData(sc, slot = "counts")
    seu_tmp <- CreateSeuratObject(counts = count_tmp, meta.data = meta)
    all_list[[i]] <- seu_tmp
}

all_merge <- merge(x = all_list[[1]], 
                       y = all_list[2:length(all_list)], 
                       add.cell.ids = sample_names, 
                       project = "chick_hip")
#save data
names(all_merge@meta.data)
table(all_merge@meta.data$orig.ident)
saveRDS(all_merge, file ="chicken_hip_scRNA_qc_merged_8_samples.rds")
