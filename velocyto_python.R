#conda activate seurat_v4
setwd("/projects/chicken_hip/expression_matrix_new2/velocyto")
rm(list=ls())
gc()

library(Seurat)
library(Matrix)
library(velocyto.R)
library(SeuratWrappers)

seudt <- readRDS("/projects/chicken_hip/expression_matrix_new2/chicken_hip_scRNA_harmony_annotation_8_samples_final.rds")
names(seudt@meta.data)

seu_cell <- Cells(seudt)
seu_cell <- gsub("_", ":", seu_cell)
seu_cell <- paste0(seu_cell,"x")
head(seu_cell)

seudt <- RenameCells(seudt, new.names = seu_cell)

WLH1 <- as.Seurat(ReadVelocity("WLH-1.loom"))
WLH2 <- as.Seurat(ReadVelocity("WLH-2.loom"))
WLH3 <- as.Seurat(ReadVelocity("WLH-3.loom"))
RJF1 <- as.Seurat(ReadVelocity("RJF-1.loom"))
RJF2 <- as.Seurat(ReadVelocity("RJF-2.loom"))
RJF3 <- as.Seurat(ReadVelocity("RJF-3.loom"))
WLH4 <- as.Seurat(ReadVelocity("WLH-4.loom"))
WLH5 <- as.Seurat(ReadVelocity("WLH-5.loom"))

all_loom <- merge(x = WLH1, y = c(WLH2, WLH3, RJF1, RJF2), merge.data = TRUE)
dim(all_loom)

all(colnames(seudt) %in% colnames(all_loom))
all(rownames(seudt) %in% rownames(all_loom))

write.csv(seu_cell, file = "cellID_8_samples_velocyto.csv", row.names = FALSE)

seu_umap <- Embeddings(seudt, reduction = "umap")
write.csv(seu_umap, file = "umap_8_samples_velocyto.csv")

seu_subtype <- as.character(seudt$Subtype)
names(seu_subtype) <- colnames(seudt)
write.csv(seu_subtype, file = "subtype_8_samples_velocyto.csv")

#conda activate velocyto
#python

import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
import os
import scanpy as sc

os.chdir("/projects/chicken_hip/expression_matrix_new2/velocyto") 
print(os.getcwd())

WLH1 = anndata.read_loom("WLH-1.loom")
WLH1
print(WLH1.var.columns)

WLH2 = anndata.read_loom("WLH-2.loom")
WLH2

WLH3 = anndata.read_loom("WLH-3.loom")
WLH3

WLH4 = anndata.read_loom("WLH-4.loom",var_names='Accession')
WLH4
print(WLH4.var.columns)

WLH5 = anndata.read_loom("WLH-5.loom")
WLH5

RJF1 = anndata.read_loom("RJF-1.loom")
RJF1

RJF2 = anndata.read_loom("RJF-2.loom")
RJF2

RJF3 = anndata.read_loom("RJF-3.loom")
RJF3

seu_cell = pd.read_csv("cellID_8_samples_velocyto.csv")
seu_umap = pd.read_csv("umap_8_samples_velocyto.csv")
seu_subtype = pd.read_csv("subtype_8_samples_velocyto.csv")

WLH1 = WLH1[np.isin(WLH1.obs.index,seu_cell["x"])]
WLH2 = WLH2[np.isin(WLH2.obs.index,seu_cell["x"])]
WLH3 = WLH3[np.isin(WLH3.obs.index,seu_cell["x"])]
WLH4 = WLH4[np.isin(WLH4.obs.index,seu_cell["x"])]
WLH5 = WLH5[np.isin(WLH5.obs.index,seu_cell["x"])]
RJF1 = RJF1[np.isin(RJF1.obs.index,seu_cell["x"])]
RJF2 = RJF2[np.isin(RJF2.obs.index,seu_cell["x"])]
RJF3 = RJF3[np.isin(RJF3.obs.index,seu_cell["x"])]

print(WLH1.obs_names[:10]) 
print(WLH1.var_names[:10]) 

print(RJF3.obs_names[:10]) 
print(RJF3.var_names[:10])

dup_names = WLH4.var_names[WLH4.var_names.duplicated()]
print(f"duplicated var_names Number: {len(dup_names)}")

adata = WLH1.concatenate(WLH2,WLH3,WLH4,WLH5,RJF1,RJF2,RJF3)

adata_index = pd.DataFrame(adata.obs.index)
adata_index = adata_index.rename(columns = {0:'Cell ID'})
adata_index = adata_index.rename(columns = {"CellID":'Cell ID'})

rep=lambda x : x.rsplit("-", 1)[0]
adata_index["Cell ID"]=adata_index["Cell ID"].apply(rep)

seu_umap = seu_umap.rename(columns = {'Unnamed: 0':'Cell ID'})
seu_umap = seu_umap[np.isin(seu_umap["Cell ID"],adata_index["Cell ID"])]
seu_umap=seu_umap.drop_duplicates(subset=["Cell ID"])
umap_ordered = adata_index.merge(seu_umap, on = "Cell ID")
umap_ordered = umap_ordered.iloc[:,1:]
adata.obsm['X_umap'] = umap_ordered.values


seu_subtype = seu_subtype.rename(columns = {'Unnamed: 0':'Cell ID'})
seu_subtype = seu_subtype[np.isin(seu_subtype["Cell ID"],adata_index["Cell ID"])]
seu_subtype=seu_subtype.drop_duplicates(subset=["Cell ID"])
seu_subtype_ordered = adata_index.merge(seu_subtype, on = "Cell ID")
seu_subtype_ordered = seu_subtype_ordered.iloc[:,1:]
adata.obs['Subtype'] = seu_subtype_ordered["x"].values

scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode = "stochastic")
scv.tl.velocity_graph(adata)

adata.write('RNA_velocyto.h5ad', compression='gzip')

