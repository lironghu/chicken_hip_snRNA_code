setwd("/projects/chicken_hip/expression_matrix_new2/miloR")
rm(list=ls())
gc()

library(Seurat)
library(ggplot2)
library(miloR)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(BiocParallel)
register(MulticoreParam(workers = 8, progressbar = TRUE)) 

seurat_obj <- readRDS("chicken_hip_scRNA_harmony_annotation_8_samples_final.rds")
names(seurat_obj@meta.data)
table(seurat_obj@meta.data$Celltype)

#miloR
scRNA_pre <- as.SingleCellExperiment(seurat_obj)
scRNA_milo <- Milo(scRNA_pre)

if("UMAP" %in% reducedDimNames(scRNA_pre)) {
  reducedDim(scRNA_milo, "UMAP") <- reducedDim(scRNA_pre, "UMAP")
} else if("umap" %in% reducedDimNames(scRNA_pre)) {
  reducedDim(scRNA_milo, "umap") <- reducedDim(scRNA_pre, "umap")
}

scRNA_milo
scRNA_milo <- buildGraph(scRNA_milo, k = 15, d = 15,reduced.dim = "PCA")
scRNA_milo <- makeNhoods(scRNA_milo, prop = 0.2, k = 15, d = 15, 
                    refined = TRUE, reduced_dims = "PCA")

plot <- plotNhoodSizeHist(scRNA_milo)
ggsave("miloR_nhood_size_histogram.pdf", plot, width = 8, height = 6)

scRNA_milo <- countCells(scRNA_milo, 
                    meta.data = as.data.frame(colData(scRNA_milo)),
                    sample="Sample") 

traj_design <- data.frame(colData(scRNA_milo))[,c("Sample", "Group")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$Sample
traj_design <- traj_design[colnames(nhoodCounts(scRNA_milo)), , drop=FALSE]
traj_design$Group<- factor(traj_design$Group, level = c("WLH", "RJF"))
traj_design

scRNA_milo <- calcNhoodDistance(scRNA_milo, d=15, reduced.dim = "PCA")
da_results <- testNhoods(scRNA_milo,design = ~ Group, design.df = traj_design)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 

scRNA_milo <- buildNhoodGraph(scRNA_milo)

#Plot single-cell UMAP
umap_1 <- plotReducedDim(scRNA_milo, dimred = "UMAP", 
                          colour_by="Group", text_by = "Celltype", 
                          text_size = 3, point_size=0.5) + guides(fill="none")
umap_2 <- plotReducedDim(scRNA_milo, dimred = "UMAP", 
                          colour_by="Sample", text_by = "Celltype", 
                          text_size = 3, point_size=0.5) + guides(fill="none")

#Plot neighbourhood graph
nh_graph <- plotNhoodGraphDA(scRNA_milo, da_results,
                                layout="UMAP",alpha=0.9)+
  theme(legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 8))   

da_results <- annotateNhoods(scRNA_milo, 
                             da_results, 
                             coldata_col = "Celltype")
str(da_results)
da_results$Celltype <- ifelse(da_results$Celltype == "nIPC-like-micro", "mIPC",
                             ifelse(da_results$Celltype == "nIPC-like-Tcell", "tIPC",
                                    da_results$Celltype))
table(da_results$Celltype)
range(da_results$SpatialFDR)

da_plot <- plotDAbeeswarm(da_results, group.by = "Celltype",alpha = 0.9)+
  labs(title = "RJF/WLH",x = NULL)+ 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 14, face = "bold"))

pdf("miloR_results_all_cell.pdf", width = 6, height = 6)
print(nh_graph)
print(da_plot)
dev.off()

