#conda activate seurat_v4
setwd("/projects/chicken_hip/expression_matrix_new2")
rm(list=ls())
gc()

library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(harmony)
library(scCustomize)
library(future)
options(future.globals.maxSize = 1e9)

merged <- readRDS("chicken_hip_scRNA_qc_merged_8_samples.rds")
merged <- merged@meta.data
merged_count <- GetAssayData(object = merged, slot = "counts")
merged <- CreateSeuratObject(counts = merged_count, meta.data = meta_5)
names(merged@meta.data)

group <- as.character(merged@meta.data$orig.ident)
group <- gsub("-[0-9]", "", group)
sample <- as.character(merged@meta.data$orig.ident)

merged <- AddMetaData(merged, group, col.name = "Group")
merged <- AddMetaData(merged, sample, col.name = "Sample")
head(merged@meta.data)

merged$Batch <- case_when(
  merged$orig.ident == "WLH-1" ~ "Batch1",
  merged$orig.ident == "WLH-2" ~ "Batch1",
  merged$orig.ident == "WLH-3" ~ "Batch1",
  merged$orig.ident == "RJF-1" ~ "Batch1",
  merged$orig.ident == "RJF-2" ~ "Batch1",
  merged$orig.ident == "RJF-3" ~ "Batch2",
  merged$orig.ident == "WLH-4" ~ "Batch2",
  merged$orig.ident == "WLH-5" ~ "Batch2",
)
names(merged@meta.data)
table(merged$Batch)

#harmony
orth <- read.csv("/projects/chicken_hip/chicken_human_Homologues_biomart.v115.csv", 
                 sep = ",",
                 header = TRUE, 
                 stringsAsFactors = FALSE, 
                 quote = "")
head(orth)

s.genes <- orth[toupper(orth$Human_gene_name) %in% toupper(cc.genes$s.genes), "Gene_stable_ID"]
g2m.genes <- orth[toupper(orth$Human_gene_name) %in% toupper(cc.genes$g2m.genes), "Gene_stable_ID"]

exp_s <- intersect(s.genes, rownames(merged))
exp_g2m <- intersect(g2m.genes, rownames(merged))

merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
merged <- CellCycleScoring(merged, s.features = exp_s, g2m.features = exp_g2m)
merged$CC.Difference <- merged$S.Score - merged$G2M.Score
merged <- ScaleData(merged,vars.to.regress = c("nCount_RNA", "percent_mt", "CC.Difference"))
merged <- RunPCA(merged, verbose = FALSE)
merged <- RunHarmony(merged, group.by.vars = "Sample")

merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30)
merged <- FindClusters(merged, resolution = seq(from = 0.5, to = 4, by = 0.5))
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30)

pdf("Seurat_Analysis_Plots_8_samples.pdf", width = 20, height = 16)
DimPlot(merged, reduction = "umap", group.by = "Group")
DimPlot(merged, reduction = "umap", group.by = "Sample")
DimPlot(merged, reduction = "umap", group.by = "Celltype_batch2")
DimPlot(merged, reduction = "umap", group.by = "Subtype_batch1")
DimPlot(merged, reduction = "umap", group.by = "Celltype_batch1")
DimPlot(merged, group.by = "RNA_snn_res.1.5", label = TRUE) + NoLegend() + ggtitle("Res 1.5")
DimPlot(merged, group.by = "RNA_snn_res.2", label = TRUE) + NoLegend() + ggtitle("Res 2")
DimPlot(merged, group.by = "RNA_snn_res.2.5", label = TRUE) + NoLegend() + ggtitle("Res 2.5")
DimPlot(merged, group.by = "RNA_snn_res.3", label = TRUE) + NoLegend() + ggtitle("Res 3")
DimPlot(merged, group.by = "RNA_snn_res.3.5", label = TRUE) + NoLegend() + ggtitle("Res 3.5")
DimPlot(merged, group.by = "RNA_snn_res.4", label = TRUE) + NoLegend() + ggtitle("Res 4")
FeaturePlot(merged, features = c("S.Score", "G2M.Score"), cols = c("#0066CC", "#FFFF00"), reduction="umap", label = FALSE, order = TRUE) + theme(plot.title = element_text(size=8))
genes_to_check <- c("ENSGALG00010022483","ENSGALG00010002750", #NPC#TOP2A#MKI67
                   "ENSGALG00010028070","ENSGALG00010016519","ENSGALG00010015393",#Pre#MEX3A#MEX3B#DCX
                   "ENSGALG00010023237","ENSGALG00010005413",#InN#GAD1#GAD2
                   "ENSGALG00010020328",#ExN#SLC17A6
                   "ENSGALG00010005544","ENSGALG00010012939","ENSGALG00010021974",#Astro#AQP4#GJA1#GFAP
                   "ENSGALG00010029758","ENSGALG00010020955",#Ependymal#DNAH9#FHAD1
                   "ENSGALG00010016804","ENSGALG00010004746",#Microglial#CSF1R#LYN
                   "ENSGALG00010023715","ENSGALG00010023743",#Tcell#CD3E#CD3D
                   "ENSGALG00010016496","ENSGALG00010007781",#Olig#MBP#PLP1
                   "ENSGALG00010013621","ENSGALG00010009400",#OPC#PDGFRA#CNTFR
                   "ENSGALG00010020959","ENSGALG00010024064",#cPOC#ZNF488#GPR17
                   "ENSGALG00010005282","ENSGALG00010024534" #Vascular#FLT1#VWF
)
custom_colors <- colorRampPalette(c( "#eeeeee", "#3954a5"))(100)
DotPlot(merged, 
        features = genes_to_check, 
        assay = 'RNA', 
        scale = TRUE,
        col.min = 0,
        group.by = "RNA_snn_res.4",
        cluster.id = TRUE
) +  
  coord_flip() +
  ggtitle("") +
  scale_color_gradientn(colors = custom_colors) +
  theme(
    panel.border = element_rect(color = "black", linewidth = 1)
  )
dev.off()

saveRDS(merged, "chicken_hip_harmony_8_samples.rds")

#marker gene
names(merged@meta.data)
Idents(merged) <- "RNA_snn_res.4"
DefaultAssay(merged) <- "RNA"

allMarkers.res4 <- FindAllMarkers(object = merged, min.cells.group = 1, min.pct = 0.25, only.pos = TRUE)
head(allMarkers.res4)
write.csv(allMarkers.res4, file = "allMarkers.res4_8_samples.csv", row.names = FALSE, quote = FALSE)

sigMarkers.res4 <- allMarkers.res4 %>%
  dplyr::select(gene, everything()) %>%
  subset(p_val_adj < 0.05 & abs(avg_log2FC) > 0.585)

top50 <- sigMarkers.res4 %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

pdf("top50_dotplot_per_cluster_res4_8_samples.pdf", width = 18, height = 14)
for (cl in sort(unique(top50$cluster))) {
  genes_cl <- top50 %>% filter(cluster == cl) %>% pull(gene)
  genes_cl <- genes_cl[genes_cl %in% rownames(merged)]
  
  if(length(genes_cl) > 0) {
    print(
      DotPlot(merged,
              features = genes_cl,
              group.by = "RNA_snn_res.4",
              cols = c("lightgrey", "red"),
              dot.scale = 4,
              scale = TRUE,
              col.min = -2.5,
              col.max = 2.5) +
        ggtitle(paste0("Cluster ", cl, "top 50 markers")) +
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              axis.text.x = element_text(angle = 45, hjust = 1)) +
        ylab("Genes") +
        xlab("Clusters")
    )
  } else {
    print(paste("Cluster", cl, "NA"))
  }
}
dev.off()

#cell type annotation
merged$Celltype <- case_when(
 merged$RNA_snn_res.4 %in% c(1,9,64) ~ "OPC",
 merged$RNA_snn_res.4 %in% c(58) ~ "cOPC",
 merged$RNA_snn_res.4 %in% c(2,3,16,14,26,41,63) ~ "Oligo",
 merged$RNA_snn_res.4 %in% c(0,11,32,50,65) ~ "Astro",
 merged$RNA_snn_res.4 %in% c(10,35,52,54) ~ "Epen",
 merged$RNA_snn_res.4 %in% c(34,49,17) ~ "VLMC",
 merged$RNA_snn_res.4 %in% c(46) ~ "Endo",
 merged$RNA_snn_res.4 %in% c(4,43,44,59) ~ "Micro",
 merged$RNA_snn_res.4 %in% c(22) ~ "Tcell",
 merged$RNA_snn_res.4 %in% c(61) ~ "IPC1",
 merged$RNA_snn_res.4 %in% c(57) ~ "IPC2",
  merged$RNA_snn_res.4 %in% c(53) ~ "Pre",
  merged$RNA_snn_res.4 %in% c(19,31,33,42,56,38,7,15,20,18,40,48,51,27,24,25,47,8,13,12,6,5,28) ~ "ExN",
  merged$RNA_snn_res.4 %in% c(45,21,36,37,23,55,29,30,39) ~ "InN",
  TRUE ~ "Misc")
names(merged@meta.data)
table(merged$Celltype)

merged$Subtype <- case_when(
 merged$RNA_snn_res.4 %in% c(1,9,64) ~ "OPC",
 merged$RNA_snn_res.4 %in% c(58) ~ "cOPC",
 merged$RNA_snn_res.4 %in% c(2,3,16,14,26,41,63) ~ "Oligo",
 merged$RNA_snn_res.4 %in% c(0,11,32,50,65) ~ "Astro",
 merged$RNA_snn_res.4 %in% c(10,35,52,54) ~ "Epen",
 merged$RNA_snn_res.4 %in% c(17) ~ "VLMC1",
 merged$RNA_snn_res.4 %in% c(49) ~ "VLMC2",
 merged$RNA_snn_res.4 %in% c(34) ~ "VLMC3",
 merged$RNA_snn_res.4 %in% c(46) ~ "Endo",
 merged$RNA_snn_res.4 %in% c(4,43,44,59) ~ "Micro",
 merged$RNA_snn_res.4 %in% c(22) ~ "Tcell",
 merged$RNA_snn_res.4 %in% c(61) ~ "IPC1",
 merged$RNA_snn_res.4 %in% c(57) ~ "IPC2",
  merged$RNA_snn_res.4 %in% c(53) ~ "Pre",
  merged$RNA_snn_res.4 %in% c(19) ~ "ExN1",
  merged$RNA_snn_res.4 %in% c(31,33,42,56) ~ "ExN2",
  merged$RNA_snn_res.4 %in% c(38,7,15,20,18,40,48,51,27) ~ "ExN3",
  merged$RNA_snn_res.4 %in% c(8,13,24,25,47) ~ "ExN4",
  merged$RNA_snn_res.4 %in% c(12,6,5,28) ~ "ExN5",
  merged$RNA_snn_res.4 %in% c(23,55) ~ "InN1",
  merged$RNA_snn_res.4 %in% c(29,30,39) ~ "InN2",
  merged$RNA_snn_res.4 %in% c(37) ~ "InN3",
  merged$RNA_snn_res.4 %in% c(21,36) ~ "InN4",
  merged$RNA_snn_res.4 %in% c(45) ~ "InN5",
  TRUE ~ "Misc")
names(merged@meta.data)
table(merged$Celltype)

p1 <- DimPlot(merged, 
              reduction = "umap",
              group.by = "Celltype",
              label = TRUE,
              label.size = 5,
              pt.size = 0.8,
              repel = TRUE) +
  ggtitle("UMAP by Subtype") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(merged, 
              reduction = "umap",
              group.by = "Subtype",
              label = TRUE,
              label.size = 5,
              pt.size = 0.8,
              repel = TRUE) +
  ggtitle("UMAP by Subtype") +
  theme(plot.title = element_text(hjust = 0.5))


p3 <- DimPlot(merged, 
              reduction = "umap",
              group.by = "Group", pt.size = 0.8)+
  ggtitle("UMAP by Group") +
  theme(plot.title = element_text(hjust = 0.5))

pdf(file = "UMAP_plots_all_cells_8_samples.pdf", width = 8, height = 6)
print(p1)
print(p2)
print(p3)
dev.off()

#remove Misc and  outlier cells
Idents(merged) <- "Subtype"
merged <- subset(merged, idents = "Misc", invert = TRUE)

umap_coords <- Embeddings(merged, reduction = "umap")
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
meta_subset <- merged@meta.data[, c("Group", "Sample", "RNA_snn_res.4",
                                           "nCount_RNA", "nFeature_RNA", 
                                           "nCount_nFeature", "percent_mt")]

all_meta <- data.frame(
  Cell = rownames(umap_coords),
  umap_coords,
  meta_subset,
  row.names = NULL
)
colnames(all_meta)[colnames(all_meta) == "RNA_snn_res.4"] <- "cluster"
head(all_meta)

rownames(all_meta) <- all_meta$Cell

outlier <- function(clusters, target_col, label, meta) {
  clusters <- setdiff(clusters, c(57, 61))
  all_outer <- data.frame()
  for (cl in clusters) {
    sub_tmp <- meta[meta[[target_col]] == cl, ]
    if (nrow(sub_tmp) == 0) next

    umap1 <- sub_tmp$UMAP_1
    q1 <- quantile(umap1, 0.25)
    q3 <- quantile(umap1, 0.75)
    iqr <- IQR(umap1)
    lower1 <- q1 - 2 * iqr
    upper1 <- q3 + 2 * iqr
    outer1 <- rownames(sub_tmp)[umap1 < lower1 | umap1 > upper1]
    
    umap2 <- sub_tmp$UMAP_2
    q1 <- quantile(umap2, 0.25)
    q3 <- quantile(umap2, 0.75)
    iqr <- IQR(umap2)
    lower2 <- q1 - 2 * iqr
    upper2 <- q3 + 2 * iqr
    outer2 <- rownames(sub_tmp)[umap2 < lower2 | umap2 > upper2]
    
    outer_cells <- union(outer1, outer2)
    if (length(outer_cells) > 0) {
      out_df <- sub_tmp[outer_cells, c("Cell", target_col, "Group", "Sample", 
                                       "UMAP_1", "UMAP_2", "nCount_RNA", 
                                       "nFeature_RNA", "percent_mt", "nCount_nFeature")]
      all_outer <- rbind(all_outer, out_df)
    }
  }
  if (nrow(all_outer) > 0) {
    write.csv(all_outer, paste0("remove_outlier_cells.", label, ".csv"), 
              quote = FALSE, row.names = FALSE)
    message(paste("Saved", nrow(all_outer), "outlier cells to CSV."))
    return(all_outer$Cell)
  } else {
    message("No outlier cells found.")
    return(character(0))
  }
}

res4 <- sort(unique(all_meta$cluster))
cells_to_remove <- outlier(clusters = res4, target_col = "cluster", label = "res4", meta = all_meta)

if (length(cells_to_remove) > 0) {
  merged_clean <- subset(merged, cells = cells_to_remove, invert = TRUE)
  message(paste("Removed", length(cells_to_remove), "outlier cells."))
} else {
  merged_clean <- merged
  message("No cells to remove.")
}

names(merged_clean@meta.data)

saveRDS(merged_clean, file = "chicken_hip_scRNA_harmony_annotation_8_samples_final.rds")

#visualization annotations
color_ct <- c("ExN1" = "#8B0000", "ExN2" = "#FFBABA", "ExN3" = "#F08080", "ExN4" = "#CD5C5C", "ExN5" = "#943C42", "InN1" = "#4682B4", "InN2" = "#4169E1", "InN3" = "#6C5CD0", "InN4" = "#0000CD", "InN5" = "#191970", "IPC1" = "#93AE0D", "IPC2" = "#78872F", "Pre" = "#DADA0B", "Astro" = "#A337B9", "Micro" = "#FF00FF", "OPC" = "#D2910D", "cOPC" = "#F9CF7B", "Oligo" = "#386C0A", "Epen" = "#02DCB8",  "Endo" = "#E1041A", "Tcell" = "#FD4B15", "VLMC1" = "#F4A460", "VLMC2" = "#D2691E", "VLMC3" = "#8B4513")
color_gp <- c("RJF" = "#A10C1B", "WLH" = "#5AB0A2")

p1 <- DimPlot(merged_clean, 
              reduction = "umap",
              group.by = "Subtype",
              cols = color_ct,
              label = TRUE,
              label.size = 5,
              pt.size = 0.8,
              repel = TRUE) +
  ggtitle("UMAP by Subtype") +
  theme(plot.title = element_text(hjust = 0.5))

table(merged_clean$Group)
p2 <- DimPlot(merged_clean, 
              reduction = "umap",
              group.by = "Group",
                        pt.size = 0.8,
              cols = color_gp) +
  ggtitle("UMAP by Group") +
  theme(plot.title = element_text(hjust = 0.5))

features <- c("ENSGALG00010020328", "ENSGALG00010023237", "ENSGALG00010005544", "ENSGALG00010000590","ENSGALG00010016496","ENSGALG00010013621")
p3 <- FeaturePlot_scCustom(seurat_object = merged_clean, features = features, order = F,num_columns = 3)

p4 <-FeaturePlot(merged, features = c("S.Score", "G2M.Score"), cols = c("#0066CC", "#FFFF00"), reduction="umap", label = FALSE, order = TRUE) + theme(plot.title = element_text(size=8))

pdf(file = "UMAP_plots_8_samples.pdf", width = 8, height = 6)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()
