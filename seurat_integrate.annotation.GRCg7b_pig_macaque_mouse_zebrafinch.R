suppressMessages({
    library(Seurat)
    library(ggplot2)
    library(qs)
    })
  
options(future.globals.maxSize = 50000 * 1024^2) ## 50GB

set.seed(123)

ortho <- read.table("../GRCg7b_pig_macaque_mouse_zebrafinch.1to1", sep="\t", header = FALSE, stringsAsFactors = FALSE, quote="")
head(ortho)
rownames(ortho) <- ortho[, 1]
dim(ortho) # 10164    14
ortho$V2 <- toupper(ortho$V2) # chicken
ortho$V4 <- toupper(ortho$V4) # pig
ortho$V7 <- toupper(ortho$V7) # macaque
ortho$V10 <- toupper(ortho$V10) # mouse
ortho$V13 <- toupper(ortho$V13) # zebrafinch

ortho <- ortho[ortho$V2!="" | ortho$V4!="" | ortho$V7!="" | ortho$V10!="" | ortho$V13!="", ]
dim(ortho) # 10164    14
ortho[ortho$V10=="", ]
ortho[ortho$V10=="PISD", ]
ortho[ortho$V10=="", ]$V10 <- ortho[ortho$V10=="", ]$V2
ortho[ortho$V4=="", ]$V4 <- ortho[ortho$V4=="", ]$V10
ortho[ortho$V7=="", ]$V7 <- ortho[ortho$V7=="", ]$V10
ortho[ortho$V13=="", ]$V13 <- ortho[ortho$V13=="", ]$V10


seudt0 <- readRDS(file = "../chicken_hip_scRNA_harmony_annotation_8_samples_final.rds")
dim(seudt0) ## 17254 40641

genes_chicken <- rownames(seudt0)
length(genes_chicken)
names(genes_chicken) <- genes_chicken

length(genes_chicken[!(genes_chicken %in% rownames(ortho))]) ## 8232

all(rownames(ortho) %in% genes_chicken) # FALSE

ortho_chicken <- ortho[rownames(ortho) %in% as.character(genes_chicken), ]
dim(ortho_chicken) ## 9022   14

dim(ortho_chicken[duplicated(ortho_chicken$V4), ]) # 0  14
dim(ortho_chicken[duplicated(ortho_chicken$V7), ]) # 2 14
dim(ortho_chicken[duplicated(ortho_chicken$V10), ]) # 0 14
dim(ortho_chicken[duplicated(ortho_chicken$V13), ]) # 0 14

genes_chicken[rownames(ortho_chicken)] <- ortho_chicken$V10
head(genes_chicken)

meta_chicken <- seudt0@meta.data
head(meta_chicken)

meta_chicken$Species <- "Chicken"

rownames(seudt0@assays$RNA@counts) = as.character(genes_chicken)

count_chicken <- GetAssayData(object = seudt0, slot = "counts")
count_chicken[1:10, 1:10]

which(rownames(count_chicken)=="")

seudt_chicken <- CreateSeuratObject(counts = count_chicken, project = "Chicken", meta.data = meta_chicken)
head(seudt_chicken@meta.data)

chicken_list <- SplitObject(seudt_chicken, split.by = "Sample")
chicken_list <- lapply(X = chicken_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})


data_mouse <- read.table("../ref/mouse_dev_Hochgerner/GSE104323_10X_expression_data_V2.tab.24185cells", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="", row.names = 1, check.names = FALSE)
dim(data_mouse) # 27933 24185
data_mouse <- data_mouse[rownames(data_mouse)!="PISD",]
rownames(data_mouse) <- toupper(rownames(data_mouse)) ## non-unique value when setting 'row.names': ??PISD??

length(intersect(rownames(data_mouse), ortho$V10)) # 9598

meta_mouse <- read.table("..//ref/mouse_dev_Hochgerner/GSE104323_metadata_barcodes_24185cells.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="")
dim(meta_mouse) #  24216    11
meta_mouse <- meta_mouse[meta_mouse$Sample.name..24185.single.cells. !="", ]
rownames(meta_mouse) <- meta_mouse$Sample.name..24185.single.cells.
meta_mouse <- meta_mouse[, -1]
meta_mouse$celltype <- meta_mouse$characteristics..cell.cluster
head(meta_mouse)
meta_mouse$Sample <- substr(rownames(meta_mouse), 1, 7)

meta_mouse <- meta_mouse[colnames(data_mouse), c("celltype", "characteristics..age", "characteristics..strain", "characteristics..sex.of.pooled.animals", "Sample")]
colnames(meta_mouse) <- c("Celltype", "Age", "Strain", "Sex_of_pooled_animals", "Sample")
meta_mouse$Species <- "Mouse"

seudt_mouse <- CreateSeuratObject(counts = data_mouse, project = "Mouse", meta.data = meta_mouse)

mouse_list <- SplitObject(seudt_mouse, split.by = "Sample")
mouse_list <- lapply(X = mouse_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})


data_macaque <- Read10X(data.dir = "../ref/human_pig_macaque_Franjic/GSE186538_Rhesus", gene.column = 1)
dim(data_macaque) # 33960 36107
rownames(data_macaque) <- toupper(rownames(data_macaque))
data_macaque[1:10, 1:10]

meta_macaque <- read.table("../ref/human_pig_macaque_Franjic/GSE186538_Rhesus/GSE186538_Rhesus_cell_meta.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="", row.names = 1)
dim(meta_macaque)
head(meta_macaque)
meta_macaque <- meta_macaque[, c("cluster", "samplename", "PMI")]
colnames(meta_macaque) <- c("Celltype", "Sample", "PMI")
meta_macaque$Species <- "Macaque"

seudt_macaque0 <- CreateSeuratObject(counts = data_macaque, meta.data = meta_macaque)

genes_macaque <- rownames(seudt_macaque0)
length(genes_macaque) # 33960
names(genes_macaque) <- genes_macaque

length(genes_macaque[!(genes_macaque %in% ortho$V7)]) ## 24173

all(ortho$V7 %in% genes_macaque) # FALSE

ortho_macaque <- ortho[ortho$V7 %in% as.character(genes_macaque), ]
dim(ortho_macaque) ## 9788   14

dim(ortho_macaque[duplicated(ortho_macaque$V4), ]) #  0 14
dim(ortho_macaque[duplicated(ortho_macaque$V7), ]) # 1 14
dim(ortho_macaque[duplicated(ortho_macaque$V10), ]) # 0 14
ortho_macaque[duplicated(ortho_macaque$V7), ]
ortho_macaque[ortho_macaque$V7 %in% c("SKP1"),]
ortho_macaque <- ortho_macaque[!(ortho_macaque$V6 %in% c("ENSMMUG00000013047")), ]
rownames(ortho_macaque) <- ortho_macaque$V7

genes_macaque[rownames(ortho_macaque)] <- ortho_macaque$V10
dup_macaque <- unique(as.character(genes_macaque[duplicated(genes_macaque)]))

genes_macaque[genes_macaque %in% dup_macaque] <- names(genes_macaque[genes_macaque %in% dup_macaque])

head(genes_macaque)
length(genes_macaque) # 33960

rownames(seudt_macaque0@assays$RNA@counts) = as.character(genes_macaque)

count_macaque <- GetAssayData(object = seudt_macaque0, slot = "counts")
count_macaque[1:10, 1:10]

seudt_macaque <- CreateSeuratObject(counts = count_macaque, project = "Macaque", meta.data = meta_macaque)

macaque_list <- SplitObject(seudt_macaque, split.by = "Sample")
macaque_list <- lapply(X = macaque_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})


data_pig <- Read10X(data.dir = "ref/human_pig_macaque_Franjic/GSE186538_Pig", gene.column = 1)
dim(data_pig) # 29886 36851
rownames(data_pig) <- toupper(rownames(data_pig))
data_pig[1:10, 1:10]

meta_pig <- read.table("..//ref/human_pig_macaque_Franjic/GSE186538_Pig/GSE186538_Pig_cell_meta.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="", row.names = 1)
dim(meta_pig)
head(meta_pig)

meta_pig <- meta_pig[, c("cluster", "samplename", "PMI")]
colnames(meta_pig) <- c("Celltype", "Sample", "PMI")
meta_pig$Species <- "Pig"

seudt_pig0 <- CreateSeuratObject(counts = data_pig, meta.data = meta_pig)

genes_pig <- rownames(seudt_pig0)
length(genes_pig) # 29886
names(genes_pig) <- genes_pig

length(genes_pig[!(genes_pig %in% ortho$V4)]) ## 20332

all(ortho$V4 %in% genes_pig) # FALSE

ortho_pig <- ortho[ortho$V4 %in% as.character(genes_pig), ]
dim(ortho_pig) ## 9554   14

dim(ortho_pig[duplicated(ortho_pig$V4), ]) #  0 14
dim(ortho_pig[duplicated(ortho_pig$V7), ]) # 1 14
dim(ortho_pig[duplicated(ortho_pig$V10), ]) # 0 14

rownames(ortho_pig) <- ortho_pig$V4

genes_pig[rownames(ortho_pig)] <- ortho_pig$V10
dup_pig <- unique(as.character(genes_pig[duplicated(genes_pig)]))

genes_pig[genes_pig %in% dup_pig] <- names(genes_pig[genes_pig %in% dup_pig])

head(genes_pig)
length(genes_pig) # 29886

sum(genes_pig %in% c("MKI67", "CENPF", "TOP2A")) # 3

rownames(seudt_pig0@assays$RNA@counts) = as.character(genes_pig)

count_pig <- GetAssayData(object = seudt_pig0, slot = "counts")
count_pig[1:10, 1:10]

seudt_pig <- CreateSeuratObject(counts = count_pig, project = "Pig", meta.data = meta_pig)

pig_list <- SplitObject(seudt_pig, split.by = "Sample")
pig_list <- lapply(X = pig_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})


data_finch <- qread("../ref/zebra_finch_Bradley/HVC_RA.qs")
dim(data_finch) # 3000 29497
head(data_finch@meta.data)
meta_finch <- data_finch@meta.data
head(meta_finch)
table(meta_finch$species, meta_finch$position2, meta_finch$orig.ident)
meta_finch$Sample <- meta_finch$orig.ident
meta_finch[meta_finch$Sample=="Channel1", ]$Sample <- paste(meta_finch[meta_finch$Sample=="Channel1", ]$position2, meta_finch[meta_finch$Sample=="Channel1", ]$Sample, sep="_")
meta_finch[meta_finch$Sample=="HVC_AV_C_seurat", ]$Sample <- "zf_HVC1"
meta_finch[meta_finch$Sample=="HVC_CTRL_C_seurat", ]$Sample <- "zf_HVC2"
meta_finch[meta_finch$Sample=="HVC_X_C_seurat", ]$Sample <- "zf_HVC3"
meta_finch[meta_finch$Sample=="SeuratProject", ]$Sample <- "zf_ra"
meta_finch$Species <- "Bengalese_finch"
meta_finch[meta_finch$species=="zf", ]$Species <- "Zebra_finch"
meta_finch <- meta_finch[, c("cluster_int_sub2", "position2", "Sample", "Species")]
colnames(meta_finch) <- c("Celltype", "Region", "Sample", "Species")
head(meta_finch)

DefaultAssay(data_finch) <- "RNA"
dim(data_finch) # 19857 29497

genes_finch <- toupper(rownames(data_finch))
length(genes_finch) # 19857
names(genes_finch) <- genes_finch

length(genes_finch[!(genes_finch %in% ortho$V13)]) ## 10114

all(ortho$V13 %in% genes_finch) # FALSE

ortho_finch <- ortho[ortho$V13 %in% as.character(genes_finch), ]
dim(ortho_finch) ## 9743   14

dim(ortho_finch[duplicated(ortho_finch$V4), ]) #  0 14
dim(ortho_finch[duplicated(ortho_finch$V7), ]) # 1 14
dim(ortho_finch[duplicated(ortho_finch$V10), ]) # 0 14
dim(ortho_finch[duplicated(ortho_finch$V13), ]) # 0 14

rownames(ortho_finch) <- ortho_finch$V13

genes_finch[rownames(ortho_finch)] <- ortho_finch$V10
dup_finch <- unique(as.character(genes_finch[duplicated(genes_finch)]))

genes_finch[genes_finch %in% dup_finch] <- names(genes_finch[genes_finch %in% dup_finch])

head(genes_finch)
length(genes_finch) # 19857

rownames(data_finch$RNA@counts) = as.character(genes_finch)

count_finch <- GetAssayData(object = data_finch, slot = "counts")
count_finch[1:10, 1:10]

seudt_finch <- CreateSeuratObject(counts = count_finch, project = "Finch", meta.data = meta_finch)

finch_list <- SplitObject(seudt_finch, split.by = "Sample")
finch_list <- lapply(X = finch_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})


all_list <- c(chicken_list, mouse_list, macaque_list, pig_list, finch_list)

all_features <- SelectIntegrationFeatures(object.list = all_list)
length(all_features)

write.table(all_features, "../integrateFeatures.GRCg7b_pig_macaque_mouse_zebrafinch.txt", sep="\t", quote=F, col.names = FALSE, row.names = FALSE)

all_list <- lapply(X = all_list, FUN = function(x) {
  x <- ScaleData(x, features = all_features, verbose = FALSE)
  x <- RunPCA(x, features = all_features, verbose = FALSE)
})

all_anchors <- FindIntegrationAnchors(object.list = all_list, anchor.features = all_features, reduction = "rpca", verbose = FALSE)

all_combined <- IntegrateData(anchorset = all_anchors, verbose = FALSE)

DefaultAssay(all_combined) <- "integrated"

all_combined <- ScaleData(all_combined, verbose = FALSE)
all_combined <- RunPCA(all_combined, npcs = 30, verbose = FALSE)
all_combined <- RunUMAP(all_combined, reduction = "pca", dims = 1:30)
all_combined <- FindNeighbors(all_combined, reduction = "pca", dims = 1:30)
all_combined <- FindClusters(all_combined, resolution = 0.1, verbose = FALSE)
all_combined <- FindClusters(all_combined, resolution = 0.2, verbose = FALSE)
all_combined <- FindClusters(all_combined, resolution = 0.3, verbose = FALSE)
all_combined <- FindClusters(all_combined, resolution = 0.4, verbose = FALSE)
all_combined <- FindClusters(all_combined, resolution = 0.5, verbose = FALSE)
all_combined <- FindClusters(all_combined, resolution = 0.6, verbose = FALSE)
all_combined <- FindClusters(all_combined, resolution = 0.7, verbose = FALSE)
all_combined <- FindClusters(all_combined, resolution = 1, verbose = FALSE)
all_combined <- FindClusters(all_combined, resolution = 1.5, verbose = FALSE)

all_meta <- cbind(all_combined@meta.data, all_combined@reductions$umap@cell.embeddings)
head(all_meta)

write.table(all_meta, "../integrate_meta_data.GRCg7b_pig_macaque_mouse_zebrafinch.txt", sep="\t", quote=F, col.names = TRUE, row.names = TRUE)

saveRDS(all_combined, file = "../integrate.GRCg7b_pig_macaque_mouse_zebrafinch.rds")

dim(all_combined) # 2000 167281

umap.cl01 <- DimPlot(all_combined, reduction = "umap", label = TRUE, repel = FALSE, group.by = "integrated_snn_res.0.1") + NoLegend()
umap.cl02 <- DimPlot(all_combined, reduction = "umap", label = TRUE, repel = FALSE, group.by = "integrated_snn_res.0.2") + NoLegend()
umap.cl03 <- DimPlot(all_combined, reduction = "umap", label = TRUE, repel = FALSE, group.by = "integrated_snn_res.0.3") + NoLegend()
umap.cl04 <- DimPlot(all_combined, reduction = "umap", label = TRUE, repel = FALSE, group.by = "integrated_snn_res.0.4") + NoLegend()
umap.cl05 <- DimPlot(all_combined, reduction = "umap", label = TRUE, repel = FALSE, group.by = "integrated_snn_res.0.5") + NoLegend()
umap.cl06 <- DimPlot(all_combined, reduction = "umap", label = TRUE, repel = FALSE, group.by = "integrated_snn_res.0.6") + NoLegend()
umap.cl07 <- DimPlot(all_combined, reduction = "umap", label = TRUE, repel = FALSE, group.by = "integrated_snn_res.0.7") + NoLegend()
umap.cl1 <- DimPlot(all_combined, reduction = "umap", label = TRUE, repel = FALSE, group.by = "integrated_snn_res.1") + NoLegend()
umap.cl1.5 <- DimPlot(all_combined, reduction = "umap", label = TRUE, repel = FALSE, group.by = "integrated_snn_res.1.5") + NoLegend()

pdf("UMAP.integrate.res01_03.pdf", w=15, h=6)
umap.cl01 + umap.cl02 + umap.cl03
dev.off()

pdf("UMAP.integrate.res04_06.pdf", w=15, h=6)
umap.cl04 + umap.cl05 + umap.cl06
dev.off()

pdf("UMAP.integrate.res07_1.5.pdf", w=15, h=6)
umap.cl07 + umap.cl1 + umap.cl1.5
dev.off()

all_combined <- readRDS(file = "integrate.GRCg7b_pig_macaque_mouse_zebrafinch.rds")

DefaultAssay(all_combined) <- "integrated"

all_meta <- cbind(all_combined@meta.data, all_combined@reductions$umap@cell.embeddings)
unique(all_meta[, c("Celltype", "Species")])
write.table(unique(all_meta[, c("Celltype", "Species")]), "celltype.integrate.GRCg7b_pig_macaque_mouse_zebrafinch.txt", sep="\t", quote=F, col.names = TRUE, row.names = FALSE)

all_meta$Celltype0 <- as.character(all_meta$Celltype)
all_meta$Celltype0[all_meta$Species == "Chicken"] <- as.character(all_meta$Subtype[all_meta$Species == "Chicken"])

#==========Excitatory neurons ==========
all_meta[all_meta$Celltype0 %in% c("GC-juv","GC-adult"), ]$Celltype0 <- "GC"
all_meta[all_meta$Celltype0 %in% c("RA_Glut-1","RA_Glut-2","RA_Glut-3"), ]$Celltype0 <- "RA-ExN"
all_meta[all_meta$Celltype0 %in% c("HVC_Glut-1","HVC_Glut-2","HVC_Glut-3","HVC_Glut-4","HVC_Glut-5"), ]$Celltype0 <- "HVC-ExN"
all_meta[all_meta$Celltype0=="MC", ]$Celltype0 <- "Mossy"
all_meta[all_meta$Celltype0 == "CR", ]$Celltype0 <- "CR"

#==========Inhibitory neurons==========
all_meta[all_meta$Celltype0 %in% c(
  "GABA","GABA-1-1","GABA-1-2","GABA-2","GABA-3","GABA-4",
  "GABA-5-1","GABA-5-2","GABA-5-3","GABA-6","GABA-7","GABA-8"
), ]$Celltype0 <- "InN"

all_meta[all_meta$Celltype0 %in% c("InN1","InN2","InN3","InN4","InN5"), ]$Celltype0 <- "InN" 
all_meta[all_meta$Celltype0 == "Immature-GABA", ]$Celltype0 <- "InN"

#==========Neural progenitor cells==========
all_meta[all_meta$Celltype0 == "nIPC-perin", ]$Celltype0 <- "nIPC"
all_meta[all_meta$Celltype0 %in% c("Pre-1","Pre-2","Pre-3","Pre-4"), ]$Celltype0 <- "Pre"
all_meta[all_meta$Celltype0 == "GABA-Pre", ]$Celltype0 <- "Pre"
all_meta[all_meta$Celltype0 == "NB", ]$Celltype0 <- "Pre"
all_meta[all_meta$Celltype0 == "RGL_young", ]$Celltype0 <- "RGL"

#==========Glial cells==========
all_meta[all_meta$Celltype0 %in% c("Astro-adult","Astro-juv"), ]$Celltype0 <- "Astro"
all_meta[all_meta$Celltype0 == "MOL", ]$Celltype0 <- "Oligo"
all_meta[all_meta$Celltype0 == "NFOL", ]$Celltype0 <- "Newly-Oligo"
all_meta[all_meta$Celltype0 == "MiCajal-Retziusoglia", ]$Celltype0 <- "Micro"
all_meta[all_meta$Celltype0 == "Ependymal", ]$Celltype0 <- "Epen"

#==========Vascular and RBC==========
all_meta[all_meta$Celltype0 == "Endothelial", ]$Celltype0 <- "Vascular"
all_meta[all_meta$Celltype0 == "Endo", ]$Celltype0 <- "Vascular"
all_meta[all_meta$Celltype0 %in% c("VLMC1","VLMC2","VLMC3","VLMC"), ]$Celltype0 <- "Vascular"
all_meta[all_meta$Celltype0 == "Vas", ]$Celltype0 <- "Vascular"
all_meta[all_meta$Celltype0 == "Mural", ]$Celltype0 <- "Vascular"
all_meta[all_meta$Celltype0 == "RBC", ]$Celltype0 <- "RBC"

#==========Immune==========
all_meta[all_meta$Celltype0 == "immune", ]$Celltype0 <- "Immune"
all_meta[all_meta$Celltype0 == "PVM", ]$Celltype0 <- "Immune"
unique(all_meta$Celltype0)

color_ct <- c("ExN1" = "#8B0000", "ExN2" = "#FFBABA", "ExN3" = "#F08080", "ExN4" = "#CD5C5C", "ExN5" = "#943C42",
              "HVC-ExN" = "#DC176D", "RA-ExN" = "#842A6C", "InN" = "#3B505F", 'Immature-Pyr'= '#C6C696',
              "IPC1" = "#93AE0D", "IPC2" = "#78872F",  'nIPC'= '#DBDB82', "Pre" = "#DADA0B", 
              'GC'= '#386cb0','Immature-GC'= '#80b1d3','CA1 Sub'= '#FFBAD9','CA2-3'= '#F08BD5','CA3-Pyr'= '#FA4E99','Mossy'= '#decbe4',
              "Astro" = "#A337B9", "Micro" = "#FF00FF", 
              "OPC" = "#D2910D", "cOPC" = "#F9CF7B", "Oligo" = "#386C0A", "Epen" = "#02DCB8",  
              "Tcell" = "#FD4B15", 'RGL'= '#575102','Cajal-Retzius'= '#68D9EB',
              'Immune'= '#BF4D24','RBC'= '#F9D3B1','Vascular'= '#F4A460'
)

all_combined <- AddMetaData(all_combined, all_meta[colnames(all_combined), ]$Celltype0, col.name = "Celltype0")

all_combined@meta.data$Celltype0 <- factor(all_combined@meta.data$Celltype0, levels = names(color_ct), order = TRUE)

all_combined <- SetIdent(all_combined, value = "Celltype0")

pdf("UMAP.integrate.species.pdf", w=5, h=4)
DimPlot(all_combined, reduction = "umap", group.by = "Species")
dev.off()

pdf("UMAP.integrate.celltype.pdf", w=5, h=4)
DimPlot(all_combined, reduction = "umap", group.by = "Celltype0", cols = color_ct, label = FALSE, label.size = 2, order = TRUE) + theme(legend.text = element_text(size = rel(0.3)), legend.key.size = unit(0.7, "mm"), legend.spacing.y = unit(0, "mm")) + guides(color = guide_legend(ncol = 1, override.aes = list(size = 1)))
dev.off()

pdf("UMAP.integrate.species_split.pdf", w=21, h=6)
DimPlot(all_combined, reduction = "umap", group.by = "Celltype0", split.by = "Species", cols = color_ct, label = TRUE, repel = FALSE, label.size = 2) + NoLegend()
dev.off()

DefaultAssay(all_combined) <- "RNA"
dim(all_combined) # 74956 190339

all_combined <- SetIdent(all_combined, value = "Celltype0")

all_species <- unique(as.character(all_meta$Species))

for(i in 1:length(all_species)){
  seudt_tmp <- subset(all_combined, subset = Species==all_species[i])
  maker_tmp <- FindAllMarkers(object = seudt_tmp, min.pct = 0.25, only.pos = TRUE, verbose = FALSE)
  write.table(maker_tmp, paste0("allMarkers.celltype.integrated.", all_species[i], ".txt"), sep="\t", quote=F, col.names = TRUE, row.names = FALSE)
}

kept_ortho <- intersect(as.character(ortho$V10), rownames(all_combined))
length(kept_ortho) # 10115

all_combined_ortho <- subset(all_combined, features = kept_ortho)
dim(all_combined_ortho) # 10115 190339

all_species <- unique(as.character(all_combined_ortho$Species))

for(i in 1:length(all_species)){
  seudt_tmp <- subset(all_combined_ortho, subset = Species==all_species[i])
  maker_tmp <- FindAllMarkers(object = seudt_tmp, min.pct = 0.25, only.pos = TRUE, verbose = FALSE)
  write.table(maker_tmp, paste0("allMarkers.celltype.integrated_ortho.", all_species[i], ".txt"), sep="\t", quote=F, col.names = TRUE, row.names = FALSE)
}

all_combined <- readRDS(file = "integrate.GRCg7b_pig_macaque_mouse_zebrafinch.rds")
