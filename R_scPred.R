####
library(Seurat)
library(scPred)
library(tidyverse) 
library(patchwork)

reference <- readRDS('./Macaque_cells.rds')

reference <- reference %>%
  getFeatureSpace("celltype1") %>%
  trainModel(allowParallel = TRUE) %>%
  trainModel(model = "mda", reclassify = "IPC")

saveRDS(reference, "reference-Macaque_cells.rds")

seurat_obj <- readRDS('./final_integrate.GRCg7b_pig_macaque_mouse_zebrafinch.rds')

seurat_obj <- scPredict(seurat_obj, reference)

seurat_obj <- RunUMAP(seurat_obj, reduction = "scpred", dims = 1:30)

print(table(seurat_obj$scpred_prediction))

saveRDS(seurat_obj, 'integrate.GRCg7b_pig_macaque_mouse_zebrafinch_scPred-refMacaque.rds')


####
library(openxlsx)

getwd()
rm(list=ls())
gc()

seurat_obj <- readRDS('integrate.GRCg7b_pig_macaque_mouse_zebrafinch_scPred-refMacaque.rds')

high_prob_cells <- seurat_obj@meta.data[seurat_obj@meta.data$scpred_IPC > 0.9, ]

result_df <- data.frame(
  Cell_Barcode = rownames(high_prob_cells),
  scpred_IPC_Probability = high_prob_cells$scpred_IPC
)


####
## Next, repeat the above process using Pig and Mouse as reference datasets. Finally, the intersection of IPC cell results predicted based on macaque, pig, and mouse reference lines was taken to obtain a high confidence cell population with cross species consistency.(venn-nIPC_cells.xlsx)



####
library(Seurat)
library(openxlsx)
library(tidydr)  
library(ggplot2)

setwd("E:/KIZ/12.RJF_project/Scpred")
getwd()
rm(list=ls())
gc()

seurat_obj <- readRDS('./final_integrate.GRCg7b_pig_macaque_mouse_zebrafinch.rds')

nIPC_cells_df <- read.xlsx("./venn-nIPC_cells.xlsx", sheet = 1)
nIPC_cell_names <- nIPC_cells_df[, 1]  

nIPC_cell_names <- as.character(nIPC_cell_names)


existing_cells <- nIPC_cell_names[nIPC_cell_names %in% colnames(seurat_obj)]
non_existing_cells <- nIPC_cell_names[!nIPC_cell_names %in% colnames(seurat_obj)]


seurat_obj$IPC_label <- "Other"  
seurat_obj$IPC_label[colnames(seurat_obj) %in% existing_cells] <- "nIPC"


label_counts <- table(seurat_obj$IPC_label)

print(label_counts)

##
celltype_cols <- c(
  "nIPC" = "#dc8139",
  "Other" = "grey"
  
)

species_order <- c("Macaque", "Mouse", "Pig", "Chicken", "Bengalese_finch", "Zebra_finch")
seurat_obj$Species <- factor(seurat_obj$Species, levels = species_order)


DimPlot(seurat_obj, reduction = "umap", group.by = "IPC_label", label = TRUE,
        cols = celltype_cols, split.by = "Species", ncol = 3, raster=FALSE) +
  theme_dr(xlength = 0.22, ylength = 0.22, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")) +
  theme(panel.grid = element_blank())

saveRDS(seurat_obj, "final_integrate.GRCg7b_pig_macaque_mouse_zebrafinch_scpred.rds")


####
setwd("E:/KIZ/12.RJF_project/Scpred")
getwd()
rm(list=ls())
gc()

seurat_obj <- readRDS('./final_integrate.GRCg7b_pig_macaque_mouse_zebrafinch_scpred.rds')

genes_to_check <- c( 'MKI67','TOP2A','CSF1R','GZMA','CD247','DCX','MEX3A','MEX3B',
                     'PAX6','SOX4','SOX11')

nIPC_cells <- subset(seurat_obj, 
                     subset = (Species == "Chicken" & Celltype %in% c("IPC1", "IPC2")) |
                       (Species == "Mouse" & Celltype %in% c("nIPC", "nIPC-perin")) |
                       (Species %in% c("Macaque", "Pig") & Celltype == "nIPC") |
                       (Species %in% c("Bengalese_finch", "Zebra_finch") & IPC_label == "nIPC"))


nIPC_cells$display_label <- ifelse(nIPC_cells$Species %in% c("Chicken", "Mouse"),
                                   paste0(nIPC_cells$Species, "_", nIPC_cells$Celltype),
                                   paste0(nIPC_cells$Species, "_nIPC"))

DotPlot(nIPC_cells, 
        features = genes_to_check,
        assay = 'RNA', 
        scale = TRUE,
        group.by = "display_label") +  
  coord_flip() + 
  scale_color_gradientn(colors = c("#e9e82c", "#c5407f")) +
  theme(
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1)  
  )


