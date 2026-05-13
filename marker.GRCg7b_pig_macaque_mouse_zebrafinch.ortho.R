suppressMessages({
    library(Seurat)
    library(ggplot2)
    library(plyr)
    library(getopt)
    })

options(future.globals.maxSize = 50000 * 1024^2) ## 50GB

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


all_combined <- readRDS(file = "../GRCg7b_pig_macaque_mouse_zebrafinch_check/integrate_pc30_nb10.GRCg7b_pig_macaque_mouse_zebrafinch.rds")

DefaultAssay(all_combined) <- "RNA"
dim(all_combined) #  74956 190339

kept_ortho <- intersect(as.character(ortho$V10), rownames(all_combined))
length(kept_ortho) # 10115

all_combined_ortho <- subset(all_combined, features = kept_ortho)
dim(all_combined_ortho) # 10115 190339

unique(all_combined_ortho$Celltype0)

all_combined_ortho <- SetIdent(all_combined_ortho, value = "Celltype0")

all_species <- unique(as.character(all_combined_ortho$Species))

for(i in 1:length(all_species)){
    seudt_tmp <- subset(all_combined_ortho, subset = Species==all_species[i])
    maker_tmp <- FindAllMarkers(object = seudt_tmp, min.pct = 0.25, only.pos = TRUE, verbose = FALSE)
    write.table(maker_tmp, paste0("allMarkers.celltype.integrated_ortho_pc30_nb10.", all_species[i], ".txt"), sep="\t", quote=F, col.names = TRUE, row.names = FALSE)
}

