##conda activate seurat_v5
setwd("/projects/chicken_hip/expression_matrix_new2/pyscenic")
rm(list=ls())
gc()

library(dplyr)
library(data.table)
library(ggplot2) 
library(SCopeLoomR)
library(SCENIC)
library(BiocParallel)
library(ComplexHeatmap)
library(pheatmap)
library(AUCell)

#-----------------------------------------------------------regulons
loom <- open_loom('RJF_SCENIC.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 

regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])
embeddings <- get_embeddings(loom)  
close_loom(loom)
rownames(regulonAUC)
head(regulonAUC)
names(regulons)
head(regulons)

regulon_df <- data.frame(
  TF = rep(names(regulons), lengths(regulons)),
  Target = unlist(regulons)
)

regulon_df$TF <- toupper(regulon_df$TF)
regulon_df$Target <- toupper(regulon_df$Target)
write.csv(regulon_df, "RJF_regulons_TF_targets.csv", row.names = FALSE)

#-----------------------------------------------------------matrix
# Read metadata
metadata <- read.csv("RJF_metadata.csv", 
                     header = TRUE,
                     stringsAsFactors = FALSE,
                     na.strings = c("", "NA", "NULL"))

head(metadata)

# Check for missing values in the X column
if(any(is.na(metadata$X))) {
  cat("Found", sum(is.na(metadata$X)), "missing values in column X\n")
  cat("Removing rows with missing X values...\n")
  metadata <- metadata[!is.na(metadata$X), ]
}

# Now set row names safely
rownames(metadata) <- metadata$X

# Verify no NAs in row names
if(any(is.na(rownames(metadata)))) {
  stop("Still have NA values in row names after cleaning!")
}

# Subset regulonAUC to match metadata
sub_regulonAUC <- regulonAUC[, match(rownames(metadata), colnames(regulonAUC))]
rownames(sub_regulonAUC) <- toupper(rownames(sub_regulonAUC))
head(sub_regulonAUC)

# Check if dimensions match
identical(colnames(sub_regulonAUC), rownames(metadata))

# Create cellTypes data frame - FIXED: use nrow instead of colnames
cellTypes <- data.frame(row.names = rownames(metadata),  # Changed from colnames(metadata)
                        celltype = metadata$Celltype)
head(cellTypes)

# Define selected resolution
selectedResolution <- "celltype" 

# Split cells by celltype - FIXED: use rownames(cellTypes)
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[, selectedResolution])

# Explore thresholds and assign cells
cells_assignment <- AUCell_exploreThresholds(sub_regulonAUC, 
                                             plotHist = FALSE, 
                                             assignCells = TRUE)

thresholds <- getThresholdSelected(cells_assignment)
cells_assignment_binary <- getAssignments(cells_assignment)

# Calculate regulon activity by celltype
regulon_activity_by_celltype <- lapply(cellsPerGroup, function(cell_group) {
  if (length(cell_group) > 0) {
    sapply(cells_assignment_binary, function(regulon_cells) {
      mean(cell_group %in% regulon_cells)
    })
  } else {
    rep(NA, length(cells_assignment_binary))
  }
})

# Convert to data frame
activity_df <- as.data.frame(regulon_activity_by_celltype)
head(activity_df)

# Calculate max activities (excluding last column if it's empty)
if(ncol(activity_df) > 1) {
  max_activities <- apply(activity_df[, 1:ncol(activity_df)], 1, max, na.rm = TRUE)
} else {
  max_activities <- activity_df[, 1]
}

summary(max_activities)
quantile(max_activities, probs = c(0.5, 0.75, 0.9, 0.95), na.rm = TRUE)

# Save results
write.csv(activity_df, "RJF_regulon_pro_by_celltype.csv", row.names = TRUE)

# Optional: Print summary information
cat("\nAnalysis completed successfully!\n")
cat("Number of cells:", nrow(metadata), "\n")
cat("Number of regulons:", nrow(activity_df), "\n")
cat("Number of cell types:", ncol(activity_df), "\n")
cat("Cell types:", paste(colnames(activity_df), collapse = ", "), "\n")