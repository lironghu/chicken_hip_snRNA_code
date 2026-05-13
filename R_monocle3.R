library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(tidydr)

#
setwd("E:/KIZ/12.RJF_project/monocle3")
rm(list=ls())
gc()

seurat_obj <- readRDS('../chicken_hip_scRNA_harmony_annotation_8samples_final.rds')
seurat_obj <- subset(seurat_obj, subset = Subtype %in% c("ExN1", "ExN2", "ExN3", "ExN4", "ExN5", "Pre"))
seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay5")

# Construct a Monocle3 object
umap_embeddings <- Embeddings(seurat_obj, reduction = "umap")
cell_metadata <- seurat_obj@meta.data
seurat_obj_new <- JoinLayers(seurat_obj[["RNA"]])
data <- as.matrix(seurat_obj_new$data)
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
reducedDims(cds)[["UMAP"]] <- umap_embeddings

# Trajectory analysis
cds <- cluster_cells(cds, reduction_method = "UMAP", k = 20, partition_qval = 0.05)
cds <- learn_graph(cds, use_partition = FALSE, close_loop = TRUE)

get_earliest_principal_node <- function(cds, cell_type = "Pre", column = "Subtype"){
  cell_ids <- which(colData(cds)[, column] == cell_type)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_node <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))]
  return(root_node)
}

root_pr_node <- get_earliest_principal_node(cds)
cds <- order_cells(cds, root_pr_nodes = root_pr_node)

# 
custom_colors <- c("ExN1" = "#8B0000", "ExN2" = "#FFBABA", "ExN3" = "#F08080", 
                   "ExN4" = "#CD5C5C", "ExN5" = "#943C42", "Pre" = "#DADA0B")

p1 <- plot_cells(cds, color_cells_by = "Subtype", label_cell_groups = F, label_leaves = F, 
                 label_branch_points = F, graph_label_size = 1.5, cell_size = 1.5) +
  scale_color_manual(values = custom_colors) + theme(legend.position = "right")

p2 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_leaves = F, 
                 label_branch_points = F, graph_label_size = 1.5, cell_size = 1.5)

p1 + p2

# 
cds$pseudotime_bin <- cut(pseudotime(cds), breaks = 5, include.lowest = TRUE)

group_proportions <- colData(cds) %>%
  as.data.frame() %>%
  group_by(pseudotime_bin, Group) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(pseudotime_bin) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

p_bar <- ggplot(group_proportions, aes(x = pseudotime_bin, y = Proportion, fill = Group)) +
  geom_col(position = "stack") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_bubble <- ggplot(group_proportions, aes(x = pseudotime_bin, y = Group)) +
  geom_point(aes(size = Count, color = Proportion), alpha = 0.8) +
  scale_size_continuous(range = c(3, 12)) +
  scale_color_viridis_c() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_bar + p_bubble