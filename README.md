# chicken_hip_snRNA_code
## Single-Cell Transcriptomic Perspectives on Hippocampal Evolution

This repository contains the code used for the primary analyses described in the manuscript, “From Wild to Domestic: Single-Cell Transcriptomic Perspectives on Hippocampal Plasticity and Evolution.”

The script `qc_process.R` performs quality control and filtering of nuclei for individual samples, whereas `harmony_umap_anno.R` is used for sample integration, dimensionality reduction, clustering, and cell-type annotation across all samples.

Differential abundance analysis was performed using `miloR.R`, while differential gene expression analysis between WLH and RJF was conducted using the custom script `wlh_rjf_deg.R`.

Pseudotime and RNA velocity analyses were performed using `monocle3.R` and `velocyto_python.R`, respectively.

Machine learning–based cell-type prediction was carried out using `scPred.R` and `XGBoost.py`.

SCENIC analyses were conducted using `pySCENIC_RJF.sh` and `pySCENIC_WLH.sh`, which invoke a series of preprocessing and visualization scripts, including `WLH_counts_metadata.R`, `RJF_counts_metadata.R`, `csv_to_loom_RJF.py`, `csv_to_loom_WLH.py`, `visualization_WLH.R`, and `visualization_RJF.R`.

Cross-species datasets from chicken, macaque, mouse, pig, and songbirds were integrated using the Seurat pipeline implemented in `seurat_integrate.annotation.GRCg7b_pig_macaque_mouse_zebrafinch.R`, while orthologous marker gene analyses were performed using `marker.GRCg7b_pig_macaque_mouse_zebrafinch.ortho.R`.
