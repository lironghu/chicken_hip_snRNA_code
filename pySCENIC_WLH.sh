#!/bin/bash

SEURAT_ENV1="/projects/software/anaconda3/envs/seurat_v5"
SCENIC_ENV2="/projects/software/anaconda3/envs/pyscenic"
DIR="/projects/chicken_hip/expression_matrix_new2/pyscenic"
REF="/projects/scenic/cisTarget_databases"

source /projects/software/anaconda3/etc/profile.d/conda.sh
conda activate "$SEURAT_ENV1"
Rscript WLH_counts_metadata.R
conda deactivate

source /projects/software/anaconda3/etc/profile.d/conda.sh
conda activate "$SCENIC_ENV2"
python csv_to_loom_WLH.py
conda deactivate

source /projects/software/anaconda3/etc/profile.d/conda.sh
conda activate "$SCENIC_ENV2"
input_loom="$DIR/WLH_counts.loom"
tfs="$REF/mm_mgi_tfs.txt"
feather="$REF/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
tbl="$REF/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
ls $tfs $feather $tbl $input_loom

#grn
pyscenic grn \
--num_workers 10 \
--output adj.WLH.tsv \
--method grnboost2 \
$input_loom \
$tfs

#cistarget
pyscenic ctx \
adj.WLH.tsv \
$feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 10 \
--mask_dropouts

#AUCell
pyscenic aucell \
$input_loom \
reg.csv \
--output WLH_SCENIC.loom \
--num_workers 5

conda deactivate

echo "====DONE===="
echo "SCENIC_results: WLH_SCENIC.loom"
echo "output_file:"
ls -lh WLH_SCENIC.loom

source /projects/software/anaconda3/etc/profile.d/conda.sh
conda activate "$SEURAT_ENV1"
Rscript visualization_WLH.R
conda deactivate

echo "====DONE===="
echo "Visualization_WLH"

