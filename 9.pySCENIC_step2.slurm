#!/bin/bash


#SBATCH --job-name=pyscenicstep2       
#SBATCH --partition=64c512g       
#SBATCH -N 1               
#SBATCH --ntasks-per-node=64   
#SBATCH --mail-type=end          
#SBATCH --mail-user=geelby@sjtu.edu.cn
#SBATCH --output=%jpyscenicstep2.out
#SBATCH --error=%jpyscenicstep2.err
#SBATCH --exclusive

module load miniconda3  
source activate pyscenic 
pyscenic ctx ./adj.csv \
    ./data/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather ./data/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
    --annotations_fname ./data/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname ./AM_filtered_scenic.loom \
    --output reg.csv \
    --mask_dropouts \
    --num_workers 64
