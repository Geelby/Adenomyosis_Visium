#!/bin/bash


#SBATCH --job-name=pyscenicstep1        
#SBATCH --partition=64c512g       
#SBATCH -N 1                
#SBATCH --ntasks-per-node=64   
#SBATCH --mail-type=end           
#SBATCH --mail-user=//////
#SBATCH --output=%jpyscenicstep1.out
#SBATCH --error=%jpyscenicstep1.err
#SBATCH --exclusive

module load miniconda3 
source activate pyscenic 
pyscenic grn \
./AM_filtered_scenic.loom \
./data/allTFs_hg38.txt \
-o adj.csv --num_workers 64













