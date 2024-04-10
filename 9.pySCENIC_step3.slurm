#!/bin/bash


#SBATCH --job-name=pyscenicstep3       
#SBATCH --partition=64c512g       
#SBATCH -N 1               
#SBATCH --ntasks-per-node=64   
#SBATCH --mail-type=end          
#SBATCH --mail-user=geelby@sjtu.edu.cn
#SBATCH --output=%jpyscenicstep3.out
#SBATCH --error=%jpyscenicstep3.err
#SBATCH --exclusive

module load miniconda3  
source activate pyscenic 
pyscenic aucell \
    ./AM_filtered_scenic.loom \
    ./reg.csv \
    --output AM_out_SCENIC.loom \
    --num_workers 64
