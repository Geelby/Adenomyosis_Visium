#!/bin/bash

#SBATCH --job-name=velocyto        
#SBATCH --partition=cpu        
#SBATCH -N 2               
#SBATCH --ntasks-per-node=40   
#SBATCH --mail-type=end           
#SBATCH --mail-user=//////
#SBATCH --output=%jvelocyto.out
#SBATCH --error=%jvelocyto.err
#SBATCH --exclusive

module load miniconda3  
source activate velocyto 

rmsk_gtf=//////datavelocyto/hg38_repeat_rmsk.gtf
cellranger_outDir=//////spaceranger_tutorial/UM/
cellranger_gtf=//////datavelocyto/refdata-gex-GRCh38-2020-A/genes/genes.gtf

velocyto run10x -m $rmsk_gtf  $cellranger_outDir $cellranger_gtf

rmsk_gtf=//////datavelocyto/hg38_repeat_rmsk.gtf
cellranger_outDir=//////spaceranger_tutorial/AM1/
cellranger_gtf=//////datavelocyto/refdata-gex-GRCh38-2020-A/genes/genes.gtf

velocyto run10x -m $rmsk_gtf  $cellranger_outDir $cellranger_gtf

rmsk_gtf=//////datavelocyto/hg38_repeat_rmsk.gtf
cellranger_outDir=//////spaceranger_tutorial/AM2/
cellranger_gtf=//////datavelocyto/refdata-gex-GRCh38-2020-A/genes/genes.gtf

velocyto run10x -m $rmsk_gtf  $cellranger_outDir $cellranger_gtf

rmsk_gtf=//////datavelocyto/hg38_repeat_rmsk.gtf
cellranger_outDir=//////spaceranger_tutorial/AM3/
cellranger_gtf=//////datavelocyto/refdata-gex-GRCh38-2020-A/genes/genes.gtf

velocyto run10x -m $rmsk_gtf  $cellranger_outDir $cellranger_gtf
