#!/bin/bash
#SBATCH --job-name=RNAseq_CNVRG_DESeq2
#SBATCH --nodes=2
#SBATCH --time=7-00:00:00
#SBATCH --account=modelscape
#SBATCH --mem-per-cpu=20G

module load arcc/1.0 gcc/12.2.0 r/4.2.2 


R CMD BATCH /project/evolgen/emcfarl2/deer_RNAseq/scripts/Streamlined_RNAseq_comparison.R
