#!/bin/bash
#SBATCH --job-name=RNAseq_CNVRG_short
#SBATCH --nodes=2
#SBATCH --time=7-00:00:00
#SBATCH --account=modelscape
#SBATCH --mem-per-cpu=20G

module load swset/2018.05 gcc/7.3.0 py-scipy/1.1.0 miniconda3/4.10.3 r/4.0.5-py27


R CMD BATCH /project/evolgen/emcfarl2/deer_RNAseq/scripts/RNAseq_Admixture_shortrun.R
