#!/bin/sh
# Grid engine options
#$ -N Build_index
#$ -cwd
#$ -l h_rt=2:00:00
#$ -M eryn.mcfarlane@ed.ac.uk
#$ -m beas
#$ -l h_vmem=20G
#$ -pe sharedmem 16

#$ -o o_files/ 
#$ -e e_files/

. /etc/profile.d/modules.sh

module load igmm/apps/python/2.7.10
module load roslin/samtools/1.9
module load igmm/apps/HISAT2/2.1.0
module load igmm/apps/gffcompare/0.11.2
module load igmm/apps/FastQC/0.11.9
module load roslin/kallisto/0.44.0
module load igmm/apps/TrimGalore/0.6.6 

##need to make an index first 

hisat2-build -p 16 GCF_910594005.1_mCerEla1.1_genomic.fna PacBio_index