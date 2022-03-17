#!/bin/sh
filename="$1"
# Grid engine options
#$ -N HTseq_count
#$ -cwd
#$ -l h_rt=4:00:00
#$ -M eryn.mcfarlane@ed.ac.uk
#$ -m beas
#$ -l h_vmem=20G
#$ -pe sharedmem 16

#$ -o o_files/ 
#$ -e e_files/

. /etc/profile.d/modules.sh

module load igmm/apps/python/2.7.10

htseq-count --format bam --order pos --mode union --stranded no --minaqual 1 --type exon --idattr gene_id ./Bamfiles/"${filename}".bam GCF_910594005.1_mCerEla1.1_genomic.gtf > ./Countfiles/"${filename}"_count.tsv