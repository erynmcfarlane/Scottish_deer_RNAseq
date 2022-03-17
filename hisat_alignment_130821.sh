#!/bin/sh
filename="$1"
# Grid engine options
#$ -N Alignment
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


### right now, I don't know where these output. But the quality afterwards seems high quality. 

##trim and quality control

trim_galore --fastqc --paired /exports/csce/eddie/biology/groups/pemberton/Eryn/transcriptomics/data_in/Transcripts/"${filename}"_1.fq.gz /exports/csce/eddie/biology/groups/pemberton/Eryn/transcriptomics/data_in/Transcripts/"${filename}"_2.fq.gz --basename "${filename}"

##align to index (already built by 'hisat_index')

hisat2 -p 8 --dta --rg PL:ILLUMINA --rg-id "${filename}" -x ./PacBio_index -1 "${filename}"_val_1.fq.gz -2 "${filename}"_val_2.fq.gz -S ./Samfiles/"${filename}".sam


##go sam to bam, summarize bam file to know what's going on
samtools sort -@ 8 -o ./Bamfiles/"${filename}".bam ./Samfiles/"${filename}".sam

find ./Bamfiles/"${filename}".bam -exec echo samtools index {} \; | sh 

samtools flagstat ./Bamfiles/"${filename}".bam > ./Summaries/"${filename}"_basicsummary.txt