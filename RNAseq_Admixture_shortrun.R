### What I need to do is 1) pull all of the RNAseq data in and 2) pull all of the admixture k=2 data in
### then, I need to put these two data sets together, and plot a PCA of how SNP_species looks. 
### Think about doing MDS rather than PCA?
### Gonna get rid of DESeq2 and EdgeR - a priori, they're not as good for an empirical study, and I want to pull out only the empirical stuff. And just write it.

### The two questions that I'm going to ask:
### 1) is there a difference in composition along the q gradient
### 2) more in line with more 'traditional DEG' studies, are there DEG between red deer and hybrids, ignoring variation in q?
### I will think more, and hard, about how to do dirichlet regression (see paper from Alex, )

###changing all the filepaths to work on teton
#/project/evolgen/emcfarl2/deer_RNAseq



library(CNVRG)
library(rstan)
options(mc.cores = parallel::detectCores())
library(shinystan)
library(ggplot2)
library(ggrepel)
#library(brmstools)

### will need to add something like this to it, to only look at the protein coding genes
read.table("./datafiles/ncbi_protein_coding_genes.tsv", header = T)->gene_list


countfiles<-list.files("./countfiles", pattern="*_count.tsv")

deernames_count<-as.factor(sapply(strsplit(countfiles, split="_"), function(x) x[1]))
tissuetype_count<-as.factor(sapply(strsplit(countfiles, split="_"), function(x) x[2]))

count<-list()

setwd("../countfiles")
for(i in 1:length(countfiles)){
  count[[i]]<-read.table(countfiles[i])[1:33296,2]
}

###works to here so far...
count_df<-matrix(unlist(count), ncol=56, byrow=TRUE)

row.names(count_df)<-read.table(countfiles[1])[1:33296,1]
count_df_coding<-count_df[which(rownames(count_df) %in% gene_list$GeneID), ]

#names(count_df)<-paste0(deernames_count, tissuetype_count)
#row.names(count_df)<-read.table(countfiles[1])[1:25518,1]
count_df_coding

setwd("../../")

colData_1<-data.frame(cbind(as.character(deernames_count), as.character(tissuetype_count)))
names(colData_1)<-c("deername", "tissue")
read.csv("Kintyre_2019_2020.csv")->phenotypes
head(phenotypes)

merge(colData_1, phenotypes, by="deername")->merge_data
length(merge_data[,1])

### read in the ADMIXTURE analysis data
read.table("all_pops_merged.ped")->deerped
deerped[,2]->deername
read.table("all_pops_merged.2.Q")->admixture
read.table("all_pops_merged.2.Q_se")->SE
cbind(deername, admixture, SE)->admixture_data
names(admixture_data)<-c("EM.code", "ad_sika", "ad_red", "ad_SE2", "ad_SE1")
admixture_data$SNP_species<-ifelse((admixture_data$ad_red+1.96*admixture_data$ad_SE1>=0.99999), "red", 
                                    ifelse((admixture_data$ad_red-1.96*admixture_data$ad_SE1<=0.00001), "sika", "hybrid"))

merge(merge_data, admixture_data, by='EM.code')->alldata
alldata[,c(2,3,8,15, 18)]->colData
colData[which(colData$deername=="SAC00002831"), ]$Weight<-c(NA, NA)
colData[which(colData$Weight=="N/A"),]$Weight<-NA


###I'm going to take out the one sika individual, SAC00002777, and two sika like individuals SAC00002805 and SAC00002817

colData_nosika<-colData[which(colData$ad_red>0.5),]
colData_nosika$treatment<-paste(colData_nosika$SNP_species, colData_nosika$tissue) ### to make one treatment, 4 possibilities, two species and two tissues
treatment<-colData_nosika$treatment
count_df_coding_nosika<-count_df_coding[,c(1:22, 25:30, 33:38, 41:56)]



cnvg_data_nosika<-cbind.data.frame(treatment, t(count_df_coding_nosika))

cnvg_data_nosika_ordered<-cnvg_data_nosika[order(cnvg_data_nosika[,1]),]

###0's aren't allowed, add one to everything
cnvg_data_nosika_ordered[,2:22927] <- 1 + cnvg_data_nosika_ordered[,2:22927]

### let's make a much smaller analysis, just to start

#cnvg_data_nosika_short<-cnvg_data_nosika_ordered[,1:200]


###Start CNVRG analysis

### run the model###
modelOut <- cnvrg_HMC(countData = cnvg_data_nosika_ordered, 
                      starts = indexer(cnvg_data_nosika_ordered$treatment)$starts, 
                      ends = indexer(cnvg_data_nosika_ordered$treatment)$ends, 
                      chains = 2, 
                      burn = 25, 
                      samples = 50, 
                      thinning_rate = 2,
                      cores = 2,
                      params_to_save = c("pi", "p"))

###check convergence, as one does

##head(rstan::summary(modelOut, pars = "pi", probs =c(0.025, 0.975))$summary)
### for the short ones, Rhat is terrible
###this takes 15-20 minutes to run and do the thing.
jpeg(file="Rhat_density.jpeg")
plot(density(rstan::summary(modelOut, pars = "pi", probs =c(0.025, 0.975))$summary[,7]),  xlab = "Rhat", ylab = "Density", main = "")
dev.off()

#shinystan::launch_shinystan(modelOut) #to look at visualizations of diagnostic parameters

point_est <- extract_point_estimate(model_out = modelOut, countData = cnvg_data_nosika_ordered) ## get point estimates out

### maybe what I want from the model is the proportions of variance explained by each species and tissue? How do I get that out?

### then I can look at pairwise differences, if I really want to? I guess I want those that are different between species in muscle, but not heart!

### differential expression
diff_abund_test <- diff_abund(model_out = modelOut, countData = cnvg_data_nosika_ordered) ### this gives pairwise differences between each of the treatments, and gives the genes that are different between them

### the two differential expression tables for heart vs heart and muscle vs muscle are:
write.csv(diff_abund_test$features_that_differed$treatment_1_vs_treatment_3, file="heart_heart_DEG.csv")
write.csv(diff_abund_test$features_that_differed$treatment_2_vs_treatment_4, file="muscle_muscle_DEG.csv")


#From vingette
#This function subtracts the posterior distribution of the pi paramater for each feature in one treatment group from the pi parameter distribution in other treatment groups. 
#The function outputs a matrix of proportions that describes the proportion of the distribution of differences that is greater than zero, for each comparison. 
#In this example, the comparison between treatment 1 and treatment 3 provided the following results
##Note:
#treatment1 = hybrid heart
#treatment2 = hybrid muscle
#treatment3 = red heart
#treatment4 = red muscle

#head(diff_abund_test)

#str(diff_abund_test)
str(diff_abund_test$features_that_differed)

### this would be the manhatten plot for heart-heart differences
###make this into a dataframe that I'm happy to work with

###volcano plot of heart - heart differences
### distribution of differences from differential abundance test
mean_ests<-colMeans(diff_abund_test$ppd[[1]][[3]])
differences<-data.frame(unlist(diff_abund_test$certainty_of_diffs[2,]))
differences<-as.numeric(differences[2:22927,])
gene_names<-as.character(names(cnvg_data_nosika_ordered)[2:22927])
heart_heart<-cbind.data.frame(gene_names,mean_ests, differences)
#heart_heart$gene_names<-gene_names
heart_heart$DE<-heart_heart$gene_names %in% diff_abund_test$features_that_differed$treatment_1_vs_treatment_3$feature_that_differed
heart_heart$delabel<-NA
heart_heart$delabel<-as.character(ifelse(heart_heart$DE == "FALSE", NA, heart_heart$gene_names))
heart_heart$differences_absolute<-ifelse(heart_heart$differences>0.5, 1-heart_heart$differences, heart_heart$differences)

png(filename="heart_volcano.png", width=4,height=4,units="in",res=1200, pointsize = 1)
###still need to do labels, including that this is treatment 1 to treatment 3
plot<-ggplot(heart_heart, aes(x=mean_ests, y=-log10(differences_absolute+0.00001), colour=as.factor(DE)))+geom_point(show.legend=FALSE)+theme_minimal()
plot2<-plot+scale_color_manual(values=c("black", "red"))+xlab("Estimated difference between hybrid heart and red deer heart")+ylab("-log10 absolute certainty in differences")+theme(axis.text=element_text(size=4),axis.title=element_text(size=8))
print(plot2)
dev.off()

###volcano plot for muscle muscle###
### distribution of differences from differential abundance test
mean_ests_muscle<-colMeans(diff_abund_test$ppd[[2]][[4]])
differences_muscle<-data.frame(unlist(diff_abund_test$certainty_of_diffs[6,]))
differences_muscle<-as.numeric(differences_muscle[2:22927,])
gene_names<-as.character(names(cnvg_data_nosika_ordered)[2:22927])
muscle_muscle<-cbind.data.frame(gene_names,mean_ests_muscle, differences_muscle)
#heart_heart$gene_names<-gene_names
muscle_muscle$DE<-muscle_muscle$gene_names %in% diff_abund_test$features_that_differed$treatment_2_vs_treatment_4$feature_that_differed
muscle_muscle$delabel<-NA
muscle_muscle$delabel<-as.character(ifelse(muscle_muscle$DE == "FALSE", NA, muscle_muscle$gene_names))
muscle_muscle$differences_absolute<-ifelse(muscle_muscle$differences>0.5, 1-muscle_muscle$differences, muscle_muscle$differences)

png(filename="muscle_volcano.png", width=4,height=4,units="in",res=1200, pointsize = 1)
###still need to do labels, including that this is treatment 1 to treatment 3
plot<-ggplot(muscle_muscle, aes(x=mean_ests_muscle, y=-log10(differences_absolute+0.00001), colour=as.factor(DE)))+geom_point(show.legend=FALSE)+theme_minimal()
plot2<-plot+scale_color_manual(values=c("black", "red"))+xlab("Estimated difference between hybrid muscle and red deer muscle")+ylab("-log10 absolute certainty in differences")+theme(axis.text=element_text(size=4),axis.title=element_text(size=8))
print(plot2)
dev.off()

### this is all a mess below this ###

###bringing in some annotation
#gbff<-read.delim("~/Downloads/GCF_910594005.1_mCerEla1.1_rna.gbff", header=F, comment.char="#")

#library(rtracklayer)
#library(zoo)
library(ape)
#gtf<-rtracklayer::import("~/Downloads/GCF_910594005.1_mCerEla1.1_genomic.gtf")
gff<-read.gff("./datafiles/GCF_910594005.1_mCerEla1.1_genomic.gff")
#gff <- read.table('./datafiles/GCF_910594005.1_mCerEla1.1_genomic.gff', sep="\t", quote="")

annotation<-data.frame(cbind(gff$gene, gff$chromosome, gff@ranges@start))
annotation[,2]<-as.numeric(na.locf(annotation[,2]))

names(annotation)<-c("gene", "chromosome", "position")

annotation$chrom_pos<-as.numeric(paste0(annotation$chromosome, ".", annotation$position))

length(tapply(annotation$position, annotation$gene, min)) ###this is what I want for position



### let's try these as only significant manhatten plots? 


### probably want to plot this as a manhattan plot? volcano plot?


rstan::summary(modelOut, pars = "pi", probs =c(0.025, 0.975))$summary[,1]
rstan::c_summary(modelOut, pars = "pi", probs =c(0.025, 0.975))$summary[,1]
### let's plot some volcano plots

##these are the pvalues 'certainty'


### I don't think this is gonna work - I think that trying to put 20K genes on the xaxis is not going to work
#heart_heart<-data.frame(diff_abund_test$certainty_of_diffs)[2,]


#muscle_muscle<-data.frame(diff_abund_test$certainty_of_diffs)[6,]

#jpeg(file="muscle_muscle.jpeg")
#plot(x=reorder(colnames(muscle_muscle[2:22927]), muscle_muscle[2:22927]), y=muscle_muscle[1,2:22927], type="p")
#dev.off()



#effect size by probability?
### this might be what I want to then ask how the diversities are changing across q?
entropies <- diversity_calc(model_out = modelOut, countData = cnvg_data_nosika_ordered, entropy_measure = 'shannon',equivalents = T)
## I don't really know what this is, or what it did. I have a nasty feeling it's just telling me how many genes there are, which isn't really a relevent stat for transcriptomics. 
### the number of genes isn't going to vary the same way that the number of species might in a microbiome study.



jpeg(file="entropies.jpeg")
plot(density(entropies[[1]][[1]]), xlab = "Entropy", ylab = "Density", main = "")
dev.off()
save.images("/project/evolgen/emcfarl2/deer_RNAseq/output/deer_dirichlet_DEG_shortrunq.RData")

#### want to pull what I need for a PCA, from the old pcaExploreAttempt_170322.R code
### basically starting over, because I can't be asked to do anything other than follow the tutorial. ###

library(pcaExplorer)

### will need to add something like this to it, to only look at the protein coding genes
read.table("~/Google Drive/Paper VI - RNAseq/Scottish_deer_RNAseq_R/ncbi_protein_coding_genes.tsv", header = T)->gene_list


countfiles<-list.files("/Volumes/GoogleDrive/My Drive/Paper VI - RNAseq/Output/Countfiles_170322/Countfiles", pattern="*_count.tsv")

deernames_count<-as.factor(sapply(strsplit(countfiles, split="_"), function(x) x[1]))
tissuetype_count<-as.factor(sapply(strsplit(countfiles, split="_"), function(x) x[2]))

tissuetype_count<-as.factor(sub(".$","\\",tissuetype_count))

count<-list()

setwd("~/Google Drive/Paper VI - RNAseq/Output/CountFiles_170322/Countfiles")
for(i in 1:length(countfiles)){
  count[[i]]<-read.table(countfiles[i])[1:33296,2]
}

###works to here so far...
count_df<-matrix(unlist(count), ncol=56, byrow=TRUE)

row.names(count_df)<-read.table(countfiles[1])[1:33296,1]

colData_1<-data.frame(cbind(as.character(deernames_count), as.character(tissuetype_count)))
names(colData_1)<-c("deername", "tissue")


### this isn't working, and I don't know why ###

### use a vlookup, I really only need the SNP species for the colData, not mass. I'm probably over thinking this.

library(tidyverse)

vlookup<-function(this, data, key, value) {
  m<-match(this, data[[key]])
  data[[value]][m]
}

colData_1$SNP_species<-vlookup(colData_1$deername, alldata, "deername", "SNP_species")
colData_1$ad_red<-vlookup(colData_1$deername, alldata, "deername", "ad_red")
merge_data_nosika<-colData_1[which(colData_1$ad_red>0.5),] ###taking out everything but red deer and red deer like hybrids
merge_data_nosika->colData
as.factor(colData$SNP_species)->colData$SNP_species
as.factor(colData$tissue)->colData$tissue
as.factor(colData$deername)->colData$deername



### need to get just the same individuals that I used in CNVRG analysis ###
count_df<-count_df[,which(deernames_count %in% colData$deername)]
count_df_coding<-count_df[which(rownames(count_df) %in% gene_list$GeneID), ]

##mostly following this vingette
##http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start

dds <- DESeqDataSetFromMatrix(countData = count_df_coding,
                              colData = colData,
                              design= ~ tissue*SNP_species) ###I don't know how this is accounting for the 2 way interaction

dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
results<-results(dds, alpha=0.05)
results[which(results$padj < 0.05),]

summary(results)
plotMA(results)

vst <- vst(dds)

str(assay(vst))
sampleDists<-dist(t(assay(vst)))


plotPCA(vst, intgroup = c("SNP_species", "tissue"))+theme_bw()+scale_color_manual(values=c("orchid4", "purple4", "indianred1", "red"), name = "Species and Tissue")

#### let's do PCA just of DEG
### I don't really know what this tells me. It's still not explaining a lot of variation, even though I know these are DEG
DEG<-as.factor(c(diff_abund_test$features_that_differed$treatment_1_vs_treatment_3$feature_that_differed, diff_abund_test$features_that_differed$treatment_2_vs_treatment_4$feature_that_differed))

count_df_coding_DEG<-count_df_coding[which(rownames(count_df_coding) %in% DEG), ]

dds_DEG <- DESeqDataSetFromMatrix(countData = count_df_coding_DEG,
                              colData = colData,
                              design= ~ tissue*SNP_species)

dds_DEG <- DESeq(dds_DEG)
#keep <- rowSums(counts(dds)) >= 5
#dds <- dds[keep,]
vst_DEG <- vst(dds_DEG, nsub=nrow(dds_DEG))

str(assay(vst_DEG))
sampleDists_DEG<-dist(t(assay(vst_DEG)))


plotPCA(vst_DEG, intgroup = c("SNP_species", "tissue"))+theme_bw()+scale_color_manual(values=c("orchid4", "purple4", "indianred1", "red"), name = "Species and Tissue")



###DESeq2 for each tissue at a time

dds_heart <- DESeqDataSetFromMatrix(countData = count_df_coding[,which(colData$tissue=="heartA")],
                              colData = colData[which(colData$tissue=="heartA"),],
                              design= ~ SNP_species) ###I don't know how this is accounting for the 2 way interaction

dds_heart <- DESeq(dds_heart)
keep <- rowSums(counts(dds_heart)) >= 5
dds_heart <- dds_heart[keep,]
results<-results(dds_heart, alpha=0.05)
heart_sig<-results[which(results$padj < 0.05),]

cbind.data.frame(rownames(results), results$log2FoldChange)->DESeq_heart

### getting totally different genes that are significant. I wouldn't expect this based on Harrison's paper. 
summary(rownames(heart_sig) %in% diff_abund_test$features_that_differed$treatment_1_vs_treatment_3$feature_that_differed)




dds_muscle <- DESeqDataSetFromMatrix(countData = count_df_coding[,which(colData$tissue=="muscle")],
                                    colData = colData[which(colData$tissue=="muscle"),],
                                    design= ~ SNP_species) ###I don't know how this is accounting for the 2 way interaction

dds_muscle <- DESeq(dds_muscle)
keep <- rowSums(counts(dds_muscle)) >= 5
dds_muscle <- dds_muscle[keep,]
results<-results(dds_muscle, alpha=0.05)
muscle_sig<-results[which(results$padj < 0.05),]
cbind.data.frame(rownames(results), results$log2FoldChange)->DESeq_muscle

summary(rownames(muscle_sig) %in% diff_abund_test$features_that_differed$treatment_2_vs_treatment_4$feature_that_differed)
### I don't understand why I get completely different DEG when I use CNVRG compared to when I use DESeq2 ###

cbind.data.frame(gene_names, mean_ests, mean_ests_muscle)->CNVRG_pointestimates
names(DESeq_muscle)<-c("gene_names", "log2fold_muscle")
names(DESeq_heart)<-c("gene_names", "log2fold_heart")

merge(DESeq_muscle, DESeq_heart, by="gene_names")->DESeq2_pointestimates
merge(DESeq2_pointestimates, CNVRG_pointestimates, by="gene_names")->DESeq2_CNVRG_pointestimates

#### to get -log2 for comparison purposes
summary(ifelse(DESeq2_CNVRG_pointestimates$mean_ests!=0, -log2(DESeq2_CNVRG_pointestimates$mean_ests), -log2(DESeq2_CNVRG_pointestimates$mean_ests+0.00001)))
cor.test(DESeq2_CNVRG_pointestimates$log2fold_heart, log2(abs(DESeq2_CNVRG_pointestimates$mean_ests)), method="spearman")

cor.test(DESeq2_CNVRG_pointestimates$log2fold_muscle, log2(abs(DESeq2_CNVRG_pointestimates$mean_ests_muscle)), method="spearman")       


jpeg(file="/Volumes/GoogleDrive/My Drive/Paper VI - RNAseq/Output/muscle_correlation_CNVRG_DESeq2.jpeg")
ggplot(DESeq2_CNVRG_pointestimates, aes(log2fold_muscle, log2(abs(mean_ests_muscle))))+geom_smooth(method="lm")+geom_point()
dev.off()


jpeg(file="/Volumes/GoogleDrive/My Drive/Paper VI - RNAseq/Output/heart_correlation_CNVRG_DESeq2.jpeg")
ggplot(DESeq2_CNVRG_pointestimates, aes(log2fold_heart, log2(abs(mean_ests))))+geom_smooth(method="lm")+geom_point()
dev.off()
