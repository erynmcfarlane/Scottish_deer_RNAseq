#### packages ####

library(CNVRG)
library(rstan)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(tidyverse)

vlookup<-function(this, data, key, value) {
  m<-match(this, data[[key]])
  data[[value]][m]
}

#### Load in all data and get it ready for the different analyses ####


### will need to add something like this to it, to only look at the protein coding genes
read.table("./datafiles/ncbi_protein_coding_genes.tsv", header = T)->gene_list

countfiles<-list.files("./countfiles", pattern="*_count.tsv")

deernames_count<-as.factor(sapply(strsplit(countfiles, split="_"), function(x) x[1]))
tissuetype_count<-as.factor(sapply(strsplit(countfiles, split="_"), function(x) x[2]))
tissuetype_count<-as.factor(sub(".$","\\",tissuetype_count))

count<-list()

setwd("./countfiles")
for(i in 1:length(countfiles)){
  count[[i]]<-read.table(countfiles[i])[1:33296,2]
}

###works to here so far...
count_df<-matrix(unlist(count), ncol=56, byrow=TRUE)

row.names(count_df)<-read.table(countfiles[1])[1:33296,1]
count_df_coding<-count_df[which(rownames(count_df) %in% gene_list$GeneID), ]

#names(count_df)<-paste0(deernames_count, tissuetype_count)
#row.names(count_df)<-read.table(countfiles[1])[1:25518,1]
str(count_df_coding) ### This can be used both for cnvrg and deseq2, I think

setwd("../")


### Pull in Admixture and Phenotype data ###
colData_1<-data.frame(cbind(as.character(deernames_count), as.character(tissuetype_count)))
names(colData_1)<-c("deername", "tissue")
read.csv("./datafiles/Kintyre_2019_2020.csv")->phenotypes
head(phenotypes)

merge(colData_1, phenotypes, by="deername")->merge_data
length(merge_data[,1])

### read in the ADMIXTURE analysis data
read.table("./datafiles/all_pops_merged.ped")->deerped
deerped[,2]->deername
read.table("./datafiles/all_pops_merged.2.Q")->admixture
read.table("./datafiles/all_pops_merged.2.Q_se")->SE
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
colData_nosika$treatment<-paste0(colData_nosika$SNP_species, colData_nosika$tissue) ### to make one treatment, 4 possibilities, two species and two tissues
treatment<-colData_nosika$treatment
count_df_coding_nosika<-count_df_coding[,c(1:22, 25:30, 33:38, 41:56)]### This can be used both for cnvrg and deseq2, I thin

## let's test it short for now, comment this out later
count_df_coding_nosika<-count_df_coding_nosika[1:1000,]

#### CNVRG Analyses ####
### data wranging specific for cnvg ###
cnvg_data_nosika<-cbind.data.frame(treatment, t(count_df_coding_nosika))
cnvg_data_nosika_ordered<-cnvg_data_nosika[order(cnvg_data_nosika[,1]),]

###This is something that I'm not at all sure about, is this appropriate to keep for all analyses?
cnvg_data_nosika_ordered[,2:length(cnvg_data_nosika_ordered)] <- 1 + cnvg_data_nosika_ordered[,2:length(cnvg_data_nosika_ordered)] ### I haven't done this for the DESeq2 analyses, how much does it matter?

## Run the actual analysis

modelOut <- cnvrg_HMC(countData = cnvg_data_nosika_ordered, 
                      starts = indexer(cnvg_data_nosika_ordered$treatment)$starts, 
                      ends = indexer(cnvg_data_nosika_ordered$treatment)$ends, 
                      chains = 2, 
                      burn = 500, 
                      samples = 1000, 
                      thinning_rate = 2,
                      cores = 10,
                      params_to_save = c("pi", "p"))

##check for convergence ##
jpeg(file="Rhat_density.jpeg")
plot(density(rstan::summary(modelOut, pars = "pi", probs =c(0.025, 0.975))$summary[,7]),  xlab = "Rhat", ylab = "Density", main = "")
dev.off()

point_est <- extract_point_estimate(model_out = modelOut, countData = cnvg_data_nosika_ordered) ## get point estimates out
### differential expression
diff_abund_test <- diff_abund(model_out = modelOut, countData = cnvg_data_nosika_ordered) ### this gives pairwise differences between each of the treatments, and gives the genes that are different between them

### the two differential expression tables for heart vs heart and muscle vs muscle are:
write.csv(diff_abund_test$features_that_differed$treatment_1_vs_treatment_3, file="heart_heart_DEG.csv")
write.csv(diff_abund_test$features_that_differed$treatment_2_vs_treatment_4, file="muscle_muscle_DEG.csv")
## CNVRG Volcano Plots ##

###volcano plot of heart - heart differences
### distribution of differences from differential abundance test
mean_ests<-colMeans(diff_abund_test$ppd[[1]][[3]])
differences<-data.frame(unlist(diff_abund_test$certainty_of_diffs[2,]))
differences<-as.numeric(differences[2:length(cnvg_data_nosika_ordered),])
gene_names<-as.character(names(cnvg_data_nosika_ordered)[2:length(cnvg_data_nosika_ordered)])
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
differences_muscle<-as.numeric(differences_muscle[2:length(cnvg_data_nosika_ordered),])
gene_names<-as.character(names(cnvg_data_nosika_ordered)[2:length(cnvg_data_nosika_ordered)])
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


#### DESeq2 Analyses ####

data.frame(treatment)->treatment.df
as.factor(treatment.df$treatment)->treatment.df$treatment

dds <- DESeqDataSetFromMatrix(countData = count_df_coding_nosika,
                              colData = treatment.df,
                              design= ~ treatment)

dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

results_heart<-results(dds, alpha=0.05, contrast=c("treatment", "hybridheart", "redheart"))
results_heart[which(results_heart$padj < 0.05),]

results_muscle<-results(dds, alpha=0.05, contrast=c("treatment", "hybridmuscle", "redmuscle"))
results_muscle[which(results_muscle$padj < 0.05),]

#### comparison between CNVRG and DESeq2 ####
cbind.data.frame(rownames(results_heart), results_heart$log2FoldChange, results_muscle$log2FoldChange)->DESeq_results
names(DESeq_results)<-c("gene_names", "heart_DESeq_log2fold", "muscle_DESeq_log2fold")

cbind.data.frame(gene_names, mean_ests, mean_ests_muscle)->CNVRG_pointestimates
merge(DESeq_results, CNVRG_pointestimates, by="gene_names")->DESeq2_CNVRG_pointestimates

cor.test(DESeq2_CNVRG_pointestimates$heart_DESeq_log2fold, DESeq2_CNVRG_pointestimates$mean_ests) ### still not at all correlated, really.

cor.test(DESeq2_CNVRG_pointestimates$muscle_DESeq_log2fold, DESeq2_CNVRG_pointestimates$mean_ests)

jpeg(file="heart_correlation.jpeg")
plot(DESeq2_CNVRG_pointestimates$heart_DESeq_log2fold, DESeq2_CNVRG_pointestimates$mean_ests)
dev.off()

jpeg(file="muscle_correlation.jpeg")
plot(DESeq2_CNVRG_pointestimates$muscle_DESeq_log2fold, DESeq2_CNVRG_pointestimates$mean_ests_muscle)
dev.off()

##this doesn't work, but also, do I want it?
#png(filename="PCA_DESeq2.png", width=4,height=4,units="in",res=1200, pointsize = 1)
#PCA<-plotPCA(vst, intgroup = c("treatment"))+theme_bw()+scale_color_manual(values=c("orchid4", "purple4", "indianred1", "red"), name = "Species and Tissue")
#print(PCA)
#dev.off()






