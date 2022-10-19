### What I need to do is 1) pull all of the RNAseq data in and 2) pull all of the admixture k=2 data in
### then, I need to put these two data sets together, and plot a PCA of how SNP_species looks. 
### Think about doing MDS rather than PCA?
### Gonna get rid of DESeq2 and EdgeR - a priori, they're not as good for an empirical study, and I want to pull out only the empirical stuff. And just write it.

### The two questions that I'm going to ask:
### 1) is there a difference in composition along the q gradient
### 2) more in line with more 'traditional DEG' studies, are there DEG between red deer and hybrids, ignoring variation in q?
### I will think more, and hard, about how to do dirichlet regression (see paper from Alex, )


library(pcaExplorer)
library(CNVRG)
library(rstan)
library(shinystan)

### will need to add something like this to it, to only look at the protein coding genes
read.table("~/Google Drive/Paper VI - RNAseq/Scottish_deer_RNAseq_R/ncbi_protein_coding_genes.tsv", header = T)->gene_list


countfiles<-list.files("~/Google Drive/Paper VI - RNAseq/Output/CountFiles_170322/Countfiles", pattern="*_count.tsv")

deernames_count<-as.factor(sapply(strsplit(countfiles, split="_"), function(x) x[1]))
tissuetype_count<-as.factor(sapply(strsplit(countfiles, split="_"), function(x) x[2]))

count<-list()

setwd("~/Google Drive/Paper VI - RNAseq/Output/CountFiles_170322/Countfiles")
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

setwd("../")

colData_1<-data.frame(cbind(as.character(deernames_count), as.character(tissuetype_count)))
names(colData_1)<-c("deername", "tissue")
read.csv("~/Google Drive/Paper VI - RNAseq/Scottish_deer_RNAseq_R/Kintyre_2019_2020.csv")->phenotypes
head(phenotypes)

merge(colData_1, phenotypes, by="deername")->merge_data
length(merge_data[,1])

### read in the ADMIXTURE analysis data
read.table("~/Google Drive/Paper VI - RNAseq/Scottish_deer_RNAseq_R/all_pops_merged.ped")->deerped
deerped[,2]->deername
read.table("~/Google Drive/Paper VI - RNAseq/Scottish_deer_RNAseq_R/all_pops_merged.2.Q")->admixture
read.table("~/Google Drive/Paper VI - RNAseq/Scottish_deer_RNAseq_R/all_pops_merged.2.Q_se")->SE
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

cnvg_data_nosika_short<-cnvg_data_nosika_ordered[,1:100]


###Start CNVRG analysis

### run the model###
modelOut <- cnvrg_HMC(countData = cnvg_data_nosika_short, 
                      starts = indexer(cnvg_data_nosika_short$treatment)$starts, #c(1,39,76,113),
                      ends = indexer(cnvg_data_nosika_short$treatment)$ends, #c(38,75,112,148),
                      chains = 2, 
                      burn = 500, 
                      samples = 1000, 
                      thinning_rate = 2,
                      cores = 1,
                      params_to_save = c("pi", "p"))

###check convergence, as one does

head(rstan::summary(modelOut, pars = "pi", probs =c(0.025, 0.975))$summary)

#shinystan::launch_shinystan(modelOut) #to look at visualizations of diagnostic parameters

point_est <- extract_point_estimate(model_out = modelOut, countData = cnvg_data_nosika_short) ## get point estimates out

### differential expression
diff_abund_test <- diff_abund(model_out = modelOut, countData = cnvg_data_nosika_short) ### this didn't work? Is it possible it only works for 2 treatments, not when I have 4? In which case, I should split up the tissues?


