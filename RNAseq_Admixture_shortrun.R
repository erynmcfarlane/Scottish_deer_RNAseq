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

setwd("../")

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
#str(diff_abund_test$features_that_differed)

### want to look at probability of differences either >0.95 or <0.05
#diff_abund_test$features_that_differed$treatment_1_vs_treatment_3 ##heart vs heart
#diff_abund_test$features_that_differed$treatment_2_vs_treatment_4 ## muscle vs muscle


### probably want to plot this as a manhattan plot? volcano plot?



rstan::summary(modelOut, pars = "pi", probs =c(0.025, 0.975))$summary[,1]
rstan::c_summary(modelOut, pars = "pi", probs =c(0.025, 0.975))$summary[,1]
### let's plot some volcano plots

##these are the pvalues 'certainty'



heart_heart<-data.frame(diff_abund_test$certainty_of_diffs)[2,]

jpeg(file="heart_heart.jpeg")
plot(x=reorder(colnames(heart_heart[2:22927]), heart_heart[2:22927]), y=heart_heart[1,2:22927], type="p")
dev.off()

muscle_muscle<-data.frame(diff_abund_test$certainty_of_diffs)[6,]

#effect size by probability?
### this might be what I want to then ask how the diversities are changing across q?
entropies <- diversity_calc(model_out = modelOut, countData = cnvg_data_nosika_ordered, entropy_measure = 'shannon',equivalents = T)
## I don't really know what this is, or what it did. I have a nasty feeling it's just telling me how many genes there are, which isn't really a relevent stat for transcriptomics. 
### the number of genes isn't going to vary the same way that the number of species might in a microbiome study.



jpeg(file="entropies.jpeg")
plot(density(entropies[[1]][[1]]), xlab = "Entropy", ylab = "Density", main = "")
dev.off()
save.images("/project/evolgen/emcfarl2/deer_RNAseq/output/deer_dirichlet_DEG_shortrunq.RData")


###think a little bit about what you want out of this model?
### Do I want a plot with the number of DEG, so pairwise comparison on the x, pairwise genes on the y?
### some kind of heatmap for each pairwise, where deeper colours are larger differences?
