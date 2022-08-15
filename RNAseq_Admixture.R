### What I need to do is 1) pull all of the RNAseq data in and 2) pull all of the admixture k=2 data in
### then, I need to put these two data sets together, and plot a PCA of how SNP_species looks. 
### Think about doing MDS rather than PCA?
### First steps can be DESeq2, just as we did in the flycatchers, using the SNP_species
### Second step could be using limma-voom
### Third step - Harrison paper
### Fourth step - corncob loop? 


### only the corn cob loop will account for variation in the admixture estimates, I think. 

###hoping for some quick and dirty PCAs for a 'here's what I'm doing talk' coming up
###famous last words, I suspect I won't be able to do anything.

library(pcaExplorer)

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


### this isn't working - go into DESeq2
library(DESeq2)
library(PCAtools)
dds <- DESeqDataSetFromMatrix(countData = count_df_coding,
                              colData = colData, 
                              design=  ~ as.factor(tissue)*as.factor(SNP_species)) 

dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vst <- vst(dds)

str(assay(vst))
sampleDists<-dist(t(assay(vst)))


###start with a heatmap
library("gplots")
library("RColorBrewer")
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$SNP_species, dds$tissue, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(poisd$dd)
heatmap(samplePoisDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )

plotPCA(vst, intgroup = "SNP_species")


(data <- plotPCA(vst, intgroup = c( "SNP_species", "tissue"), returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))
qplot(PC1, PC2, color=SNP_species, shape=tissue, data=data) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

pca <- prcomp(t(assay(dds)))

PC1=pca$x[,1]
PC2=pca$x[,2]
PC3=pca$x[,3]
PC4=pca$x[,4]

qplot(PC1, PC2, color=SNP_species, shape=tissue, data=data) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))


mds <- data.frame(cmdscale(samplePoisDistMatrix))
mds <- cbind(mds, colData(vst))
qplot(X1,X2,color=SNP_species,shape=tissue,data=mds)


###rlog because the data are heterschedastic

#rld<-rlog(dds) ### vst is already doing this - don't need this!

### let's do some actual differential expression, at least between tissues?
#dds <- DESeq(dds) - already did this, this is the analysis
results<-results(dds) ###not really sure what this is doing, there's some pairwise thing that isn't working. 
### just for point of view, let's use just the muscle (where we expect the difference) and look at red like and sika liek animals. 
colData$species_like<-ifelse(colData$ad_red>0.95, 'red like', ifelse(colData$ad_red<0.5, 'sika like', 'mid hybrid'))
plot_2A<-ggplot(colData, aes(x=reorder(deername, ad_red), y=ad_red, colour=SNP_species))+geom_point(size=1)
plot_2A<-plot_2A+labs(x="Deer", y="Proportion Red Deer")+scale_colour_manual(values=c("purple", "red", "blue"))+theme(panel.grid=element_blank(), panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), legend.key=element_blank())
plot2A<-plot_2A+theme(axis.line=element_line(colour="black", size=0.5))
plot2A<-plot_2A+geom_hline(yintercept=c(0, 0.05, 0.95, 1), colour="grey")
plot2A+guides(colour=guide_legend(title = "SNP species"))

#### Things that I've learned - most of the animals are very red deer like hybrids. Might be better to take out the 3 sika like individuals, and only look at the other deer?
colData[which(colData$ad_red>0.5), ]->colData_2
count_df_coding[,c(as.numeric(unlist(rownames(colData_2))))]->count_df_coding_2
dds_2 <- DESeqDataSetFromMatrix(countData = count_df_coding_2,
                              colData = colData_2, 
                              design=  ~ tissue+ad_red:tissue) 

dds <- DESeq(dds_2)
res<-results(dds)
#https://support.bioconductor.org/p/126713/
#the interpretation of the log2 fold change is the fold change in expression for one unit of the variable

res[which(res$padj<0.05),]

### my feeling is that each tissue should be done separtately, that's how deseq might work?



##### Step 2 using Limma-Voom

## need to transpose the count data, and stick it on the coldata?
### let's do the voom step first. This 'Transform count data to log2-counts per million (logCPM), estimate the mean-variance relationship and use this to compute appropriate observation-level weights. The data are then ready for linear modelling.'
### following this tutorial
library(edgeR)

counts_edgeR<-count_df_coding
colnames(counts_edgeR)<-as.factor(paste0(colData$deername, colData$tissue))
do<-DGEList(counts_edgeR)
d0<-calcNormFactors(do)
cutoff<-3
drop<-which(apply(cpm(d0), 1, max)<cutoff)
d<-d0[-drop,]
dim(d) ### number of genes left
snames<-colData$deername
species<-colData$SNP_species
tissue<-colData$tissue

group<-interaction(species, tissue)

plotMDS(d, col=as.numeric(group))


mm<-model.matrix(~0+group)
y<-voom(d, mm, plot=T) ### this gives a negative quadratic shape, which sounds bad? stack over flow basically says 'do more filtering!'
### https://stats.stackexchange.com/questions/160255/voom-mean-variance-trend-plot-how-to-interpret-the-plot
### I  don't entirely understand what this plot is showing, but it does seem like we want a positive quadratic, not a negatic one. 
<<<<<<< HEAD
### another tutorial with more explanation of the mean-variance trend https://f1000research.com/articles/5-1408
#### is the variance going down at the end because there are so few genes that are so highly expressed?
### It also seems like a marked pattern like this is evidence of lots of biological variance, which is good, that's kind of what we're going for. 
=======

>>>>>>> 544282bc5d2ac935165bb4b891487f850b8b0f24
