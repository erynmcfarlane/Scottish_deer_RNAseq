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
merge_data[,c(1,2,8,9)]->colData
colData[which(colData$deername=="SAC00002831"), ]$Weight<-c(NA, NA)
colData[which(colData$deername=="SAC00002831"), ]$Species<-c("UNK", "UNK")
colData[which(colData$Weight=="N/A"),]$Weight<-NA
colData$Species<-as.factor(colData$Species)

### this isn't working - go into DESeq2
library(DESeq2)
library(PCAtools)
dds <- DESeqDataSetFromMatrix(countData = count_df_coding,
                              colData = colData,
                              design= ~ tissue*Species) ## taking into account the two factor design?
###pull the phenotypic species from Fraser

dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vst <- vst(dds)

### from a different tutorial
#https://www.bioconductor.org/help/course-materials/2014/SeattleOct2014/B02.1.1_RNASeqLab.html#diagnostic
### start again with the vst data

str(assay(vst))
sampleDists<-dist(t(assay(vst)))

###start with a heatmap
library("gplots")
library("RColorBrewer")
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$Species, dds$tissue, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(poisd$dd)
heatmap.2( samplePoisDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )

plotPCA(vst, intgroup = "Species")


(data <- plotPCA(vst, intgroup = c( "Species", "tissue"), returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))
qplot(PC1, PC2, color=Species, shape=tissue, data=data) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

pca <- prcomp(t(assay(dds)))

PC1=pca$x[,1]
PC2=pca$x[,2]
PC3=pca$x[,3]
PC4=pca$x[,4]

qplot(PC1, PC2, color=Species, shape=tissue, data=data) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))


mds <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mds, colData(vst))
qplot(X1,X2,color=Species,shape=tissue,data=mds)
