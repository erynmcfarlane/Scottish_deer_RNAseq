###hoping for some quick and dirty PCAs for a 'here's what I'm doing talk' coming up
###famous last words, I suspect I won't be able to do anything.

library(pcaExplorer)

### will need to add something like this to it, to only look at the protein coding genes
read.table("~/Google Drive/Paper VI - RNAseq/Scottish_deer_RNAseq_R/ncbi_protein_coding_genes.tsv", header = T)->gene_list


##read in all of the data, in heart and muscle separately
###heart

heart_files<-list.files("~/Google Drive/Paper VI - RNAseq/Output/CountFiles_170322/Countfiles", pattern="*heartA_count.tsv")

deernames_heart<-as.factor(sapply(strsplit(heart_files, split="_"), function(x) x[1]))

heart<-list()

setwd("~/Google Drive/Paper VI - RNAseq/Output/CountFiles_170322/Countfiles")
for(i in 1:length(heart_files)){
  heart[[i]]<-read.table(heart_files[i])[1:33296,2]
}

heart_df<-data.frame(matrix(unlist(heart), ncol=length(deernames_heart), byrow=TRUE), stringsAsFactors = FALSE)
names(heart_df)<-deernames_heart
row.names(heart_df)<-read.table(heart_files[1])[1:33296,1]
heart_df_coding<-heart_df[which(rownames(heart_df) %in% gene_list$GeneID), ]
setwd("../")

###muscle

muscle_files<-list.files("~/Google Drive/Paper VI - RNAseq/Output/CountFiles_170322/Countfiles", pattern="*muscleA_count.tsv")

deernames_muscle<-as.factor(sapply(strsplit(muscle_files, split="_"), function(x) x[1]))

muscle<-list()

setwd("~/Google Drive/Paper VI - RNAseq/Output/CountFiles_170322/Countfiles")
for(i in 1:length(muscle_files)){
  muscle[[i]]<-read.table(muscle_files[i])[1:33296,2]
}

###works to here so far...
muscle_df<-data.frame(matrix(unlist(muscle), ncol=length(deernames_muscle), byrow=TRUE), stringsAsFactors = FALSE)
names(muscle_df)<-deernames_muscle
row.names(muscle_df)<-read.table(muscle_files[1])[1:33296,1]
muscle_df_coding<-muscle_df[which(rownames(muscle_df) %in% gene_list$GeneID), ]


setwd("../")


### just kidding, maybe I want everything in one count matrix!

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
count_df

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
dds <- DESeqDataSetFromMatrix(countData = count_df,
                              colData = colData,
                              design= ~ tissue*Species) ## taking into account the two factor design?
###pull the phenotypic species from Fraser

dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vst <- vst(dds)
p <- pca(vst)
biplot(p)
screeplot(p, axisLabSize = 18, titleLabSize = 22)

resultsNames(dds)


### from a different tutorial
#https://www.bioconductor.org/help/course-materials/2014/SeattleOct2014/B02.1.1_RNASeqLab.html#diagnostic
### start again with the vst data

head(assay(vst))
sampleDists<-dist(t(assay(vst)))

###start with a heatmap
library("gplots")
library("RColorBrewer")
rownames(sampleDistMatrix)<-paste(vst)

sampleDistMatrix <- as.matrix(sampleDists)

#rownames(sampleDistMatrix)<-paste(vst$tissue, vst$Species, sep="-")
rownames(sampleDistMatrix)<-vst$Species
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

hc <- hclust(sampleDists)

heatmap.2(sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )

plotPCA(vst, intgroup = "Species")
