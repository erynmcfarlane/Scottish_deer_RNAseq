###hoping for some quick and dirty PCAs for a 'here's what I'm doing talk' coming up
###famous last words, I suspect I won't be able to do anything.

library(pcaExplorer)

##read in all of the data, in heart and muscle separately
###heart

heart_files<-list.files("./TSV_files/", pattern="*heartA_CerEla1.1.tsv")

deernames_heart<-as.factor(sapply(strsplit(heart_files, split="_"), function(x) x[1]))

heart<-list()

setwd("./TSV_files/")
for(i in 1:length(heart_files)){
  heart[[i]]<-read.table(heart_files[i])[1:25518,2]
}

###works to here so far...
heart_df<-data.frame(matrix(unlist(heart), ncol=length(deernames_heart), byrow=TRUE), stringsAsFactors = FALSE)
names(heart_df)<-deernames_heart
row.names(heart_df)<-read.table(heart_files[1])[1:25518,1]

setwd("../")

###muscle

muscle_files<-list.files("./TSV_files/", pattern="*muscleA_CerEla1.1.tsv")

deernames_muscle<-as.factor(sapply(strsplit(muscle_files, split="_"), function(x) x[1]))

muscle<-list()

setwd("./TSV_files/")
for(i in 1:length(muscle_files)){
  muscle[[i]]<-read.table(muscle_files[i])[1:25518,2]
}

###works to here so far...
muscle_df<-data.frame(matrix(unlist(muscle), ncol=length(deernames_muscle), byrow=TRUE), stringsAsFactors = FALSE)
names(muscle_df)<-deernames_muscle
row.names(muscle_df)<-read.table(muscle_files[1])[1:25518,1]
muscle_df

setwd("../")


### just kidding, maybe I want everything in one count matrix!

countfiles<-list.files("./TSV_files/", pattern="*CerEla1.1.tsv")

deernames_count<-as.factor(sapply(strsplit(countfiles, split="_"), function(x) x[1]))
tissuetype_count<-as.factor(sapply(strsplit(countfiles, split="_"), function(x) x[2]))

count<-list()

setwd("./TSV_files/")
for(i in 1:length(countfiles)){
  count[[i]]<-read.table(countfiles[i])[1:25518,2]
}

###works to here so far...
count_df<-matrix(unlist(count), ncol=56, byrow=TRUE)
#names(count_df)<-paste0(deernames_count, tissuetype_count)
#row.names(count_df)<-read.table(countfiles[1])[1:25518,1]
count_df

setwd("../")

colData<-data.frame(cbind(as.character(deernames_count), as.character(tissuetype_count)))
names(colData)<-c("deername", "tissue")


### this isn't working - go into DESeq2
library(DESeq2)
library(PCAtools)
dds <- DESeqDataSetFromMatrix(countData = count_df,
                              colData = colData,
                              design= ~ deername + tissue)
###pull the phenotypic species from Fraser

dds <- DESeq(dds)
vst <- assay(vst(dds))
p <- pca(vst, metadata = colData, removeVar = 0.1)
biplot(p)
screeplot(p, axisLabSize = 18, titleLabSize = 22)

