###batch correction
library(edgeR)
library(DESeq2)
setwd("/Users/amitsingh/Desktop/Light-dark/")
meta<-read.delim("meta_1.txt",header=T,sep="\t")
sampleTable<-data.frame(sampleName=meta$condition, fileName=meta$fileName, condition=meta$sampletype,time=meta$time,sampletype=meta$sampletype,samplelabel=meta$samplelabel)
#### load full data 
fullData <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory = ".",design = ~ 1)
norm<-counts(fullData)
dark<-norm[,c(29:52)]
y <- DGEList(counts=dark)
#y$batch<-batch
##here two library
ba<-factor(c(rep(1,times=12),rep(0,times=12)))
####Filter non-expressed genes
A <- aveLogCPM(y)
y2 <- y[A>1,]
##Then normalize and compute log2 counts-per-million
y2 <- calcNormFactors(y2)
CPM <- cpm(y2, log=TRUE, prior.count=5)
##Then remove batch correct
CPMc <- removeBatchEffect(CPM, batch)
#####plot before Batch correction 
sampleDists <- dist(t(CPM))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix)<-paste(colnames(CPM), sep="-")
pheatmap(sampleDistMatrix, cluster_col=F,cluster_row=F)
#####plot after Batch correction
sampleDists <- dist(t(CPMc))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix)<-paste(colnames(CPMc), sep="-")
pheatmap(sampleDistMatrix, cluster_col=F,cluster_row=F)

