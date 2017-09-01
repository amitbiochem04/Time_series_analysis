library(DESeq2)
library(DESeq2)
library(dplyr)
library(reshape2)
library(annotate)
library(scales)
library(calibrate)
library(openxlsx)
CairoFonts(
  regular="Arial:style=Regular",
  bold="Arial:style=Bold",
  italic="Arial:style=Italic",
  bolditalic="Arial:style=Bold Italic,BoldItalic",
  symbol="Symbol")

#####load the meta file 
setwd("/Users/amitsingh/Desktop/Light-dark/")
meta<-read.delim("meta_1.txt",header=T,sep="\t")
sampleTable<-data.frame(sampleName=meta$condition, fileName=meta$fileName, condition=meta$sampletype,time=meta$time,sampletype=meta$sampletype,samplelabel=meta$samplelabel)
#### load full data 
fullData <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory = ".",design = ~ 1)
fullData <- estimateSizeFactors(fullData)
####normalize the fullData
norm<-counts(fullData, normalized=TRUE)
###save normalized data
save(fullData,file="fullData.rda")
#########plot distance matrix 
#sampleDists <- dist(t(log2(norm+1)))
#sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix)<-paste(colnames(norm), sep="-")
#CairoPDF(file="samaple_cluster",width=15/1.54,height=15/1.54)
#pheatmap(sampleDistMatrix, cluster_col=F,cluster_row=F)
#dev.off()

CairoPDF(file="samaple_PCA",width=10/1.54,height=10/1.54)

####
#mymat<-as.matrix(cnts.norm[!grepl('NA', rownames(cnts.norm)), ])
#p5 <- labdsv::pca(t(mymat),dim=4,cor=F)
#mysum <- summary(p5)
#plot(p5$score[,1:2],pch=19,ylab=paste("PC2 (",round(100*mysum[2,2],2),"%)",sep=""),type="n", xlab=paste("PC1 (",round(100*mysum[2,1],2),"%)",sep=""))
#loads <- as.data.frame(p5$loadings[,1:2])
#loads$sym<-rownames(loads)
#loads$len <- sqrt(loads[,1]^2+loads[,2]^2)
#loads <- loads[order(abs(loads[,4]),decreasing=T),]
#num <- 5
#arrows(0,0,200*loads[1:num,1],200*loads[1:num,2],col="grey",length=1)
#text(200*loads[1:num,1],200*loads[1:num,2],loads[1:num,3],cex=0.5,col="blue")
#text(p5$score[,1],p5$score[,2],rownames(p5$score),pos=3)
#lines(p5$score[1:14,1],p5$score[1:14,2],col="red",type="o",pch=19)
#lines(p5$score[15:28,1],p5$score[15:28,2],col="green",type="o",pch=19)
#lines(p5$score[29:41,1],p5$score[29:41,2],col="blue",type="o",pch=19)
#lines(p5$score[42:54,1],p5$score[42:54,2],col="black",type="o",pch=19)
#lines(p5$score[54:60,1],p5$score[54:60,2],col="yellow",type="o",pch=19)
#legend("topleft",c("WT","DFrq","D_rep1","Drep2","Drep3"),fill=c("red","green","blue","black","yellow"),inset=0.02)
#write.xlsx(data.frame(loads), "pca_load.xlsx")
#dev.off()
###simple plot 
##plot(p5$scores[, 1], p5$scores[, 2])
#text(p5$scores[, 1], p5$scores[, 2], labels=rownames(p5$scores), cex= 0.7)
 


#####contrast between light and dark 
dds <- fullData[,fullData$samplelabel=="dl-wt" & fullData$time %in% c(4,6,8,10,11.5)]
design(dds) <- ~ condition
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "light", "dark"))
ressig<- subset(res, padj <0.1)
up<-subset(res, padj <0.1 & log2FoldChange >0)
down<-subset(res, padj <0.1 & log2FoldChange <0)
ressig<-subset(res, padj < 0.1& abs(log2FoldChange)>1))
#####
col_1<-alpha("gray50", 0.5)
col_2<-alpha("red", 0.5)
res<-na.omit(res)
CairoPDF(file="MA.plotlight.pdf",width=10/1.54,height=10/1.54)
with(res, plot(baseMean, log2FoldChange, pch = 16, cex = 0.45, log = "x",xlab = "mean of normalized counts", ylab=expression(log[2]~fold~change),main="light inducible gene",col=col_1))
with(subset(res, padj < 0.1& abs(log2FoldChange)>1), points(baseMean, log2FoldChange, col = col_2,pch = 16,cex = 0.45))
axis (2, at = c(-1, 1, axis (2)))
abline(a=1, b=0, col=alpha ("black", 0.9), lty=2)
abline(a=-1, b=0, col=alpha ("black", 0.9), lty=2)
marked_genes<-c(rownames(ressig)[1:5])
#marked_genes<-c(head(res$name)[1:3],tail(res$name)[1:3])
marked_genes <- res[rownames(res) %in% marked_genes,]
points(marked_genes$baseMean, marked_genes$log2FoldChange, col="steelblue4", cex = 0.45, pch = 16)
with (marked_genes, text (baseMean, log2FoldChange, labels = rownames(marked_genes), pos = 2, col="steelblue4", cex = 0.5))

#legend("topleft", legend = c("FDR<0.1", "|LFC|>1"), pch = 16, col = c("col_1","col_2"))
dev.off()

anno<-read.delim("Ncrassa.genes.FungiDB.txt",header=T,sep="\t")
anno<-anno[,1:2]
###

keep<-list()
for (i in rownames(ressig)){
  idx<-grep(i,anno$X.Gene.ID.)
  idx<-anno[idx,]
  keep[[i]]<-idx
} 
df<-melt(keep,id=c("X.Gene.ID.","X.Gene.Name.or.Symbol."))
df<-df[,1:2]
colnames(df)<-c("ens.id","symbol")
#df<-df[!duplicated(df), ]
ressig<-as.data.frame(cbind(ressig,df[!duplicated(df), ]))
####save the file 
write.xlsx(ressig,file="Light_Indicube.xlsx",asTable = TRUE)

