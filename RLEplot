Rle.rco1<-import("/Users/amit/Desktop/Axel/csp1_ana/s_7_sequence_RCO1_sort.ChIP.bw",as="Rle")
gtf<-import("~/Desktop/genome/Neurospora_crassa.NC12.34.gtf")
##import fast file 
Neurospora.fasta<-FaFile("~/Desktop/genome/Neurospora_crassa.NC12.dna.toplevel.fa")
###gene of interest
a<-c("NCU02265","NCU03967")
gene<-gtf[((mcols(gtf)$gene_id %in% a))]
##gene1<-gene[((mcols(gene)$type =="type"))]
gene1<-gene[((mcols(gene)$type =="gene"))]
##last<-flank(gene1,width = 1500,both = TRUE)
last<-(gene1+3500)
names(last)<-mcols(gene1)[[5]]
##
rco1.Profiles<-Rle.rco1[last]
for (i in 1:length(rco1.Profiles)){
  if(as.vector(strand(last))[i]=="-"){
   rco1.Profiles[[i]]<-rev(rco1.Profiles[[i]])
  }
}
for(i in 1:length(rco1.Profiles)) {
pdf(file = paste(names(last)[i],".rco1_3500.pdf",sep=""))
plot(rco1.Profiles[[i]],type="l",xlab="Base",ylab="Reads",main=names(last)[i])
#axis(side=1,at=c("0","1500", length(rco1.Profiles.smooth[[i]])-1500, length(Random.NT.smooth[[i]])),labels=c("-1500","TSS","ATG","+1500"))
#axis(side=2)
#lines(Random.Ser5.smooth[[i]],type="l",col="red",lwd=5)
#legend("topleft",c("RNApolII_NT_35","RNApolII_Ser5P_35"),fill=c("black","red"))
#textxy(df$dese2, df$qpcr, labs=df$gene_name, cex=1)
dev.off()
} 
Rle.csp1<-import("/Users/amit/Desktop/Axel/csp1_ana/s_6_sequence_GTR__csp_sort.ChIP.bw",as="Rle")
csp1.Profiles<-Rle.csp1[last]
for (i in 1:length(csp1.Profiles)){
  if(as.vector(strand(last))[i]=="-"){
   csp1.Profiles[[i]]<-rev(csp1.Profiles[[i]])
  }
}

####

for(i in 1:length(csp1.Profiles)) {
pdf(file = paste(names(last)[i],".csp1.Profiles_3500.pdf",sep=""))
plot(csp1.Profiles[[i]],type="l",xlab="Base",ylab="Reads",main=names(last)[i])
#axis(side=1,at=c("0","1500", length(rco1.Profiles.smooth[[i]])-1500, length(Random.NT.smooth[[i]])),labels=c("-1500","TSS","ATG","+1500"))
#axis(side=2)
#lines(Random.Ser5.smooth[[i]],type="l",col="red",lwd=5)
#legend("topleft",c("RNApolII_NT_35","RNApolII_Ser5P_35"),fill=c("black","red"))
#textxy(df$dese2, df$qpcr, labs=df$gene_name, cex=1)
dev.off()
} 
rco1.Profiles.smooth<-S4Vectors::runmean(rco1.Profiles,101,endrule ="constant")
for(i in 1:length(rco1.Profiles.smooth)) {
pdf(file = paste(names(last)[i],".plo2profile.pdf",sep=""))
plot(rco1.Profiles.smooth[[i]],type="l",xlab="Base",ylab="Reads",main=names(last)[i])
#axis(side=1,at=c("0","1500", length(rco1.Profiles.smooth[[i]])-1500, length(Random.NT.smooth[[i]])),labels=c("-1500","TSS","ATG","+1500"))
#axis(side=2)
#lines(Random.Ser5.smooth[[i]],type="l",col="red",lwd=5)
#legend("topleft",c("RNApolII_NT_35","RNApolII_Ser5P_35"),fill=c("black","red"))
#textxy(df$dese2, df$qpcr, labs=df$gene_name, cex=1)
dev.off()
} 


