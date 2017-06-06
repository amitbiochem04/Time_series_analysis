

library(biomaRt)

##get all the 3'utrse
##do-d1

a<-read.delim("do_d1_up.csv",header=T,sep=",")
eg<-a$isoform_id  
ensembl=useMart("ensembl",dataset="drerio_gene_ensembl")
utr <- getSequence(id=eg, type="ensembl_transcript_id", seqType="3utr", mart=ensembl)
outfile <- file("do_d1_up.fa", "w")
for (i in 1:nrow(utr)) {
    h = paste(c(">", utr[i,2]), collapse="")
    writeLines(h, outfile)
    writeLines(utr[i,1], outfile)
}
close(outfile)


####for human 

a<-read.delim("d2-d4_up_homer.csv",header=T,sep=",")
a$ensembl<-a$isoform_id
genemap <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene","ensembl_transcript_id","refseq_mrna"),filters = "ensembl_gene_id",values = a$ensembl,mart = ensembl)
idx <- match(a$ensembl, genemap$ensembl_gene_id)
a$transc <- genemap$ensembl_transcript_id[idx]
a$ref<- genemap$refseq_mrna[idx]
eg<-a$ensembl
utr <- getSequence(id=eg, type="ensembl_transcript_id", seqType="3utr", mart=ensembl)
outfile <- file("d2-d4_up_3utr.fa", "w")
for (i in 1:nrow(utr)) {
    h = paste(c(">", utr[i,2]), collapse="")
    writeLines(h, outfile)
    writeLines(utr[i,1], outfile)
}
close(outfile)



fungi = useMart(biomart="fungal_mart",host="fungi.ensembl.org",dataset="NC12")

fungi = useMart(biomart="fungal_mart",host="fungi.ensembl.org",dataset="ncrassa_eg_gene")
ensembl_id="NCU15815"

gb <-getBM(attributes=c("ensembl_gene_id","external_gene_id","entrezgene"),filters = "ensembl_gene_id", values=ensembl_id, mart=fungi)  
attributes <- listAttributes(fungi)
head(attributes)

ensembl_id=c("NCU11424","NCU15815")

gb <-getBM(attributes=c("ensembl_gene_id","external_gene_name"),filters = "ensembl_gene_id", values=id$name, mart=fungi)






