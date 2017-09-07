library( hwriter )
###get the motif pattern form  from http://prosite.expasy.org/PS00028
##use the sequnce to creat motif.meme form online fimo webpage
##follow following command
fimo --oc fimo_overlap_mus_tans --parse-genomic-coord motifs.meme ~/Desktop/genome/Neurospora_crassa.NC12.pep.all.fa
fimo<-read.delim("fimo.txt",header=T,sep="\t")
anno<-read.delim("all_gene_name.txt",header=T,sep="\t")
last<-anno[,c(10:11,17)]
#idx<-unique(fimo$sequence.name)
#idx<-duplicated(fimo$sequence.name)
#fimo<-fimo[idx,]
fimo<-fimo[!duplicated(fimo$sequence.name), ]
row.names(fimo)<-fimo$sequence.name
#idx<-unique(t$sequence.name)
#fimo<-t[idx,]
idx<-last[last$transcript_id %in%fimo$sequence.name,]
last<-idx[!duplicated(idx$transcript_id),]
row.names(final)<-final$transcript_id
sig <- merge(fimo, final, by=0, all=TRUE)
sig<-sig[,c(3:10,12,11)]
###write in html file 
maintable <- hwrite(sig, 
                     table.style = "border-collapse: collapse;",
                     cellpadding = 5,row.bgcolor='#ffdc98',
                     center=TRUE,onmouseover="this.bgColor='#ffaaaa'", 
                     onmouseout="this.bgColor='white'", bgcolor='white')

hwrite(maintable, page="table.html")


