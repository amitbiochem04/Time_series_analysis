library(dplyr)
library(ggplus)
library(tidyr)
#####plot all ressignificant gene in one pdf file 
pdf('light_inducibale_gene.pdf')
p<- data.frame(colData(fullData),t(counts(fullData, normalized=TRUE)[rownames(ressig),]))[-58,]%>%
gather( "gene", "ncount", -(1:5) ) %>%
ggplot + geom_line(aes( x=time, y=log2(ncount+1), col=samplelabel)) +facet_wrap(~gene)
p<- p + xlab( "Time[h]") + ylab("Log2(count)+1") + theme_bw()
p<-p+geom_point(aes(x=time,y=log2(ncount+1),col=samplelabel))
facet_multiple( plot = p, facets = 'gene', ncol = 4, nrow = 4)
dev.off()

######plot single gene to check 
p<-data.frame( colData(fullData),t(counts(fullData, normalized=TRUE)["NCU00552", , drop=FALSE]))[-58,]%>%
gather( "gene", "ncount", -(1:5) )%>%
ggplot +geom_line(aes(x=time, y=log2(ncount+1), col=samplelabel)) + facet_wrap(~gene)
p+ xlab("Time[h]") + ylab("Log2(count)+1") + theme_bw()

###### plot in one gene in one png file, removed one sample for fulldata
for( i in row.names(ressig))
 {
  p<-data.frame(colData(fullData),t(counts(fullData, normalized=TRUE)[i, , drop=FALSE]) )[-58,] %>%
    gather( "gene", "ncount", -(1:5) )%>%
    ggplot +geom_line(aes(x=time, y=log2(ncount+1), col=samplelabel)) + facet_wrap(~gene)
  p<- p+ xlab("Time[h]") + ylab("Log2(count)+1")
  p<-p+theme_bw()+coord_fixed(ratio=.5)
  p<-p+geom_point(aes(x=time,y=log2(ncount+1),col=samplelabel))
  p<-p+ scale_x_continuous(breaks = c(0,2, 4, 6,8,10,12,14,16,18,20,22)) 
  p<-p+ ylim(0,16)
 p<-p+theme(axis.text = element_text (face="bold", size=14), axis.title=element_text(size=18,face="bold"))
 p<-p+theme(legend.title = element_text(size=16, face="bold"),legend.text = element_text(size = 16, face = "bold"))
 p<-p+theme(strip.text = element_text(size=25,face="bold"))
 #theme(legend.title = element_text(colour="blue", size=16, face="bold"))
ggsave(paste('~/Desktop/Light-dark/images/',i , '.png', sep=''), width=20, height=15, p)
  print(i)
}


########test html plot 
 HTMLoutput=file.path(".","output.html")
 graph1="light_inducibale_gene.png"
 HTMLInsertGraph(graph1,file=HTMLoutput,caption="Sample discrete distribution plot")

#######pdf(file='ressig_11hr.pdf', height=10, width=12)
xaxis=c(18,20,22,0,1,2,4,6,8,10,11,11.5,12,16)

pdf(file='sig_gene.pdf', height=10, width=12)
xaxis=c(0,1,2,4,6,8,10,11,11.5,12,16,18,20,22)

par(mfrow=c(3,3))


for(i in 1:nrow(idx))
{
    dd_wt=idx[i,c(38:40,29:37)]
    del_dd=idx[i,c(70:72,61:69)]
    #dd_wt=idx[i,29:40]
    #dd_pub=idx[i,41:52]
    #dd_df=idx[i,53:60]
    plot(xaxis,log2(dd_wt)+1,type="l",ylab="norm.count(log2+1)",xlab="Time(hours)",xlim=c(0,22),main=idx$ens.id[i],ylim=c(5,18))
    points(xaxis,log2(del_dd)+1,pch=18,col="blue",lwd = 5)
    lines(xaxis,log2(dd_wt)+1,col="blue")
    
    lines(xaxis,log2(del_df)+1,col="red")
    points(xaxis,log2(del_df)+1,pch=18,col="red",lwd = 5)
    
    
  #  x=c(0,2,4,6,8,10,12,14,16,18,20,22)
   # lines(x,log2(dd_wt)+1,col="green",,main=rownames(idx)[i])
    #points(x,log2(dd_wt)+1,pch=18,col="green",lwd = 5)
    
   # z=c(0,2,4,6,8,10,12,14,16,18,20,22)
    #lines(z,log2(dd_pub)+1,col="black",,main=rownames(idx)[i])
    #points(z,log2(dd_pub)+1,pch=18,col="black",lwd = 5)
    
    #t=c(0,2,4,6,8,10,12,14)
    #lines(t,log2(dd_df)+1,col="orange",,main=rownames(idx)[i])
    #points(t,log2(dd_df)+1,pch=18,col="orange",lwd = 5)
    
    legend("topright",c("dd-wt","led_dd"),fill=c("blue","red"))
}
dev.off()

####heatmap 

idx<-match(dd$ens.id,rownames(fullData))


heatmap.2(as.matrix( dd[,c(1:12)]), col=greenred(30), trace="none",
         Colv=FALSE, Rowv=FALSE, dendrogram = "none",scale="row",
           colsep=c(12)
           ,rowsep = c(98:165))







pdf(file='sig_gene.pdf', height=10, width=12)
#xaxis=c(0,1,2,4,6,8,10,11,11.5,12,16,18,20,22)

par(mfrow=c(3,3))
#x = seq( 0,22, by=2)
  xaxis=c(0,2,4,6,8,10,12,14,16,18,20,22)
#plot(x,y, xaxt='n')

#axis(1, at=x,labels=c(18 ,20 ,22 , 0  ,2  ,4  ,6 , 8 ,10 ,12, 14,16))
#del=idx[i,c(70:72,61:69)]
#dd_test=idx[c(38:40,29:37)]
for(i in 1:nrow(idx))
  {
jpeg(file=paste("/Users/amit/Desktop/Light-dark/rythim_gene_analysis/dark_dark_rep/image/", rownames(idx)[i],
                ".jpeg", sep=""),width = 4, height = 4, units = 'in', res = 300)
dd_wt=idx[i,c(38:40,29:37)]
csp1_dd=idx[i,c(70:72,61:69)]
  #gene=rownames(dd_wt)[i]
#dd_wt=idx[i,29:40]
  #dd_pub=idx[i,41:52]
  #dd_df=idx[i,53:60]
  plot(xaxis,dd_wt,type="l",ylab="count(log2+1)",xlab="Time(hr)",main=rownames(idx)[i],xaxt='n',ylim=c(0,22))
  points(xaxis,dd_wt,pch=18,col="blue",lwd = 1)
  lines(xaxis,dd_wt,col="blue")
  
lines(xaxis,csp1_dd,col="red")
  points(xaxis,csp1_dd,pch=18,col="red",lwd = 2)
  axis(1, at=x,labels=c(18,20,22,0,2,4,6, 8,10,12,14,16))
  
  #  x=c(0,2,4,6,8,10,12,14,16,18,20,22)
  # lines(x,log2(dd_wt)+1,col="green",,main=rownames(idx)[i])
  #points(x,log2(dd_wt)+1,pch=18,col="green",lwd = 5)
  
  # z=c(0,2,4,6,8,10,12,14,16,18,20,22)
  #lines(z,log2(dd_pub)+1,col="black",,main=rownames(idx)[i])
  #points(z,log2(dd_pub)+1,pch=18,col="black",lwd = 5)
  
  #t=c(0,2,4,6,8,10,12,14)
  #lines(t,log2(dd_df)+1,col="orange",,main=rownames(idx)[i])
  #points(t,log2(dd_df)+1,pch=18,col="orange",lwd = 5)
  
legend("topright",c("dd-wt","csp1_dd"),fill=c("blue","red"))

for(i in 1:nrow(idx))
 {
   jpeg(file=paste("~/Desktop/Light-dark/rythim_gene_analysis/dark_dark_rep/image/", rownames(idx)[i],
                   ".jpeg", sep=""), width = 3, height = 4, units = 'in', res = 300)
   dd_wt=idx[i,c(38:40,29:37)]
   csp1_dd=idx[i,c(70:72,61:69)]
   #gene=rownames(dd_wt)[i]
   #dd_wt=idx[i,29:40]
   #dd_pub=idx[i,41:52]
   #dd_df=idx[i,53:60]
   plot(xaxis,dd_wt,type="l",ylab="count(log2+1)",xlab="Time(hr)",main=rownames(idx)[i],xaxt='n',ylim=c(0,22))
   points(xaxis,dd_wt,pch=20,col="blue",lwd = 0.5)
   lines(xaxis,dd_wt,col="blue")
  
   lines(xaxis,csp1_dd,col="red")
   points(xaxis,csp1_dd,pch=20,col="red",lwd = 0.5)
   axis(1, at=xaxis,labels=c(18,20,22,0,2,4,6, 8,10,12,14,16))
   #  x=c(0,2,4,6,8,10,12,14,16,18,20,22)
# lines(x,log2(dd_wt)+1,col="green",,main=rownames(idx)[i])
  #points(x,log2(dd_wt)+1,pch=18,col="green",lwd = 5)  
  # z=c(0,2,4,6,8,10,12,14,16,18,20,22)
   #lines(z,log2(dd_pub)+1,col="black",,main=rownames(idx)[i])
   #points(z,log2(dd_pub)+1,pch=18,col="black",lwd = 5)
   
  #t=c(0,2,4,6,8,10,12,14)
 #lines(t,log2(dd_df)+1,col="orange",,main=rownames(idx)[i])
   #points(t,log2(dd_df)+1,pch=18,col="orange",lwd = 5)
   
   legend("topright",c("dd-wt","csp1_dd"),fill=c("blue","red"),x.intersp=0.5,y.intersp=1,cex=0.5, bty="n")
   
   #jpeg(file = paste(rownames(idx)[i], '.jpeg', sep = ''))
     dev.off()
  print(i)
 }

#jpeg(file = paste(rownames(idx)[i], '.jpeg', sep = ''))
dev.off()
print(i)
}


#paste("/Users/amit/Desktop/Light-dark/rythim_gene_analysis/dark_dark_rep/test/", rownames(idx)[i], ".jpeg", sep="")

#p<- data.frame(colData(fullData),t(counts(fullData, normalized=TRUE)[,]))[-58,]%>%
 # gather( "gene", "ncount", -(1:5) ) %>%
  #ggplot + geom_line(aes( x=time, y=log2(ncount+1), col=samplelabel)) +facet_wrap(~gene)
#p<- p + xlab( "Time[h]") + ylab("Log2(count)+1") + theme_bw()
#p<-p+geom_point(aes(x=time,y=log2(ncount+1),col=samplelabel))
#facet_multiple( plot = p, facets = 'gene', ncol = 4, nrow = 4)

#ggplot(cbind( as.data.frame(colData(fullData)), ncount=log2(counts(fullData,normalized=TRUE)["NCU00175",] )+1))+geom_line(aes(x=time,y=ncount,col=samplelabel))
#d<-ggplot(p+geom_line(aes(x=time,y=ncount,col=samplelabel)))




