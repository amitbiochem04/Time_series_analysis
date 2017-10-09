####norm data of dl-wt 
data<- log2( norm[,c(29:40)]+1 )
x<-colData(fullData)$time[29:40]
#x = seq( 0,, by=2)
approx <- as.data.frame( mat.or.vec( nrow( data), 16))
pval<- as.vector(c())
amplitute <- as.vector(c())
phase <- as.vector(c())
for( i in 1:nrow( data) ){
y <- as.numeric( data[ i,1:12 ])
fit1 <- lm( y ~ x + cos( x/24*2*pi ) + sin( x/24*2*pi ) )
fit0 <- lm( y ~ x )
amplitute[i]<- sqrt( coef(fit1)[3]^2 + coef(fit1)[4]^2 )
phase[i]<- -atan2( coef(fit1)[4], coef( fit1)[3] )
model<-anova( fit0, fit1 )
pval[i]<-model$`Pr(>F)`[2]
approx <- as.data.frame( cbind( data,pval,amplitute,phase) )
approx$p.adj<-p.adjust( approx$pval, method = "BH",n = length( approx$pval))
print(i)
}
###
sig<- subset( approx,p.adj< 0.1)
sig<-sig[ order( sig[,"p.adj"] ), ]

keep<-list()
for (i in rownames( sig)){
  idx<-grep(i,anno$X.Gene.ID.)
  idx<-anno[ idx,]
  keep[[i]]<-idx
} 
df<-melt(keep,id=c("X.Gene.ID.","X.Gene.Name.or.Symbol."))
df<-df[,1:2]
colnames(df)<-c("ens.id","symbol")
#df<-df[!duplicated(df), ]
sig<-as.data.frame( cbind( sig,df[!duplicated(df), ]))



####plot
dd<-subset(sig,p.adj < 0.1)
dd<-dd[ order( dd[,"phase"] ), ]
#dd$phase<-abs(dd$phase)
mydata<-dd[,c(1:11)]
category=c(rep("Day_gene", 50),rep("Morning_and_evening_gene", 14), rep("Morning_gene",16 ))
mydf<-data.frame(row.names(mydata),category)
CairoPDF(file="Rythmic_gene_with 0h p.adj < 0.1",width=10/1.54,height=10/1.54)
#heatmap.2(as.matrix( dd[,c(1:11)]),Colv=F,col=greenred(30),trace="none",scale="row",dendrogram = "none")
heatmap.2(as.matrix(mydata), col=greenred(30), trace="none",
          Colv=FALSE, Rowv=FALSE, dendrogram = "none",scale="row",rowsep = c(50,64),
          keysize=0.5,labRow=NA,density.info="none",
        RowSideColors=col1[as.numeric(mydf$category)] )
#legend(x="center", legend=c("genes","genes2","Morning"),col=col1[as.numeric(mydf$category)], pch=15)
#legend(x="topleft", legend=c("genes","genes2","Morning"), col=c("red","blue","black"), pch=15)

#legend("topright",      
 #      legend = category,
  #     col = col1[as.numeric(mydf$category)], lty= 1,             
   #    lwd = 5,cex=.7)











####
sig<-sig[,c(25:30)]
page <- openPage( "ryth_gene_final.html",
                  head = paste( sep="\n",
                                "<script>",
                                "   function set_image( name ) {",
                                "      document.getElementById( 'plot' ).setAttribute( 'src', 'images/' + name + '.png' );",
                                "   }",
                                "</script>" ) )
cat(file=page,
    '<table><tr><td style="vertical-align:top"><div style="height:800px; overflow-y:scroll">' )
#####
bgcolor1=matrix( c( '#FAEBD7'), nr=10, nc=6 )
bgcolor2=matrix( c( '#E0EEEE'), nr=793, nc=6 )
hwrite(sig, border=NULL, page=page, bgcolor=rbind( bgcolor1,bgcolor2),
onmouseover = sprintf( "set_image( '%s' );",sig$ens.id))
cat( file=page,
     '</div></td><td style="vertical-align:top"><img id="plot" width="700px"></td></tr></table>' )
closePage(page)
browseURL( "ryth_gene_final.html" )

######
dd <- fullData[,fullData$samplelabel=="dd-wt"]
dd$phase <- as.integer( dd$time) / 22 * 2*pi
#colData( dd)$condition<-as.factor( c( "0h","2h","4h", "6h","8h","10h","12h","14h","16h","18h","20h","22h") )
colData( dd)$samplelabel<-as.factor(( dd)$samplelabel)
design = model.matrix( ~ dd$condition*( dd$time + cos( dd$phase) + sin( dd$phase) ))
#design = model.matrix( ~ dd$samplelabel*( dd$time + cos( dd$phase) + sin( dd$phase) ))


#y<- DGEList( counts = assay( dd))

#design <- ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi )
#table(rowSums(y$counts==0)==9)

keep <- rowSums( cpm( y) > 0.5) >= 2
table( keep)
#keep
y<-y[keep, ,keep.lib.sizes=FALSE]
y <- calcNormFactors( y)
v <- voom( y, design, plot=TRUE)
efit <- eBayes( vfit)

#########final 
mat <- norm[,c(29:52)]
#####take avrergae 
newmat <- as.data.frame( mat.or.vec( nrow( mat), 12))
newmat[,1] = apply(mat[,c(1,13)],1,mean)
newmat[,2] = apply(mat[,c(2,14)],1,mean)
newmat[,3] = apply(mat[,c(3,15)],1,mean)
newmat[,4] = apply(mat[,c(4,16)],1,mean)
newmat[,5] = apply(mat[,c(5,17)],1,mean)
newmat[,6] = apply(mat[,c(6,18)],1,mean)
newmat[,7] = apply(mat[,c(7,19)],1,mean)
newmat[,8] = apply(mat[,c(8,20)],1,mean)
newmat[,9] = apply(mat[,c(9,21)],1,mean)
newmat[,10] = apply(mat[,c(10,22)],1,mean)
newmat[,11] = apply(mat[,c(11,23)],1,mean)
newmat[,12] = apply(mat[,c(12,24)],1,mean)
rownames(newmat)<-rownames(mat)
colnames(newmat)<-c("0hr","2hr","4hr","6hr","8hr","10hr","12hr","14hr","16hr","18hr","20hr","22hr")
save(newmat,file="dark_rep_mean.rda")
newmat<-log2(newmat[,c(1:12)]+1)

x<-colData(fullData)$time[c(30:40)]
#data<-newmat[,c(1:12)]
##
approx <- as.data.frame( mat.or.vec( nrow( newmat), 15))
pval<- as.vector( c() )
amplitute <- as.vector( c() )
phase <- as.vector( c() )
 
 for( i in 1:nrow( newmat ))  {
  y <- as.numeric( newmat[i,2:12 ] )
  fit1 <- lm( y ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi ) )
  fit0 <- lm( y ~ x )
  #fit1<-glm.nb( y ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi ) )
  #fit0 <-glm.nb( y ~ x )
  amplitute[i]<- sqrt( coef(fit1 )[3]^2 + coef( fit1 )[4]^2 )
  #phase[i]<- atan2( coef(fit1)[4], coef( fit1)[3] )
  phase[i]=-atan2( coef( fit1 )[4], coef( fit1 )[3] )*22/( 2*pi )
  model<-anova( fit0, fit1 )
  #model<-anova( fit0, fit1, test="Chisq" )
  #pval[i]<-model$`Pr(>F)`[2]
  pval[i]<-model$`Pr(>F)`[2]
  #pval[i]<-model[["Pr(Chi)"]][2]

  print(i)
 }

approx <- as.data.frame( cbind(newmat, pval, amplitute, phase ) )
#approx$p.adj<-p.adjust( pval, method = "BH",n = length( pval ) )
approx$p.adj<-p.adjust( approx$pval, method = "BH",n = length( approx$pval ) )
#approx$p.adj_homel<-p.adjust( approx$pval, method = "hommel",n = length( approx$pval ) )
#approx$p.adj_homel<-p.adjust( approx$pval, method = "bonferroni",n = length( approx$pval ) )

#######
#sig<- subset(approx,p.adj< 0.1)
sig<-approx[ order(approx[,"p.adj"] ), ]
keep<-list()
for (i in rownames( sig)){
  idx<-grep(i,anno$X.Gene.ID.)
  idx<-anno[ idx,]
  keep[[i]]<-idx
} 
df<-melt(keep,id=c("X.Gene.ID.","X.Gene.Name.or.Symbol."))
df<-df[,1:2]
colnames(df)<-c("ens.id","symbol")
df<-df[!duplicated(df), ]
sig<-as.data.frame( cbind( sig,df[!duplicated(df), ]))
save(sig,file="ryth_dd_rep_with_out_0_log.rda")

write.xlsx(sig, "rythmicgene_with_out0.xlsx",row.names=TRUE,asTable = TRUE)
####
dd<-subset(sig,p.adj < 0.1)
dd<-dd[ order( dd[,"phase"] ), ]
#dd$phase<-abs(dd$phase)
mydata<-dd[,c(1:11)]
category=c(rep("Day_gene", 50),rep("Morning_and_evening_gene", 14), rep("Morning_gene",16 ))
mydf<-data.frame(row.names(mydata),category)
col1<-c("red","blue","black")
CairoPDF(file="Rythmic_gene_with out 0h",width=10/1.54,height=10/1.54)
#heatmap.2(as.matrix( dd[,c(1:11)]),Colv=F,col=greenred(30),trace="none",scale="row",dendrogram = "none")
heatmap.2(as.matrix(mydata), col=greenred(30), trace="none",
          Colv=FALSE, Rowv=FALSE, dendrogram = "none",scale="row",rowsep = c(50,64),
         labRow=NA,density.info="none",
        RowSideColors=col1[as.numeric(mydf$category)],breaks = seq(-2, 2,length.out = 30+1) )
dev.off()

######


x<-colData(fullData)$time[c(29:40)]
#data<-newmat[,c(1:12)]
##
approx <- as.data.frame( mat.or.vec( nrow( newmat), 15))
pval<- as.vector( c() )
amplitute <- as.vector( c() )
phase <- as.vector( c() )
 
 for( i in 1:nrow( newmat ))  {
  y <- as.numeric( newmat[i,1:12 ] )
  fit1 <- lm( y ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi ) )
  fit0 <- lm( y ~ x )
  #fit1<-glm.nb( y ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi ) )
  #fit0 <-glm.nb( y ~ x )
  amplitute[i]<- sqrt( coef(fit1 )[3]^2 + coef( fit1 )[4]^2 )
  #phase[i]<- atan2( coef(fit1)[4], coef( fit1)[3] )
  phase[i]=-atan2( coef( fit1 )[4], coef( fit1 )[3] )*22/( 2*pi )
  model<-anova( fit0, fit1 )
  #model<-anova( fit0, fit1, test="Chisq" )
  #pval[i]<-model$`Pr(>F)`[2]
  pval[i]<-model$`Pr(>F)`[2]
  #pval[i]<-model[["Pr(Chi)"]][2]

  print(i)
 }

approx <- as.data.frame( cbind(newmat, pval, amplitute, phase ) )
#approx$p.adj<-p.adjust( pval, method = "BH",n = length( pval ) )
approx$p.adj<-p.adjust( approx$pval, method = "BH",n = length( approx$pval ) )
#approx$p.adj_homel<-p.adjust( approx$pval, method = "hommel",n = length( approx$pval ) )
#approx$p.adj_homel<-p.adjust( approx$pval, method = "bonferroni",n = length( approx$pval ) )

#######
#sig<- subset(approx,p.adj< 0.1)
sig<-approx[ order(approx[,"p.adj"] ), ]
keep<-list()
for (i in rownames( sig)){
  idx<-grep(i,anno$X.Gene.ID.)
  idx<-anno[ idx,]
  keep[[i]]<-idx
} 
df<-melt(keep,id=c("X.Gene.ID.","X.Gene.Name.or.Symbol."))
df<-df[,1:2]
colnames(df)<-c("ens.id","symbol")
df<-df[!duplicated(df), ]
sig<-as.data.frame( cbind( sig,df[!duplicated(df), ]))
save(sig,file="ryth_dd_rep_with_0_log.rda")
write.xlsx(sig, "rythmicgene.xlsx",row.names=TRUE,asTable = TRUE)
#####plot
dd<-subset(sig,p.adj < 0.1)
#dd<-dd[ order( dd[,"phase"] ), ]
dd<-dd[ order( dd$phase , decreasing=FALSE),]
#dd$phase<-abs(dd$phase)
mydata<-dd[,c(1:12)]
category=c(rep("Day_gene", 91),rep("Morning_and_evening_gene", 88), rep("Morning_gene",38 ))
mydf<-data.frame(row.names(mydata),category)
col1<-c("red","blue","black")
CairoPDF(file="Rythmic_gene_with_0h",width=10/1.54,height=10/1.54)
#heatmap.2(as.matrix( dd[,c(1:11)]),Colv=F,col=greenred(30),trace="none",scale="row",dendrogram = "none")
heatmap.2(as.matrix(mydata), col=greenred(30), trace="none",
          Colv=FALSE, Rowv=FALSE, dendrogram = "none",scale="row",rowsep = c(91,179),
         labRow=NA,density.info="none",
        RowSideColors=col1[as.numeric(mydf$category)],breaks = seq(-2, 2,length.out = 30+1))
dev.off()

###########



#######
index=(apply(newmat,1,sum)>20) 

row_sub = apply(newmat, 1, function(row) all(row !=0 ))
data<-newmat[index,]

approx <- as.data.frame( mat.or.vec( nrow(data), 15))
pval<- as.vector( c() )
amplitute <- as.vector( c() )
phase <- as.vector( c() )
x<-colData(fullData)$time[c(29:40)]
for( i in 1:nrow( data ))  {
  y <- as.numeric( data[i,1:12 ] )
  #fit1 <- lm( y ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi ) )
  #fit0 <- lm( y ~ x )
  fit1<-glm.nb( y ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi ) )
  fit0 <-glm.nb( y ~ x )
  amplitute[i]<- sqrt( coef(fit1 )[3]^2 + coef( fit1 )[4]^2 )
  #phase[i]<- -atan2( coef(fit1)[4], coef( fit1)[3] )
  phase[i]=atan2( coef( fit1 )[4], coef( fit1 )[3] )*22/( 2*pi )
  model<-anova( fit0, fit1 )
  model<-anova( fit0, fit1, test="Chisq" )
  #pval[i]<-model$`Pr(>F)`[2]
  pval[i]<-model$`Pr(>F)`[2]
  #pval[i]<-model[["Pr(Chi)"]][2]
#  approx <- as.data.frame( cbind( data, pval, amplitute, phase ) )
 # approx$p.adj<-p.adjust( approx$pval, method = "BH",n = length( approx$pval ) )
  #approx$p.adj_by<-p.adjust( approx$pval, method = "BY",n = length( approx$pval ) )
  #approx$p.adj_homel<-p.adjust( approx$pval, method = "hommel",n = length( approx$pval ) )
  #approx$p.adj_homel<-p.adjust( approx$pval, method = "bonferroni",n = length( approx$pval ) )
  print(i)
}

#######
sig<- subset(approx,p.adj< 0.1)
sig<-approx[ order(approx[,"p.adj"] ), ]
keep<-list()
for (i in rownames( sig)){
  idx<-grep(i,anno$X.Gene.ID.)
  idx<-anno[ idx,]
  keep[[i]]<-idx
} 
df<-melt(keep,id=c("X.Gene.ID.","X.Gene.Name.or.Symbol."))
df<-df[,1:2]
colnames(df)<-c("ens.id","symbol")
#df<-df[!duplicated(df), ]
sig<-as.data.frame( cbind( sig,df[!duplicated(df), ]))
save(sig,file="ryth_dd_rep_with__out_0_log.rda")
                
########################
                
                y1<-c(NA,NA,NA,8,8,9,7,4,8);
x<-c("fail","fail","pass","fail","pass","fail","fail","fail","pass")
y2<-c(NA,8.14,8,9.26,8,10,6.13,3.10,9.13);
lib<-c(0,1,1,1,0,0,0,1,1)
##dat2$scode[dat2$sex=="M"]<-"1" 
#dat2$scode[dat2$sex=="F"]<-"0" 
 sex<-c(0,1,1,1,0,0,0,1,0)
 model.omit <- lm(y2 ~ y1, na.action = na.omit);
 resid(model.omit)
 ##here you see that there model fit from 4 to 9 observation 
 ##but you need to use following function which ll not take na value in account
 #and will ignore the value
 model.exclude <- lm(y2 ~ y1+lib+sex,na.action = na.exclude)
 model.exclude2 <- lm(y2 ~ y1,na.action = na.exclude)
 final<-anova( model.exclude, model.exclude2)
 ##
#resid( model.exclude)
 #logLik(model.exclude)
 #summary(model.exclude)$r.squared
 #anova(model.exclude)
 ### 
 #coef(summary(model.exclude))[, "Pr(>|t|)"]
 ##
#plot(model.exclude)
#summary(model.exclude)$coefficients[2,4]
#abline(model.exclude$coef[1],model.exclude$coef[2],model.exclude$coef[3],model.exclude$coef[4])



#####

###here is the real model 
x<-dss
#x = seq( 0,, by=2)
approx <- as.data.frame( mat.or.vec( nrow( data), 52))
model_pval<- as.vector(c())
sex<- as.vector(c())
batch<- as.vector(c())
Log2FC<-as.vector(c())
FC<-as.vector(c())
for( i in 1:nrow( data) ){
  y <- as.numeric( data[ i,1:52 ])
  ###here is the full model 
  fit1 <- lm( y ~ x + Extraction_Method+colData(rnaSE1), na.action = na.exclude)
  ###here is the minmal model
  fit0 <- lm( y ~ x, na.action = na.exclude)
###lets look for limma fold chnage how that calulate 
##here i took predicted value 
B<-mean(predict(fit1),na.rm=TRUE)
A<-mean(i,na.rm=TRUE)
Log2FC[i] = log2(B) - log2(A)
FC[i] = 2 ^ log2FC
###to get final pvlaue we can use anova ##
###here you should comapre two model by anovan that would be the final model information
model<-anova( fit0, fit1 )
###get the p value for the model, this pvalue we are actually looking
model_pval[i]<-model$`Pr(>F)`[2]
##if you like to comapre or tell some words about the age has any effect or not 
##then you have to have a pvalue for that so we follwoing way.
batch_pval[i]<-summary( fit )$coefficients[ 3,4 ]
sex_pval[i]<-summary( fit )$coefficients[4,4]
####just put all the data in a same data file so that you can also check for expression wise 
approx <- as.data.frame( cbind( data,Log2FC,FC,model_pval, batch_pval, sex_pval) )
###now adjust your pvalue
approx$model_p.adj<-p.adjust( approx$model_pval, method = "BH",n = length( approx$model_pval))
approx$batch_pval<-p.adjust( approx$batch_pval, method= "BH",n = length( approx$batch_pval))
approx$sex_pval<-p.adjust( approx$sex_pval, method = "BH",n = length( approx$sex_pval))
#### we can also try any other p.val adjust ment just to check if you like 
##duniya mein jitna correction method hoga utna lele
#approx$model_p.adj_homel<-p.adjust( approx$pval, method = "hommel",n = length( approx$pval ) )
#approx$$model_p.adj_homel<-p.adjust( approx$pval, method = "bonferroni",n = length( approx$pval ) )
print(i)
}
###select significant gene from the model that we fitted now 
##select p.adj as your interest 

final_gene<-subset(approx,model_p.adj<0.1)

###now you can annotate your gene that you select and get entrezid, put in edgre R 
library(edgeR)
tr<-rownames(final_gene)
###get the entrez gene from ensebmle
rownames(final_gene) <- as.character(mget(rownames(final_gene),org.Hs.egENSEMBL2EG,ifnotfound=NA))
keg <- kegga(tr, species="Hs")
###we need only up regulated pathways so we write "up", n=10 we need top 10
final<- topKEGG(keg, n=10, truncate=34, sort="up")
##plot them based on pvale
write.csv(final,file="keegpathway.csv")
#####plot in heatmap 
##put your code according to your column
a<-read.delim("keegpathway.csv",sep=",",header=T)
p<-ggplot(data=,aes(reorder(x=X,-log10(P.Up)),y=-log10(P.Up )))+ coord_flip() +geom_bar(colour="black", stat="identity",width=.3)+theme_minimal()+ylab(expression(paste('-log10',italic(q),'value')))+ggtitle("Up_GOTerm EvsL")
p <- p + xlab("Keegpathway")+theme(axis.text = element_text(size =13),axis.title=element_text(size=13))+theme(axis.ticks = element_line(size = 2))
p<-p+theme(axis.title.x = element_text(face="bold",size=13))
p
dev.off()





##extact all stat value 

#library(broom)
##to check summary 
#tidy(fit1)
###to check f-stat
#glance(fit1)
###to get p-vlaue
#glance(fit1)$p.value

https://www.nature.com/ng/journal/v32/n4s/pdf/ng1032.pdf
