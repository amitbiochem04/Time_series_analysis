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
fit1 <- lm( y ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi ) )
fit0 <- lm( y ~ x )
amplitute[i]<- sqrt( coef(fit1)[3]^2 + coef(fit1)[4]^2 )
phase[i]<- atan2( coef(fit1)[4], coef( fit1)[3] )
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

