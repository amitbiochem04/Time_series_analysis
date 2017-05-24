
data<- log2( norm[,c(29:40)]+1 )
x<-colData(fullData)$time[29:40]
#x = seq( 0,, by=2)
approx <- as.data.frame( mat.or.vec( nrow( data), 16))
pval<- as.vector(c())
amplitute <- as.vector(c())
phase <- as.vector(c())
for( i in 1:nrow( data) ){
y <- as.numeric( data[ i,1:12 ])
fit1 <- glm.nb( y ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi ) )
fit0 <- glm.nb( y ~ x )
amplitute[i]<- sqrt( coef(fit1)[3]^2 + coef(fit1)[4]^2 )
#phase[i]<- atan2( coef(fit1)[4], coef( fit1)[3] )
model<-anova( fit0, fit1 )
p=-1*atan2( coef(fit1)[4], coef( fit1)[3] )*22/( 2*pi )
  	if (p<0) {
		p=-1*p
	} else {
		p=22-p
	}
phase[i]<-p
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



