https://www.biostars.org/p/214100/

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, padj<.01 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(res, padj<.01 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.8))


pdf(paste("plots/", "volcano_", s1, "-", s2, ".pdf", sep=""))
t<-quantile(-log10(de.filtered$pval),0.999)
plot(de.filtered$log2Fold[(-log10(de.filtered$pval))<t[[1]]],-log10(de.filtered$pval[(-log10(de.filtered$pval))<t[[1]]]), pch=16, cex=0.6, col = ifelse((de.filtered$padj[(-log10(de.filtered$pval))<t[[1]]]<=0.05 & abs(de.filtered$log2Fold[(-log10(de.filtered$pval))<t[[1]]])>=1), alpha(color_de, 0.5), alpha(color_dots, 0.4)), ylab="-log10 pvalue", xlab="log2 fold change", main=paste("Comparison of samples ", s2, "vs", s1, sep=" "), ylim=c(0,(t[[1]])+5))

#Define special tickmarks at -1 and +1 on x-axis
axis (1, at = c(-1, 1, axis (1)))

points(de.filtered$log2Fold[-log10(de.filtered$pval)>=t[[1]]], rep(t[[1]]+1, length(de.filtered$pval[-log10(de.filtered$pval)>=t[[1]]])), pch=17, col=alpha(color_de, 0.5), cex=0.7)
abline(b=0, v=1, col=alpha ("black", 0.9))
abline(b=0, v=-1, col=alpha ("black", 0.9))
legend("topright", "Genes with a log2Fold change > 1 and padj <0.05", x.intersp=0.3, pch=20, col=alpha(color_de, 0.5), cex=0.6, horiz=TRUE)
dev.off()

