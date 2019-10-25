###################
# batch correction
###################

library(edgeR)
library(DESeq2)

plot_batch <- function(cpm, outfpath = NULL) {
  sampleDistMatrix <- as.matrix(dist(t(cpm)))
  rownames(sampleDistMatrix) <- paste(colnames(cpm), sep='-")
  if (!is.null(outfpath)) {
    svg(outfpath)
    pheatmap(sampleDistMatrix,
             cluster_col = FALSE,
             cluster_row = FALSE)
    dev.off()
  } else {
    pheatmap(sampleDistMatrix,
             cluster_col = FALSE,
             cluster_row = FALSE)
  }
}

batch_correct <- function(meta.fpath, out_dir = NULL) {
  meta <- read.delim(meta.fpath, header = T, sep = "\t")
  sampleTable <- data.frame(sampleName = meta$condition,
                            fileName   = meta$fileName,
                            condition  = meta$sampletype,
                            time       = meta$time,
                            sampletype = meta$sampletype,
                            samplelabel = meta$samplelabel)
  de <- DEeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = dirpath,
                                  design = ~ 1)
  norm <- counts(de)
  dark <- norm[, c(29:52)]
  y  <- DGEList(counts=dark)
  ba <- factor(c(rep(1, times=12),
                 rep(0, times=12)))
  A  <- aveLogCPM(y)
  y2 <- y[A>1, ]
  y2 <- calcNormFactors(y2)
  CPM <- cpm(y2, log=TRUE, prior.count=5)
  CPMc <- removeBatchEffect(CPM, batch)
  if (is.null(out_dir)) {
    out_dir <- dirname(dirpath)
  }
  fname <- basename(meta.fpath)
  out.fpath <- file.path(out_dir, paste0(fname, ".before.batch.correction.svg"))
  plot_batch(CPM, outfpath = out.fpath)
  out.fpath <- file.path(out_dir, paste0(fname, ".after.batch.correction.svg"))
  plot_batch(CPMc, outfpath = out.fpath)
}

batch_correct(meta.fpath = "/Users/amitsingh/Desktop/Light-dark/meta_1.txt")
