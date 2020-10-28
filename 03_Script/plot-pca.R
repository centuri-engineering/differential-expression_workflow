
## @knitr plot-pca
library("DESeq2")

# load deseq2 data
dds <- readRDS(snakemake@input[["rds"]])

# this gives log2(n + 1)
ntd <- normTransform(dds)

# obtain normalized counts
counts <- rlog(dds, blind=FALSE)
#pdf(snakemake@output[[1]])
plotPCA(ntd, intgroup=snakemake@params[["pca_labels"]])
#dev.off()