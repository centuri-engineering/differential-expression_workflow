## @knitr heatmap-count

# load deseq2 data
dds <- readRDS(RDS)
# this gives log2(n + 1)
ntd <- normTransform(dds)
# matrix of normalized values
# Variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)
# Regularized log transformation
rld <- rlog(dds, blind=FALSE)

# 30 genes with top and bottom loadings
pcaobj <- prcomp(t(SummarizedExperiment::assay(rld)))
hi_load <- hi_loadings(pcaobj, topN = 30, exprTable = log2(counts(dds, normalized=TRUE)) + 1)

df <- as.data.frame(colData(dds)[,c("condition","type")])

pheatmap(hi_load, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main="log2 values")

## @knitr heatmap-distance
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rownames(colData( dds))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

## @knitr plot-pca
name = rownames(colData( dds))
# PCA using the 30 genes with the top and bottom loadings calculated with the Regularized log transformation
se <- SummarizedExperiment( hi_load, colData= colData( dds))
pca_and_plot = plotPCA( DESeqTransform( se ), intgroup = "condition") + aes(label=name)
ggplotly(pca_and_plot)
