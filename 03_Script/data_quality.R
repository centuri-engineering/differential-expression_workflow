## @knitr heatmap-count
# load deseq2 data
dds <- readRDS(RDS)
# this gives log2(n + 1)
ntd <- normTransform(dds)
# matrix of normalized values
# Variance stabilizing transformation
vsd <- vst(dds, blind=TRUE)
# Regularized log transformation
rld <- rlog(dds, blind=TRUE)

# 30 genes with top and bottom loadings
pcaobj <- prcomp(t(SummarizedExperiment::assay(rld)))
hi_load <- hi_loadings(pcaobj, topN = 50, exprTable = log2(counts(dds, normalized=TRUE)) + 1)

df <- as.data.frame(colData(dds)[,c("condition","type")])

pheatmap(hi_load, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main="log2 values")
#ggsave("heatmap.pdf",
#  plot = p,
#  width = 11, height = 8
#)

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
# PCA using the 50 genes with the top and bottom loadings calculated with the Regularized log transformation
#se <- SummarizedExperiment( hi_load, colData= colData( dds))
#plotPCA( DESeqTransform( se ), intgroup = "condition") + aes(label=name)

# PCA using all gene value in rlog with different component
df <- cbind(df, pcaobj$x)
ggplot(df) + geom_point(aes(x=PC1, y=PC2, color=condition))
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color=condition))

## @knitr plot-gene
gene_name=snakemake@params[["gene_name"]]

d <- plotCounts(dds, gene=gene_name, intgroup="condition", returnData=TRUE)
# Plotting the gene normalized counts, using the samplenames (rownames of d as labels)
ggplot(d, aes(x = condition, y = count, color = condition)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle(gene_name) +
  theme(plot.title = element_text(hjust = 0.5))
