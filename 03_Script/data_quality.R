
## @knitr plot-pca
library("DESeq2")
library("plotly")

# load deseq2 data
dds <- readRDS("../05_Output/08_deseq2_init/all.rds")
#dds <- readRDS(RDS)
# this gives log2(n + 1)
ntd <- normTransform(dds)

# obtain normalized counts
counts <- rlog(dds, blind=FALSE)
pdf(file.path( OUTPUT_DIR, STEP, "pca.pdf"))
plot = plotPCA(ntd)
ggplotly(plot)
dev.off()
