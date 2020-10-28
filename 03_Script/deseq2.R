# #################################################################
# 
#                               DESeq2
#                          
# #################################################################

library(DESeq2)
library(readr)

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

#input_cts=snakemake@input[["cts"]]
project=snakemake@params[["project"]]
samples=snakemake@params[["samples"]]

# Rename column name of the count matrix as coldata
# colData and countData must have the same sample order
cts <- read.delim(snakemake@input[["cts"]], header=FALSE, comment.char="#", quote="",sep="\ ")

for (i in 1:length(project)) {
  cts[1,] <- lapply(cts[1,], sub, pattern = paste(project[[i]],"_",sep=""), replacement = "")
  cts[1,] <- lapply(cts[1,], sub, pattern = paste("_",samples[[i]],"\\.bam",sep=""), replacement = "")
}

# Format codData and countData for DESeq2
cts2 <- cts[,-1]
rownames(cts2) <- cts[,1]
cts <- cts2[-1,]
colnames(cts) <- cts2[1,]
cts <- data.matrix(cts)

coldata_read <- read.delim(snakemake@input[["coldata"]], header=TRUE, comment.char="#", quote="")
coldata <- coldata_read[,-1]
rownames(coldata) <- coldata_read[,1]
coldata$condition <- factor(coldata_read$condition)
coldata$type <- factor(coldata_read$type)

#write.table(coldata, file = snakemake@output[["count_modif"]], quote = FALSE, row.names = FALSE, col.names = FALSE)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 10, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)
# Add from the Lionel workflow
dds <- estimateSizeFactors( dds)

saveRDS(dds, file=snakemake@output[["rds"]])

#sessionInfo()