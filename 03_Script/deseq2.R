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

# Loading the parameters
project=snakemake@params[["project"]]
samples=snakemake@params[["samples"]]
ref_level=snakemake@params[["ref_level"]]
normalized_counts_file=snakemake@output[["normalized_counts_file"]]

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

# Create the DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

# Remove uninformative columns (to do when filter not already done with the CPM threshold)
#dds <- dds[ rowSums(counts(dds)) > 10, ]

# Specifying the reference level
dds$condition <- relevel(dds$condition, ref = ref_level)

# DESeq : Normalization and preprocessing (counts divided by sample-specific size factors
# determined by median ratio of gene counts relative to geometric mean per gene)
dds <- DESeq(dds, parallel=parallel)
# To save the object in a file for later use
saveRDS(dds, file=snakemake@output[["rds"]])

# Already done in the DESeq function
#dds <- estimateSizeFactors( dds)

# Save the normalized data matrix
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file=normalized_counts_file, sep="\t", quote=F, col.names=NA)
