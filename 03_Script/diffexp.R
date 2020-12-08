## @knitr diffexp

ref_level=snakemake@params[["ref_level"]]
lfcshrink_type=snakemake@params[["lfcshrink_type"]]

## Ploduce the MA plot and volcanoplot with matrices 

cat("We use ",lfcshrink_type ," lfcshrink type as adaptive t prior shrinkage estimator for ranking and visualization.\n")
cat("The Sample of reference for the pairwise comparison is ",ref_level,".\n")

# Automation of the comparisons with the reference condition established in (config.yaml for deseq2.R) 
comparison <- resultsNames(dds)
comparison_df <- as.data.frame(comparison)[-1,]
for (i in 1:length(comparison_df)) {
    outputname <- as.character(lapply(comparison_df[[i]], sub, pattern = "condition_", replacement = ""))
    condition <- strsplit(comparison_df[[i]], split = "_vs_")
    condition <- as.character(lapply(condition[[1]], sub, pattern = "condition_", replacement = ""))
    results <- results(dds, contrast=c("condition",condition))
    # Log fold change shrinkage for visualization and ranking (shrink fold changes for lowly expressed genes)
    rlog_results <- lfcShrink(dds, coef = comparison_df[[i]], res=results, type=lfcshrink_type)
    # p-values and adjusted p-values
    rlog_results <- rlog_results[order(rlog_results$padj),]
    write.table(as.data.frame(rlog_results), file= paste("../05_Output/09_differential_expression/",outputname,"_all_genes_stats.tsv", sep=""), quote = FALSE, sep = "\t")

    ## MA-plot
    #pdf(file = paste("../05_Output/Figue/",outputname,"_maplot.pdf"), width = 8, height = 9, pointsize = 10)
    xlim <- c(500,5000); ylim <- c(-2,2)
    plotMA(rlog_results, xlim=xlim, ylim=ylim, main=comparison_df[[i]])
    idx <- identify(rlog_results$baseMean, rlog_results$log2FoldChange)
    #dev.off()

    ## Volcano plot
    FC <- 0.584963
    p <- 10e-3

    keyvals <- rep('grey75', nrow(rlog_results))
    names(keyvals) <- rep('NS', nrow(rlog_results))

    keyvals[which(abs(rlog_results$log2FoldChange) > FC & rlog_results$padj > p)] <- 'grey50'
    names(keyvals)[which(abs(rlog_results$log2FoldChange) > FC & rlog_results$padj > p)] <- 'log2FoldChange'

    keyvals[which(abs(rlog_results$log2FoldChange) < FC & rlog_results$padj < p)] <- 'grey25'
    names(keyvals)[which(abs(rlog_results$log2FoldChange)  < FC & rlog_results$padj < p)] <- '-Log10Q'

    keyvals[which(rlog_results$log2FoldChange < -FC & rlog_results$padj < p)] <- 'blue2'
    names(keyvals)[which(rlog_results$log2FoldChange  < -FC & rlog_results$padj < p)] <- 'Signif. down-regulated'
    sdr <- subset(rlog_results, (rlog_results$log2FoldChange  < -FC) & (rlog_results$padj < p))
    write.table(as.data.frame(sdr), file= paste("../05_Output/09_differential_expression/",outputname,"_signif-down-regulated.txt", sep=""), quote = FALSE, sep = "\t")

    keyvals[which(rlog_results$log2FoldChange > FC & rlog_results$padj < p)] <- 'red2'
    names(keyvals)[which(rlog_results$log2FoldChange > FC & rlog_results$padj < p)] <- 'Signif. up-regulated'
    sur <- subset(rlog_results, (rlog_results$log2FoldChange > FC) & (rlog_results$padj < p))
    write.table(as.data.frame(sur), file= paste("../05_Output/09_differential_expression/",outputname,"_signif-up-regulated.txt", sep=""), quote = FALSE, sep = "\t")

    unique(keyvals)
    unique(names(keyvals))

    ## Volcanoplot with the adjusted p-value
    volcanoplot_padj <- EnhancedVolcano(rlog_results,
    lab = rownames(rlog_results),
    x = 'log2FoldChange',
    y = "padj",
    xlim = c(-6,6),
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~adjusted~italic(P)),
    title = "",
    subtitle = "",
    axisLabSize = 12,
    titleLabSize = 15,
    pCutoff = 10e-3,
    FCcutoff = 0.6,
    pointSize = 1.5,
    labSize = 2,
    colCustom = keyvals,
    colAlpha = 1,
    legendPosition = 'bottom',
    legendLabSize = 10,
    legendIconSize = 1.5,
    #DrawConnectors = TRUE,
    #widthConnectors = 0.2,
    colConnectors = 'grey50',
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    border = 'full',
    borderWidth = 0.7,
    borderColour = 'black')    
    plot(volcanoplot_padj)
 }

## @knitr toppadj

mutant_level=snakemake@params[["mutant_level"]]
stat_file=paste("../05_Output/09_differential_expression/",mutant_level,"_vs_",ref_level,"_all_genes_stats.tsv", sep="")
nbpval=snakemake@params[["nbpval"]]
mutant_level_file <- read_delim(stat_file, "\t", escape_double = FALSE, trim_ws = TRUE)

# Top adjusted p-value
top_adjpval <- mutant_level_file[1:nbpval,1]
# Keep the normalized count values of these genes
dds <- readRDS(RDS)
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
normalized_counts <- rownames_to_column(normalized_counts, "Gene")
topnorm <- data.frame()
for (i in 1:nrow(top_adjpval)){
    genename=as.character(top_adjpval[i,1])
    ligne=subset(normalized_counts,normalized_counts$Gene==genename)
    topnorm <- rbind(topnorm,ligne)
}
topnorm_length=length(topnorm)
gathered_topnorm <- topnorm %>% gather(colnames(topnorm)[2:topnorm_length], key = "samplename", value = "normalized_counts")

# Counts colored by sample group
coldata_read <- read.delim(paste("../",snakemake@input[["coldata"]],sep=""), header=TRUE, comment.char="#", quote="")
coldata <- coldata_read[,-3]
names(coldata)[names(coldata) == "project"] <- "samplename"
gathered_topnorm <- inner_join(coldata, gathered_topnorm)

cat(paste("\nWe plot the normalized count values for the top ",nbpval," differentially expressed genes (by padj values).\n",sep=""))

cat(paste("\nThese genes are selected with the padj-value obtened from the comparison ",mutant_level," vs ", ref_level,"\n", sep=""))

p=ggplot(gathered_topnorm) +
        geom_point(aes(x = Gene, y = normalized_counts, color = condition)) +
        scale_y_log10() +
        xlab("Genes") +
        ylab("Normalized Counts") +
        ggtitle(paste("Top ", nbpval, " Significant DE Genes",sep="")) +
        theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	theme(plot.title = element_text(hjust = 0.5))
ggplotly(p)
