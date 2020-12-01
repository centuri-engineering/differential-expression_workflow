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
    write.table(as.data.frame(rlog_results), file= paste("../05_Output/09_differential_expression/",outputname,".csv", sep=""))

    ## MA-plot
    #pdf(file = paste(outputname,"maplot.pdf"), width = 9, height = 9, pointsize = 10)
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
    write.table(as.data.frame(sdr), file= paste("../05_Output/09_differential_expression/",outputname,"_signif-down-regulated.csv", sep=""))

    keyvals[which(rlog_results$log2FoldChange > FC & rlog_results$padj < p)] <- 'red2'
    names(keyvals)[which(rlog_results$log2FoldChange > FC & rlog_results$padj < p)] <- 'Signif. up-regulated'
    sur <- subset(rlog_results, (rlog_results$log2FoldChange > FC) & (rlog_results$padj < p))
    write.table(as.data.frame(sur), file= paste("../05_Output/09_differential_expression/",outputname,"_signif-up-regulated.csv", sep=""))

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
    title = "Volcano plot with the adjusted p-values",
    subtitle = "",
    titleLabSize = 15,
    pCutoff = 10e-3,
    FCcutoff = 0.6,
    pointSize = 1.5,
    labSize = 2.5,
    colCustom = keyvals,
    colAlpha = 1,
    legendPosition = 'bottom',
    legendLabSize = 10,
    legendIconSize = 3.0,
    #DrawConnectors = TRUE,
    widthConnectors = 0.2,
    colConnectors = 'grey50',
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    border = 'partial',
    borderWidth = 1.5,
    borderColour = 'black')
    #pdf(paste(outputname,"volcano.pdf"), width = 9, height = 9, pointsize = 10)
    plot(volcanoplot_padj)
    #dev.off()

    ## Volcanoplot with the p-value
    # keyvals <- rep('grey75', nrow(rlog_results))
    # names(keyvals) <- rep('NS', nrow(rlog_results))

    # keyvals[which(abs(rlog_results$log2FoldChange) > FC & rlog_results$pvalue > p)] <- 'grey50'
    # names(keyvals)[which(abs(rlog_results$log2FoldChange) > FC & rlog_results$pvalue > p)] <- 'log2FoldChange'

    # keyvals[which(abs(rlog_results$log2FoldChange) < FC & rlog_results$pvalue < p)] <- 'grey25'
    # names(keyvals)[which(abs(rlog_results$log2FoldChange)  < FC & rlog_results$pvalue < p)] <- '-Log10Q'

    # keyvals[which(rlog_results$log2FoldChange < -FC & rlog_results$pvalue < p)] <- 'blue2'
    # names(keyvals)[which(rlog_results$log2FoldChange  < -FC & rlog_results$pvalue < p)] <- 'Signif. down-regulated'

    # keyvals[which(rlog_results$log2FoldChange > FC & rlog_results$pvalue < p)] <- 'red2'
    # names(keyvals)[which(rlog_results$log2FoldChange > FC & rlog_results$pvalue < p)] <- 'Signif. up-regulated'

    # unique(keyvals)
    # unique(names(keyvals))

    # volcanoplot_pval <- EnhancedVolcano(rlog_results,
    # lab = rownames(rlog_results),
    # x = 'log2FoldChange',
    # y = 'pvalue',
    # xlim = c(-6,6),
    # xlab = bquote(~Log[2]~ 'fold change'),
    # ylab = bquote(~-Log[10] ~ italic(P)),
    # title = "Volcano plot with the p-values",
    # subtitle = "",
    # titleLabSize = 15,
    # pCutoff = 10e-3,
    # FCcutoff = 0.6,
    # pointSize = 1.5,
    # labSize = 2.5,
    # colCustom = keyvals,
    # colAlpha = 1,
    # legendPosition = 'bottom',
    # legendLabSize = 10,
    # legendIconSize = 3.0,
    # #DrawConnectors = TRUE,
    # widthConnectors = 0.2,
    # colConnectors = 'grey50',
    # gridlines.major = FALSE,
    # gridlines.minor = FALSE,
    # border = 'partial',
    # borderWidth = 1.5,
    # borderColour = 'black')
    # plot(volcanoplot_pval)
 }
