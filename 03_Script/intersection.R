# #################################################################
# 
#        Presence or absence of differential expressed 
#            genes in different comparisons
#              
# #################################################################

file_1 = snakemake@params[["gene_list_1"]]
file_2 = snakemake@params[["gene_list_2"]]
file_1=gsub("^.*/", "", file_1)
file_2=gsub("^.*/", "", file_2)

gene_list_1 = read.delim(file = snakemake@params[["gene_list_1"]],
                           header = FALSE,
                           stringsAsFactors = FALSE)
                           
gene_list_2 = read.delim(file = snakemake@params[["gene_list_2"]],
                           header = FALSE,
                           stringsAsFactors = FALSE)

## Data frame with common genes and corresponding fold change
# Remove the second column (basemean)
gene_list_1 = gene_list_1[,-2]
gene_list_2 = gene_list_2[,-2]

gene_list_1$foldchange_2=gene_list_2[,2][match(gene_list_1[,1], gene_list_2[,1])]
gene_list_1$pvalue_2=gene_list_2[,5][match(gene_list_1[,1], gene_list_2[,1])]
gene_list_1[,3]=gene_list_1[,5]
gene_list_1[,4]=gene_list_1$foldchange_2
gene_list_1[,5]=gene_list_1$pvalue_2

gene_list=gene_list_1[!is.na(gene_list_1[,4]),]
gene_list <- gene_list[-1,]
colnames(gene_list)[1] = "gene_id"
colnames(gene_list)[2] = paste("log2foldchange",file_1,sep = " ")
colnames(gene_list)[3] = paste("padj",file_1,sep = " ")
colnames(gene_list)[4] = paste("log2foldchange",file_2,sep = " ")
colnames(gene_list)[5] = paste("padj",file_2,sep = " ")

write.table(gene_list[,1:5], file =snakemake@params[["common_genes"]],row.names=FALSE,col.names = TRUE,sep="\t", quote = FALSE)

## Foldgene of the common genes in different selected comparison

stat_file_1 = snakemake@params[["stat_file_1"]]
stat_file_1=gsub("^.*/", "", stat_file_1)
stat_file_2 = snakemake@params[["stat_file_2"]]
stat_file_2=gsub("^.*/", "", stat_file_2)

stat_1 = read.delim(file = snakemake@params[["stat_file_1"]],
                           header = FALSE,
                           stringsAsFactors = FALSE)
                           
stat_2 = read.delim(file = snakemake@params[["stat_file_2"]],
                           header = FALSE,
                           stringsAsFactors = FALSE)

stat_1 = stat_1[,-2]
stat_2 = stat_2[,-2]

gene_list$foldchange_stat_1=stat_1[,2][match(gene_list$gene_id, stat_1[,1])]
gene_list$pvalue_1=stat_1[,5][match(gene_list$gene_id, stat_1[,1])]

gene_list$foldchange_stat_2=stat_2[,2][match(gene_list$gene_id, stat_2[,1])]
gene_list$pvalue_2=stat_2[,5][match(gene_list$gene_id, stat_2[,1])]

gene_list[,2]=gene_list$foldchange_stat_1
gene_list[,3]=gene_list$pvalue_1
gene_list[,4]=gene_list$foldchange_stat_2
gene_list[,5]=gene_list$pvalue_2

colnames(gene_list)[2] = paste("log2foldchange",stat_file_1,sep = " ")
colnames(gene_list)[3] = paste("padj",stat_file_1,sep = " ")
colnames(gene_list)[4] = paste("log2foldchange",stat_file_2,sep = " ")
colnames(gene_list)[5] = paste("padj",stat_file_2,sep = " ")

write.table(gene_list[,1:5], file =snakemake@params[["stat_common_genes"]],row.names=FALSE,col.names = TRUE,sep="\t", quote = FALSE)

# File to complete the snakemake rule
intersection_output=snakemake@output[["intersection_output"]]
writeLines(c("intersection step done"), intersection_output)