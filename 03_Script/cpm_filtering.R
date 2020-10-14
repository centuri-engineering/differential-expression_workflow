# #################################################################
# 
#        Calcul the cpm (Counts Per Million) and
#            filter the low expressed genes
#              
# #################################################################

#library("cpm")

files_count <- as.data.frame(snakemake@input[["count"]])

# Dataframe with the gene count for each sample
dataframe_total_count <- data.frame()
# Add column with the gene name
gene_name <- read.delim(as.character(files_count[[1]]), comment.char="#")
dataframe_total_count <- gene_name["Geneid"]
for (i in 1:length(files_count)) {
  file_count <- as.character(files_count[[i]])
  count <- read.delim(file_count, comment.char="#")
  dataframe_fd_mouse <- cbind(count[7,])
}

#cts <- read.delim(snakemake@input[["count"]], comment.char="#")
#saveRDS(cts, file=snakemake@output["count_mtx"])
class(count)

