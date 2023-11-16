This workflow performs differential expression analysis on paired-end RNA-seq data.

Materials and methods
---------------------

We follow the data analysis procedure performed in Ghosh et al. (doi: 10.1016/j.compbiolchem.2020.107239). Version of tools and genomes have been updated and parameter used have been adapted according to our data.

The quality of the raw reads were assessed using FastQC v0.12.1 toolkit (Andrews, 2010). Adapters and low quality reads were trimmed using Trimmomatic v0.39 (Bolger et al., 2014).

HiSat2 v2.2.1 (Kim et al., 2015) was then used for mapping the raw reads to the rreference genome Drosophila melanogaster (RefSeq:GCF_000001215.4). The expression for each gene were evaluated using featureCounts from the Subread v2.0.1 package (Liao et al., 2019). 

The low expressed genes which did not have more than one count per million reads (1CPM) in at least three samples within each datasets were removed from further analysis. The raw counts were then normalized and used for differential expression testing using DESeq2 v1.28.0 (Love et al., 2014). 

For differential RNA-seq analysis, genes with a Fold change between −2 to 2 and an adjusted P-value (Padj) > 0.01 were not considered as significantly expressed. Volcanoplot was generated with the software EnhancedVolcano.