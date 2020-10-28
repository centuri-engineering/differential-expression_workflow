This workflow performs differential expression analysis on paired-end RNA-seq data. 
After adapter removal with Cutadapt, reads were mapped and gene counts were generated with Hisat and featurecounts.
Gene counts of replicated were summed up.
Integrated normalization and differential expression analysis was conducted with DESeq2 following standard procedure as outlined in the manual.