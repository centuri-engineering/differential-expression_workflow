SAMPLES = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())
# Remove the run with bad fastqc quality
#SAMPLES_FILTERED = SAMPLES.remove(config["filtering"]["rmrun"])
index = config["ref"]["index"]
annotation = config["ref"]["annotation"]

# ----------------------------------------------
# HISAT2-build: Indexing a reference genome
# ----------------------------------------------

# rule hisat_build:
#   output:
#     config["ref"]["index"].1.ht2,
#     config["ref"]["index"].2.ht2,
#     config["ref"]["index"].3.ht2,
#     config["ref"]["index"].4.ht2,
#     config["ref"]["index"].5.ht2,
#     config["ref"]["index"].6.ht2,
#     config["ref"]["index"].7.ht2,
#     config["ref"]["index"].8.ht2

#   input:
#     config["ref"]["genome"]

#   conda: 
#     "../02_Container/hisat2.yaml"

#   shell:
#     "hisat2-build {input} {output}" 

# ----------------------------------------------
# HISAT2: alignment of NGS reads to a population of genomes
# ----------------------------------------------
# rule hisat:
#   output:
#     bam = expand( "05_Output/05_hisat/{samples}.bam", samples=SAMPLES)
    
#   input:
#     fastq1 = expand( "05_Output/02_trimmomatic/{samples}_1.trimmed.fastq", samples=SAMPLES),
#     fastq2 = expand( "05_Output/02_trimmomatic/{samples}_2.trimmed.fastq", samples=SAMPLES)

#   conda: 
#     "../02_Container/hisat2.yaml"

#   params:
#     sam = expand( "05_Output/05_hisat/{samples}.sam", samples=SAMPLES)

#   shell:
#     """
#     fastq1=({input.fastq1})
#     fastq2=({input.fastq2})
#     sam=({params.sam})
#     bam=({output.bam})
#     len=${{#fastq1[@]}}
#     for (( i=0; i<$len; i++ ))
#       do hisat2 -p 12 -x {index} -1 ${{fastq1[$i]}} -2 ${{fastq2[$i]}} -S ${{sam[$i]}}
#         samtools sort ${{sam[$i]}} > ${{bam[$i]}}
#         rm ${{sam[$i]}}
#     done
#     """


# ----------------------------------------------
# featureCounts: read count/summarization program
# ----------------------------------------------

# rule featureCounts:
#   output:
#     count = expand( "05_Output/06_featurecounts/{samples}_count.txt", samples=SAMPLES)
    
#   input:
#     {annotation},
#     bam = expand( "05_Output/05_hisat/{samples}.bam", samples=SAMPLES)

#   conda: 
#     "../02_Container/featureCounts.yaml"

#   shell:
#     """
#     bam=({input.bam})
#     count=({output.count})
#     len=${{#bam[@]}}
#     for (( i=0; i<$len; i++ ))
#       do featureCounts -T 12 -p -t exon -g gene_id -a {annotation} -o ${{count[$i]}} ${{bam[$i]}}
#     done
#     """

# ----------------------------------------------
# cpm: Counts Per Million and filtering
# ----------------------------------------------

rule cpm_filtering:
  output:
    count_mtx = "05_Output/07_cpm/count_matrix.txt"
    #cpm = "05_Output/07_cpm/cpm_matrix.txt"
    
  input:
    count = expand( "05_Output/06_featurecounts/{samples}_count.txt", samples=SAMPLES)
  conda: 
    "../02_Container/cpm.yaml"

  script:
    "../03_Script/cpm_filtering.R"