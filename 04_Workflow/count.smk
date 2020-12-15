SAMPLES = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())

index = config["ref"]["index"]
annotation = config["ref"]["annotation"]

# ----------------------------------------------
# HISAT2-build: Indexing a reference genome
# ----------------------------------------------
rule hisat_build:
  output:
    expand( "{index}.1.ht2", index=index),
    expand( "{index}.2.ht2", index=index),
    expand( "{index}.3.ht2", index=index),
    expand( "{index}.4.ht2", index=index),
    expand( "{index}.5.ht2", index=index),
    expand( "{index}.6.ht2", index=index),
    expand( "{index}.7.ht2", index=index),
    expand( "{index}.8.ht2", index=index)

  input:
    config["ref"]["genome"]

  conda: 
    "../02_Container/hisat2.yaml"

  shell:
    """
    index=({config.ref.index})
    "hisat2-build {input} ${{index}}"
    """

# ----------------------------------------------
# HISAT2: alignment of NGS reads to a population of genomes
# ----------------------------------------------
rule hisat:
  output:
    bam = expand( "05_Output/05_hisat/{samples}.bam", samples=SAMPLES)
    
  input:
    fastq1 = expand( "05_Output/02_trimmomatic/{samples}_1.trimmed.fastq", samples=SAMPLES),
    fastq2 = expand( "05_Output/02_trimmomatic/{samples}_2.trimmed.fastq", samples=SAMPLES)

  conda: 
    "../02_Container/hisat2.yaml"

  params:
    sam = expand( "05_Output/05_hisat/{samples}.sam", samples=SAMPLES)

  shell:
    """
    fastq1=({input.fastq1})
    fastq2=({input.fastq2})
    sam=({params.sam})
    bam=({output.bam})
    len=${{#fastq1[@]}}
    for (( i=0; i<$len; i++ ))
      do hisat2 -p 12 -x {index} -1 ${{fastq1[$i]}} -2 ${{fastq2[$i]}} -S ${{sam[$i]}}
        samtools sort ${{sam[$i]}} > ${{bam[$i]}}
        rm ${{sam[$i]}}
    done
    """

# ----------------------------------------------
# featureCounts: read count/summarization program
# ----------------------------------------------

rule featureCounts:
  output:
    count = expand( "05_Output/06_featurecounts/{samples}_count.txt", samples=SAMPLES)
    
  input:
    {annotation},
    bam = expand( "05_Output/05_hisat/{samples}.bam", samples=SAMPLES)

  conda: 
    "../02_Container/featureCounts.yaml"

  shell:
    """
    bam=({input.bam})
    count=({output.count})
    len=${{#bam[@]}}
    for (( i=0; i<$len; i++ ))
      do featureCounts -T 12 -p -t exon -g gene_id -a {annotation} -o ${{count[$i]}} ${{bam[$i]}}
    done
    """

# ----------------------------------------------
# cpm: Counts Per Million and filtering
# ----------------------------------------------

rule cpm_filtering:
  output:
    count_df = report("05_Output/07_cpm/count.txt", caption="../report/count.rst", category="02 Count matrices"),
    output_filter_count = report("05_Output/07_cpm/count_filtered.txt", caption="../report/count_filtered.rst", category="02 Count matrices"),
    cpm = report("05_Output/07_cpm/cpm_filtered.txt", caption="../report/cpm_filtered.rst", category="02 Count matrices"),
    lengths = "05_Output/07_cpm/lengths.txt"

  input:
    expand( "05_Output/06_featurecounts/{samples}_count.txt", samples=SAMPLES)
    
  params:
    thresh_cpm = config["filtering"]["thresh_cpm"],
    thresh_sample = config["filtering"]["thresh_sample"],
    rmrun_list = config["filtering"]["rmrun"],
    path = "05_Output/06_featurecounts"

  conda: 
    "../02_Container/cpm.yaml"

  script:
    "../03_Script/cpm_filtering.R"
