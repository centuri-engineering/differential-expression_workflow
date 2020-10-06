SAMPLES = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())
# Remove the run with bad fastqc quality
SAMPLES_FILTERED = SAMPLES.remove(config["filtering"]["rmrun"])
index = config["ref"]["index"]

# ----------------------------------------------
# HISAT2-build: Indexing a reference genome
# ----------------------------------------------

rule hisat_build:
  output:
    config["ref"]["index"].1.ht2,
    config["ref"]["index"].2.ht2,
    config["ref"]["index"].3.ht2,
    config["ref"]["index"].4.ht2,
    config["ref"]["index"].5.ht2,
    config["ref"]["index"].6.ht2,
    config["ref"]["index"].7.ht2,
    config["ref"]["index"].8.ht2

  input:
    config["ref"]["genome"]

  conda: 
    "../02_Container/hisat2.yaml"

  shell:
    "hisat2-build {input} {output}" 

# ----------------------------------------------
# HISAT2: alignment of NGS reads to a population of genomes
# ----------------------------------------------
rule hisat:
  output:
    "05_Output/05_hisat/hisat.sam"

  input:
    fastq1 = expand( "05_Output/02_trimmomatic/{samples}_1.trimmed.fastq,", samples=samples.itertuples(),samples=SAMPLES_FILTERED),
    fastq2 = expand( "05_Output/02_trimmomatic/{samples}_2.trimmed.fastq,", samples=samples.itertuples(),samples=SAMPLES_FILTERED)

  conda: 
    "../02_Container/hisat2.yaml"

  shell:
    "hisat2 -x index -1 {input.fastq1} -2 {input.fastq2} -S 05_Output/05_hisat/hisat.sam " 