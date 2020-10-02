SAMPLES = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())

# ----------------------------------------------
# HISAT2: alignment of NGS reads to a population of genomes
# ----------------------------------------------
rule hisat:
  output:
    expand( "05_Output/04_hisat2/{samples}_{ext}.trimmed_fastqc.html", samples=SAMPLES, ext=["1","2"])

  input:
    expand( "05_Output/02_trimmomatic/{samples}_{ext}.trimmed.fastq", samples=SAMPLES, ext=["1","2"])

  conda: 
    "../02_Container/hisat2.yaml"

  shell:
    "fastqc --outdir 05_Output/03_fastqc_qual2/ {input}"