
# ----------------------------------------------
# FastQC to check the reads quality
# ----------------------------------------------
rule fastqc:
  output:
    expand( "05_Output/01_fastqc/{samples.project}_{samples.condition}_{samples.sample}_fastqc.html", samples=samples.itertuples()),
    expand( "05_Output/01_fastqc/{samples.project}_{samples.condition}_{samples.sample}_fastqc.zip", samples=samples.itertuples())

  input:
    expand( "00_RawData/{samples.project}_{samples.condition}_{samples.sample}.fastq", samples=samples.itertuples())

  conda: 
    "../02_Container/fastqc.yaml"

  shell:
    "fastqc --outdir 05_Output/01_fastqc/ {input}"


