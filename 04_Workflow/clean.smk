
# ----------------------------------------------
# FastQC to check the reads quality
# ----------------------------------------------
# rule fastqc:
#   output:
#     expand( "05_Output/01_fastqc/{samples.project}_{samples.condition}_{samples.sample}_{samples.paired}_fastqc.html", samples=samples.itertuples()),
#     expand( "05_Output/01_fastqc/{samples.project}_{samples.condition}_{samples.sample}_{samples.paired}_fastqc.zip", samples=samples.itertuples())

#   input:
#     expand( "00_RawData/{samples.project}_{samples.condition}_{samples.sample}_{samples.paired}.fastq", samples=samples.itertuples())

#   conda: 
#     "../02_Container/fastqc.yaml"

#   shell:
#     "fastqc --outdir 05_Output/01_fastqc/ {input}"


# ----------------------------------------------
# Trimmomatic: trimming reads and removing adapter sequences
# ----------------------------------------------
rule trimmomatic:
  input:
    expand( "00_RawData/{samples.project}_{samples.condition}_{samples.sample}_1.fastq", samples=samples.itertuples()),
    expand( "00_RawData/{samples.project}_{samples.condition}_{samples.sample}_2.fastq", samples=samples.itertuples())

  output:
    expand( "{samples.project}_{samples.condition}_{samples.sample}_1.trimmed.fastq", samples=samples.itertuples()),
    expand( "{samples.project}_{samples.condition}_{samples.sample}_2.trimmed.fastq", samples=samples.itertuples()),
    expand( "{samples.project}_{samples.condition}_{samples.sample}_1un.trimmed.fastq", samples=samples.itertuples()),
    expand( "{samples.project}_{samples.condition}_{samples.sample}_2un.trimmed.fastq", samples=samples.itertuples())

  conda: 
    "../02_Container/trimmomatic.yaml"

  shell:
    "trimmomatic PE -threads 4 {input} {output} LEADING:3 TRAILING:3 MINLEN:36"

