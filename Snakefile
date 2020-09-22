import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="06_Schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index(["project", "condition", "sample"], drop=False)
validate(samples, schema="06_Schemas/samples.schema.yaml")

# ----------------------------------------------
# Target rules
# ----------------------------------------------
rule all:
  input:
    expand( "05_Output/01_fastqc/{samples.project}_{samples.condition}_{samples.sample}_fastqc.html", samples=samples.itertuples()),
    expand( "05_Output/01_fastqc/{samples.project}_{samples.condition}_{samples.sample}_fastqc.zip", samples=samples.itertuples())

# ----------------------------------------------
# setup singularity 
# ----------------------------------------------

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

# ----------------------------------------------
# Load rules 
# ----------------------------------------------

include: "04_Workflow/fastqc.smk"



# ----------------------------------------------
# FastQC to check the reads quality
# ----------------------------------------------
# rule fastqc:
#   output:
#     expand( "05_Output/01_fastqc/{sample}_fastqc.html", sample = SAMPLE_ID),
#     expand( "05_Output/01_fastqc/{sample}_fastqc.zip", sample = SAMPLE_ID)

#   input:
#     expand( "00_RawData/{sample}.fastq", sample = SAMPLE_ID)

#   conda: 
#     "02_Container/fastqc.yaml"

#   shell:
#     "fastqc --outdir 05_Output/01_fastqc/ {input}"