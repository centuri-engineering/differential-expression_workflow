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
    expand( "{samples.project}_{samples.condition}_{samples.sample}_1.trimmed.fastq", samples=samples.itertuples()),
    expand( "{samples.project}_{samples.condition}_{samples.sample}_2.trimmed.fastq", samples=samples.itertuples()),
    expand( "{samples.project}_{samples.condition}_{samples.sample}_1un.trimmed.fastq", samples=samples.itertuples()),
    expand( "{samples.project}_{samples.condition}_{samples.sample}_2un.trimmed.fastq", samples=samples.itertuples())

# ----------------------------------------------
# setup singularity 
# ----------------------------------------------

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

# ----------------------------------------------
# Load rules 
# ----------------------------------------------

include: "04_Workflow/clean.smk"
