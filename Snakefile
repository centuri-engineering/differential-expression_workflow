import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="06_Schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index(["project", "condition", "sample", "fq1", "fq2"], drop=False)
validate(samples, schema="06_Schemas/samples.schema.yaml")

# ----------------------------------------------
# Target rules
# ----------------------------------------------
rule all:
  input:
    expand( "05_Output/03_fastqc/{samples}_{ext}_fastqc.trimmed.html", samples=samples.itertuples(),samples=SAMPLES, ext=["1","2"]),
    expand( "05_Output/03_fastqc/{samples}_{ext}_fastqc.trimmed.zip", samples=samples.itertuples(),samples=SAMPLES, ext=["1","2"])


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
