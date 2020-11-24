import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="06_Schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index(["project", "condition", "sample"], drop=False)
validate(samples, schema="06_Schemas/samples.schema.yaml")

coldata = pd.read_table(config["coldata"]).set_index(["condition"], drop=False)
validate(coldata, schema="06_Schemas/coldata.schema.yaml")

# ----------------------------------------------
# Target rules
# ----------------------------------------------
SAMPLES = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())
ref_level = config["diffexp"]["ref_level"]

rule all:
  input:
    # expand( "05_Output/03_fastqc/{samples}_{ext}.trimmed_fastqc.html", samples=SAMPLES, ext=["1","2"]),
    # "05_Output/07_cpm/count.tsv",
    # "05_Output/07_cpm/cpm_filtered.tsv",
    # "05_Output/07_cpm/count_filtered.txt",
    # "05_Output/08_deseq2_init/all.rds",
    "05_Output/09_differential_expression/diffexp.html",
    expand("05_Output/09_differential_expression/{coldata.condition}_vs_{ref_level}.csv",coldata=coldata.itertuples(),ref_level=ref_level)
# ----------------------------------------------
# setup singularity 
# ----------------------------------------------

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
#singularity: "docker://continuumio/miniconda3"


# ----------------------------------------------
# setup report
# ----------------------------------------------

report: "report/workflow.rst"

# ----------------------------------------------
# Load rules 
# ----------------------------------------------

include: "04_Workflow/clean.smk"
include: "04_Workflow/count.smk"
include: "04_Workflow/diffexp.smk"
