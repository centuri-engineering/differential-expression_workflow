import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="06_Schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index(["project", "condition", "sample"], drop=False)
validate(samples, schema="06_Schemas/samples.schema.yaml")

coldata = pd.read_table(config["coldata"]).set_index(["project", "condition", "type"], drop=False)
validate(coldata, schema="06_Schemas/coldata.schema.yaml")

condition = pd.read_table(config["condition"]).set_index(["condition"], drop=False)
validate(coldata, schema="06_Schemas/condition.schema.yaml")

# ----------------------------------------------
# Target rules
# ----------------------------------------------
SAMPLES = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())
ref_level = config["diffexp"]["ref_level"]
index = config["ref"]["index"]
annotation = config["ref"]["annotation"]
RUN_ID = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())

rule all:
  input:
    expand( "05_Output/02_trimmomatic/{samples}_{ext}.trimmed.fastq", samples=SAMPLES, ext=["1","2"]),
    expand( "05_Output/02_trimmomatic/{samples}_{ext}un.trimmed.fastq", samples=SAMPLES, ext=["1","2"]),
    expand( "05_Output/03_fastqc/{samples}_{ext}.trimmed_fastqc.html", samples=SAMPLES, ext=["1","2"]),
    expand( "05_Output/03_fastqc/{samples}_{ext}.trimmed_fastqc.zip", samples=SAMPLES, ext=["1","2"]),
    expand( "{index}.1.ht2", index=index),
    expand( "{index}.2.ht2", index=index),
    expand( "{index}.3.ht2", index=index),
    expand( "{index}.4.ht2", index=index),
    expand( "{index}.5.ht2", index=index),
    expand( "{index}.6.ht2", index=index),
    expand( "{index}.7.ht2", index=index),
    expand( "{index}.8.ht2", index=index),
    expand( "05_Output/05_hisat/{samples}.bam", samples=SAMPLES),
    expand( "05_Output/06_featurecounts/{samples}_count.txt", samples=SAMPLES),
    "05_Output/07_cpm/count.txt",
    "05_Output/07_cpm/count_filtered.txt",
    "05_Output/07_cpm/cpm_filtered.txt",
    "05_Output/08_deseq2_init/all.rds",
    "05_Output/08_deseq2_init/normalized_counts.tsv",
    "05_Output/09_differential_expression/diffexp.html",
    expand("05_Output/09_differential_expression/{condition.condition}_vs_{ref_level}_all_genes_stats.tsv",condition=condition.itertuples(),ref_level=ref_level),
    expand("05_Output/09_differential_expression/{condition.condition}_vs_{ref_level}_signif-up-regulated.txt", condition=condition.itertuples(), ref_level=ref_level),
    expand("05_Output/09_differential_expression/{condition.condition}_vs_{ref_level}_signif-down-regulated.txt", condition=condition.itertuples(), ref_level=ref_level)

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
# Impose rule order for the execution of the workflow 
# ----------------------------------------------

ruleorder: trimmomatic > fastqc_trimmed > hisat_build > hisat > featureCounts > cpm_filtering > deseq2_init > diffexp

# ----------------------------------------------
# Load rules 
# ----------------------------------------------

include: "04_Workflow/clean.smk"
include: "04_Workflow/count.smk"
include: "04_Workflow/diffexp.smk"
