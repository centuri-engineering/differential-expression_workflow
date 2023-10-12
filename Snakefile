# Docker container based on a minimal Ubuntu installation that includes conda-forge's mambaforge installer.
# container: "docker://condaforge/mambaforge"

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

##### Set variables ####
ROOTDIR = os.getcwd()
RAWDATA = srcdir("00_RawData/")
REF = srcdir("01_Reference/")
CONTAINER = srcdir("02_Container/")
SCRIPTDIR = srcdir("03_Script/")
ENVDIR = srcdir("04_Workflow/")
OUTPUTDIR = srcdir("05_Output/")
REPORT = srcdir("07_Report/")

# ----------------------------------------------
# Target rules
# ----------------------------------------------
SAMPLES = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())
RUN =  config["run"]["type"].split(',')
EXT = config["run"]["ext"]
ref_level = config["diffexp"]["ref_level"]
genome = config["ref"]["genome"]
index = config["ref"]["index"]
annotation = config["ref"]["annotation"]
RUN_ID = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())

rule all:
  input:
    fastqc_output="05_Output/01_fastqc/fastqc_output.txt",
    trimmomatic_output= OUTPUTDIR + "02_trimmomatic/trimmomatic_output.txt",
    fastqc_trimmed_output= OUTPUTDIR + "03_fastqc/fastqc_trimmed_output.txt",
    multiqc_trimmed_output = OUTPUTDIR + "03_fastqc/trimmed_multiqc.txt",
    index1 = expand( OUTPUTDIR + "{index}.1.ht2", index=index),
    index2 = expand( OUTPUTDIR + "{index}.2.ht2", index=index),
    index3 = expand( OUTPUTDIR + "{index}.3.ht2", index=index),
    index4 = expand( OUTPUTDIR + "{index}.4.ht2", index=index),
    index5 = expand( OUTPUTDIR + "{index}.5.ht2", index=index),
    index6 = expand( OUTPUTDIR + "{index}.6.ht2", index=index),
    index7 = expand( OUTPUTDIR + "{index}.7.ht2", index=index),
    index8 = expand( OUTPUTDIR + "{index}.8.ht2", index=index),
    hisat_output=OUTPUTDIR + "05_hisat/hisat.txt",
    featureCounts_output = OUTPUTDIR + "06_featurecounts/featurecounts.txt",
    cpm_filtering_output = OUTPUTDIR + "07_cpm/cpm_filtering.txt",
    deseq2_init_output = OUTPUTDIR + "08_deseq2_init/deseq2_init.txt",
    diffexp_output = OUTPUTDIR + "09_differential_expression/diffexp.txt"

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

# ruleorder: trimmomatic > fastqc_trimmed > hisat_build > hisat > featureCounts > cpm_filtering > deseq2_init > diffexp

# ----------------------------------------------
# Load rules 
# ----------------------------------------------

include: ENVDIR + "clean.smk"
include: ENVDIR + "count.smk"
include: ENVDIR + "diffexp.smk"
