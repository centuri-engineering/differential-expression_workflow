
# ----------------------------------------------
# Define some contants relative to project
# ----------------------------------------------
PROJECT_NAME="ibdm_rattier_rnaseq"
EXPERIMENT_NAME="fastqc"
SAMPLE_ID = ["NG-11683_Control_1_1_lib167668_5189_8_1.fastq.bz2","NG-11683_Control_1_1_lib167668_5189_8_2.fastq.bz2","NG-11683_Control_1_1_lib167668_5230_1_1.fastq.bz2","NG-11683_Control_1_1_lib167668_5230_1_2.fastq.bz2","NG-11683_Control_1_2_lib167669_5189_8_1.fastq.bz2","NG-11683_Control_1_2_lib167669_5189_8_2.fastq.bz2", "NG-11683_Control_1_3_lib167670_5189_8_1.fastq.bz2","NG-11683_Control_1_3_lib167670_5189_8_2.fastq.bz2","NG-11683_Control_1_3_lib167670_5230_1_1.fastq.bz2","NG-11683_Control_1_3_lib167670_5230_1_2.fastq.bz2","NG-11683_Control_2_1_lib167671_5189_8_1.fastq.bz2", "NG-11683_Control_2_1_lib167671_5189_8_2.fastq.bz2","NG-11683_Control_2_2_lib167672_5189_4_1.fastq.bz2","NG-11683_Control_2_2_lib167672_5189_4_2.fastq.bz2","NG-11683_Control_2_3_lib167673_5189_4_1.fastq.bz2","NG-11683_Control_2_3_lib167673_5189_4_2.fastq.bz2","NG-11683_Mutant_1_1_lib167674_5189_4_1.fastq.bz2","NG-11683_Mutant_1_1_lib167674_5189_4_2.fastq.bz2","NG-11683_Mutant_1_2_lib167675_5189_4_1.fastq.bz2" , "NG-11683_Mutant_1_2_lib167675_5189_4_2.fastq.bz2","NG-11683_Mutant_1_3_lib167676_5189_8_1.fastq.bz2","NG-11683_Mutant_1_3_lib167676_5189_8_2.fastq.bz2","NG-11683_Mutant_1_3_lib167676_5230_1_1.fastq.bz2","NG-11683_Mutant_1_3_lib167676_5230_1_2.fastq.bz2"]



# ----------------------------------------------
# The default target rule
# ----------------------------------------------
rule all:
  input:
    final_report = expand( "05_Output/01_fastqc/{project}_{experiment}_{sample}.html", project = PROJECT_NAME, experiment = EXPERIMENT_NAME, sample = SAMPLE_ID)

# ----------------------------------------------
# FastQC to check the reads quality
# ----------------------------------------------
rule fastqc:
  input:
    expand( "00_RawData/{sample}.fastq.bz2", sample = SAMPLE_ID)

  output:
    expand( "05_Output/01_fastqc/{project}_{experiment}_{sample}.html", project = PROJECT_NAME, experiment = EXPERIMENT_NAME, sample = SAMPLE_ID)

  singularity: 
    "02_Container/rna_velocity.img"

  shell:
    "fastqc --outdir 05_Output/01_fastqc/ {input}"








