# Differential Expression Workflow: RNA-seq analyse

## Author

Thomas Vannier (@metavannier), https://centuri-livingsystems.org/t-vannier/

## About

This workflow performs an RNA-seq analyses from the sequencing output data to the differential expression analyses. 
It run into a docker container (see Dockerfile) including a general conda environment (see envfair.yaml).
Each snakemake rules call a specific conda environment. In this way you can easily change/add tools for each step if necessary. 

3 steps for the analyses:
- clean.smk : 
- count.smk :
- differential_exp.smk :

## Usage

### Step 1: Install workflow

You can use this workflow by downloading and extracting the latest release. If you intend to modify and further extend this workflow or want to work under version control, you can fork this repository.

We would be pleased if you use this workflow and participate in its improvement. If you use it in a paper, don't forget to give credits to the author by citing the URL of this repository and, if available, its DOI (see above).

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files:
- sample.tsv and coldata.tsv to indicate the samples, run, condition, etc. for the analyse.
- config.yaml indicating the parameters to use.

### Step 3: Execute workflow

You just need to have Docker installed on your computer.

- Load the docker environment and run the workflow by using these commands  :
docker build -t de_workflow /home/thomas/project/rattier/ibdm_rattier_rnaseq/
docker run -it -v ${PWD}:/de_workflow de_workflow

- Then execute the workflow locally via
cd de_workflow/
snakemake --use-conda -R

### Step 4: Investigae results 

After successful execution, you can create a self-contained interactive HTML report with all results via:

snakemake --report report.html


