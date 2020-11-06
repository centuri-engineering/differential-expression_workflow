# ibdm_rattier_rnaseq
B4| IBDM | Epithelial integrity controls the differentiation of hiPSCs into Primitive Streak by regulating TGF-b receptor accessibility

Applicant: Rosanna Dono & Diane Rattier

Institute: IBDM

Engineer: Thomas Vannier

Quick Start :
- Fill in the sample.tsv and coldata.tsv files to indicate the samples to analyse.
- Complete the config.yaml indicating the file and parameters to use.

Then you are ready to load the environment and run the workflow by using these commands in the  :
docker build -t ibdm_rattier_rnaseq /home/thomas/project/rattier/ibdm_rattier_rnaseq/
docker run -it -v ${PWD}:/ibdm_rattier_rnaseq ibdm_rattier_rnaseq
cd ibdm_rattier_rnaseq/
conda activate envfair
snakemake --use-conda --use-singularity -R
snakemake --report report.html


About this workflow :
This workflow run into a docker environment including conda.
The snakemake, python and panda used to run this workflow are in a conda environment (envfair.yaml)
Each snakemake rules call a specific conda environment. In this way you can easily change the tools using in the rules if necessary. 

3 steps for the analyses:
- clean.smk : 
- count.smk :
- differential_exp.smk :