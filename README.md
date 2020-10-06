# ibdm_rattier_rnaseq
B4| IBDM | Epithelial integrity controls the differentiation of hiPSCs into Primitive Streak by regulating TGF-b receptor accessibility

Applicant: Rosanna Dono & Diane Rattier

Institute: IBDM

Engineer: Thomas Vannier

docker run -i -t -v ${PWD}:/ibdm_rattier_rnaseq continuumio/miniconda3
cd ibdm_rattier_rnaseq/
conda env create -n envfair -f 02_Container/envfair.yaml
conda activate envfair
snakemake --use-conda -R