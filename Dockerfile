FROM continuumio/miniconda3

ADD /envfair.yaml /tmp/envfair.yaml

# Pull the environment name out of the envfair.yaml
RUN conda env create -f /tmp/envfair.yaml -n envfair

ENV PATH /opt/conda/envs/envfair/bin:$PATH
# PATH pour lire le unsquashfs de singularity mais ne marche pas (voir si modification du .conf?)
#ENV PATH /opt/conda/envs/envfair/etc/singularity:$PATH
# PATH for r packages
ENV PATH /opt/conda/envs/envfair/lib/R/library:$PATH
ENV PATH /opt/conda/envs/envfair/lib/python3.7/site-packages:$PATH

# Activation de directement l'environnement
RUN echo "source activate envfair" > ~/.bashrc
# Run the workflow
#RUN echo "snakemake --use-conda --use-singularity -R" >> ~/.bashrc



