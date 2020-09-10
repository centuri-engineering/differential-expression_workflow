# Based on rocker - https://github.com/rocker-org/rocker-versioned
FROM rocker/tidyverse:3.6.0

MAINTAINER Thomas Vannier (thomas.vannier@univ-amu.fr)

#############
# ibdm_rattier_rnaseq
#############


# ##########
# R PACKAGES 
# ##########

#### Figures & layout management
# ggplot2
RUN R -e 'install.packages( "ggplot2")'
RUN R -e 'install.packages( "cowplot")'        # plot_grid, themes, ...
RUN R -e 'install.packages( "ggpubr")'         # add_summary, geom_signif, ...
RUN R -e 'install.packages( "ggrepel")'        # geom_text_repel, geom_label_repel
RUN R -e 'install.packages( "gridExtra")'      # grid.arrange, ...
RUN R -e 'install.packages( "patchwork")'      # +/ operators for ggplots

# plotly
RUN R -e 'install.packages( "plotly")'

# general
RUN R -e 'install.packages( "gplots")'         # heatmap.2
RUN R -e 'install.packages( "heatmaply")'      # heatmaply (interactive)
RUN R -e 'BiocManager::install( "iheatmapr")'  # iheatmap (interactive, uses plotly)
RUN R -e 'install.packages( "pheatmap")'       # pheatmap
RUN R -e 'install.packages( "venn")'           # venn
RUN R -e 'install.packages( "VennDiagram")'    # venn.diagram
RUN R -e 'BiocManager::install("qvalue")'
RUN R -e 'install.packages( "fastcluster")'    # fastcluster
RUN R -e 'install.packages( "pvclust")'        # pvclust
RUN R -e 'install.packages("factoextra")'      # prcomp (ACP)

#### Reporting
RUN R -e 'install.packages( "DT")'             # datatable
RUN R -e 'install.packages( "pander")'         # pander


#### General
RUN R -e 'install.packages( "funr")'           # get_script_path
RUN R -e 'install.packages( "reshape")'        # melt
RUN R -e 'install.packages( "readr")'          # write_csv2
RUN R -e 'install.packages( "stringr")'        # str_to_lower
RUN R -e 'install.packages( "tidyverse)")'     # select


#### Custom

# ####################################
# INSTALLING BASIC PACKAGES
# ####################################

# -- Install vim
RUN apt-get update \
  && apt-get -y install vim

# ##################
# OTHER DEPENDENCIES
# ##################

#### Python
# pip
RUN apt-get update && apt-get install -y \
    python-pip \
    && pip install --upgrade --user pip

# ##################
# FASTQC
# ##################

# fastqc requires java
RUN apt-get update && apt-get install -y \
  curl \
  unzip \
  perl \
  openjdk-8-jre-headless

# Installs fastqc from compiled java distribution into /opt/FastQC
ENV FASTQC_URL http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
ENV FASTQC_VERSION 0.11.9
ENV FASTQC_RELEASE fastqc_v${FASTQC_VERSION}.zip
ENV DEST_DIR /opt/

# Make destination directory
RUN mkdir -p $DEST_DIR

# Download & extract archive - Repo includes binaries for linux
WORKDIR /tmp

# Do this in one command to avoid caching the zip file and its removal in separate layers
RUN curl -SLO ${FASTQC_URL}/${FASTQC_RELEASE} && unzip ${FASTQC_RELEASE} -d ${DEST_DIR} && rm ${FASTQC_RELEASE}

# Make the wrapper script executable
RUN chmod a+x ${DEST_DIR}/FastQC/fastqc

# Include it in PATH
ENV PATH ${DEST_DIR}/FastQC:$PATH


