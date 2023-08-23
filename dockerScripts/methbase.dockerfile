# Dockerfile for Methylation Pipeline Base Packages
# 
# Author: Jonathan Serrano
# Organization: NYU Langone Health
# Date Created: August 23, 2023
# Last Updated: August 23, 2023
# Version: 1.0.0
# Description: Loads Base R 4.3.1 environment for all packages required for the Methylation Pipeline
# Copyright NYULH 2023

FROM rocker/verse:4.3.1

# Set JAVA_HOME environment variable
ENV JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64

# Install system dependencies, R packages, Github packages, and Bioconductor packages
RUN apt-get update && apt-get install -y \
    libmysqlclient-dev libgdal-dev libproj-dev gcc g++ libsqlite3-dev sqlite3 apt-transport-https gfortran libatlas-base-dev libbz2-dev libicu-dev liblzma-dev libpango-1.0-0 libpangocairo-1.0-0 libtcl8.6 libtiff5 libtk8.6 libudunits2-dev libx11-6 libxt6 zlib1g-dev openjdk-11-jdk libglu1-mesa-dev littler \
    && ln -s /usr/local/lib/R/site-library/littler/examples/installBioc.r /usr/local/bin/installBioc.r \
    && install2.r --error --skipinstalled --ncpus=-1 -d TRUE \
    tinytest ncdf4 sf deldir XML leaflet \
    ade4 randomForest glmnet tidyverse gridExtra knitr Rtsne pander plotly BiocManager compiler crayon fs htmltools R.utils beepr rmarkdown mlr RPMM \
    && install2.r --error --skipinstalled --ncpus=-1 --repos=https://rspatial.r-universe.dev -d TRUE terra \
    && installGithub.r -d TRUE mwsill/IlluminaHumanMethylationEPICmanifest \
    && R -e "install.packages(c('https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz', 'https://cran.r-project.org/src/contrib/Archive/FField/FField_0.1.0.tar.gz'), repos = NULL, dependencies = TRUE, verbose = TRUE, type = 'source', ask = FALSE)" \
    && installBioc.r --error --skipinstalled --ncpus=-1 -d TRUE \
    impute MethylAid BiocVersion lumi methylumi conumee minfi IlluminaHumanMethylation450kmanifest wateRmelon BiocParallel \
    bumphunter locfit iterators foreach Biostrings XVector SummarizedExperiment DelayedArray matrixStats \
    Biobase GenomicRanges GenomeInfoDb IRanges S4Vectors BiocGenerics FDb.InfiniumMethylation.hg19 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*
