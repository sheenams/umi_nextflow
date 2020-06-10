# Using a biocontainer base image
# Please see below for further details : 
# https://github.com/BioContainers/containers/blob/master/biocontainers/1.1.0/Dockerfile
FROM biocontainers/biocontainers:v1.1.0_cv1

################ Create appropriate directory structure for code to run ################
USER root
RUN mkdir -p /srv/qgen/code && \
<<<<<<< HEAD
    mkdir -p /srv/qgen/bin/downloads 
=======
    mkdir -p /srv/qgen/bin/downloads
>>>>>>> master

################ Update package repository and install dependencies using apt-get ################
RUN apt-get -y update && \
    apt-get -y install r-base libcurl4-openssl-dev libssl-dev libxml2-dev && \
    apt-get clean && apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

################ Install various version specific 3rd party tools ################
RUN conda install bedtools=2.25.0

################ Install python modules ################
## Install some modules with conda
RUN conda install scipy pysam=0.9.0
    
################ R packages ################
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages(c('plyr','tidyverse','magrittr','data.table'))"

<<<<<<< HEAD
################ Install source code ################
ADD smc2-af0bd7.tar.gz /srv/qgen/code/
ADD data.tar.gz /srv/qgen/

################ Update Environment Variables ################
ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python2.7/site-packages/:/srv/qgen/code/smc2-af0bd7/
=======
################ Update Environment Variables ################
ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python2.7/site-packages/:/srv/qgen/code/qiaseq-dna/
>>>>>>> master

