FROM --platform=linux/amd64 ubuntu:18.04

MAINTAINER Guanliang Meng <linzhi2012[AT]gmail[DOT].com>, ZFMK

ENV DEBIAN_FRONTEND=noninteractive


# Install required packages
RUN apt-get update
RUN apt-get install -y  wget bzip2 bash


# install anaconda
RUN mkdir /app

RUN if [ ! -d /app/anaconda ]; then \
    wget -c https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh -O /app/anaconda.sh && \
    bash /app/anaconda.sh -b -p /app/anaconda && \
    rm -rf /app/anaconda.sh ; fi

RUN apt-get clean && \
    apt-get autoclean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# set anaconda path
ENV PATH="/app/anaconda/bin:$PATH"

RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install mamba -n base -c conda-forge && \
    mamba install -c bioconda mitoz=3.4 && \
    mamba clean -y -a && \
    conda clean -y -a

# install NCBI taxonomy database
RUN cd / ;  python3 -c 'from ete3 import NCBITaxa; ncbi = NCBITaxa()' ; rm -rf taxdump.tar.gz

ENV LC_ALL=C
ENV PATH=/app/anaconda/bin:$PATH

VOLUME /project
WORKDIR /project

# https://app.deepsource.com/directory/analyzers/docker/issues/DOK-DL4005
# Use bash as the default shell
SHELL ["/bin/bash", "-c"]

CMD ["/bin/bash"]
