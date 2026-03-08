# Base image
FROM continuumio/miniconda3:latest

# Metadata
LABEL maintainer="Shohei Noda"
LABEL description="Container image with BWA and Samtools for reference indexing"

# Channel setup
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# Installation of BWA and Samtools
RUN conda install -y \
    bwa=0.7.17 \
    samtools=1.19 \
    && conda clean -a -y

# Bash as default command
CMD ["/bin/bash"]

# Build the Docker image
# docker build --platform linux/amd64 -t bwa_samtools:latest .
# docker tag bwa_samtools:latest shohei0212/bwa_samtools:1.0.0
# docker push shohei0212/bwa_samtools:1.0.0