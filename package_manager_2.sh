#!/bin/bash

# Create a new conda environment
conda create --name bioconda-env --yes

# Activate the new environment
conda activate bioconda-env

# Add bioconda and conda-forge channels
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Install packages
conda install --yes \
    bowtie2 \
    cutadapt \
    samtools \
    picard \
    minimap2 \
    bwa-mem2 \
    trimmomatic \
    fastp \
    bedtools \
    delly \
    bcftools \
    gatk4 \
    snp-sites \
    fasttree \
    raxml-ng \
    iqtree

# Install bwa-index from the grst channel
conda install --yes -c grst bwa-index
