# Use a base image with Miniconda installed
FROM continuumio/miniconda3

# Create the working directory inside the container
WORKDIR /app

# Create a new Conda environment and install packages
RUN conda create --name bioconda-env --yes && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict && \
    conda install --name bioconda-env --yes \
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
        iqtree && \
    conda install --name bioconda-env --yes -c grst bwa-index && \
    # Clean up to reduce image size
    conda clean --all --yes && \
    rm -rf /var/cache/apt/*

# Set environment variable to activate the Conda environment by default
ENV PATH /opt/conda/envs/bioconda-env/bin:$PATH

# Set the default command to run a shell
CMD ["/bin/bash"]
