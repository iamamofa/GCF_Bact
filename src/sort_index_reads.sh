#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_sam> <output_dir>"
    exit 1
fi

input_sam=$1
output_dir=$2

# Check if input file exists
if [ ! -f "$input_sam" ]; then
    echo "Error: Input SAM file '$input_sam' not found!"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Function to install Conda
install_conda() {
    echo "Installing Conda..."
    if [ "$(uname -s)" = "Linux" ]; then
        curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda
    elif [ "$(uname -s)" = "Darwin" ]; then
        curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
        bash Miniconda3-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda
    elif [[ "$OS" == "Windows_NT" ]]; then
        echo "Please install Miniconda manually from https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    else
        echo "Unsupported OS for automatic Conda installation."
        exit 1
    fi
    export PATH="$HOME/miniconda/bin:$PATH"
    conda init
    source ~/.bashrc
}

# Function to install Samtools using Conda
install_samtools() {
    echo "Installing Samtools..."
    conda create -n samtools_env -c bioconda samtools -y
    conda activate samtools_env
    echo "Samtools installation complete."
}

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    install_conda
else
    echo "Conda is already installed."
fi

# Activate the Conda environment and install Samtools
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate samtools_env || conda activate
if ! command -v samtools &> /dev/null; then
    install_samtools
else
    echo "Samtools is already installed."
fi

# Convert SAM to BAM, sort, and index the reads
echo "Converting SAM to BAM and sorting..."
if samtools view -Sb "$input_sam" | samtools sort -o "${output_dir}/sorted_reads.bam"; then
    echo "Indexing sorted BAM file..."
    if samtools index "${output_dir}/sorted_reads.bam"; then
        echo "Sorting and indexing complete. BAM and index files stored in $output_dir"
    else
        echo "Error: Samtools indexing failed!"
        exit 1
    fi
else
    echo "Error: Samtools conversion and sorting failed!"
    exit 1
fi
