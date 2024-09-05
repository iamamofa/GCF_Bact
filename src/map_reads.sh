#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_fastq> <reference> <output_dir>"
    exit 1
fi

input_fastq=$1
reference=$2
output_dir=$3

# Check if input file and reference file exist
if [ ! -f "$input_fastq" ]; then
    echo "Error: Input FASTQ file '$input_fastq' not found!"
    exit 1
fi

if [ ! -f "$reference" ]; then
    echo "Error: Reference file '$reference' not found!"
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

# Function to install BWA using Conda
install_bwa() {
    echo "Installing BWA..."
    conda create -n bwa_env -c bioconda bwa -y
    conda activate bwa_env
    echo "BWA installation complete."
}

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    install_conda
else
    echo "Conda is already installed."
fi

# Activate the Conda environment and install BWA
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bwa_env || conda activate
if ! command -v bwa &> /dev/null; then
    install_bwa
else
    echo "BWA is already installed."
fi

# Map reads using BWA MEM
echo "Mapping reads..."
if bwa mem "$reference" "$input_fastq" > "${output_dir}/aligned_reads.sam"; then
    echo "Mapping complete. SAM file stored in $output_dir"
else
    echo "Error: BWA mapping failed!"
    exit 1
fi
