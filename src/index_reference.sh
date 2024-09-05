#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fasta> <output_dir>"
    exit 1
fi

input_fasta=$1
output_dir=$2

# Check if input file exists
if [ ! -f "$input_fasta" ]; then
    echo "Error: Input FASTA file '$input_fasta' not found!"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Define the prefix for the indexed files
index_prefix="$output_dir/indexed_reference"

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

# Function to install Bowtie2 using Conda
install_bowtie2() {
    echo "Installing Bowtie2..."
    conda create -n bowtie2_env -c bioconda bowtie2 -y
    conda activate bowtie2_env
    echo "Bowtie2 installation complete."
}

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    install_conda
else
    echo "Conda is already installed."
fi

# Activate the Conda environment and install Bowtie2
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bowtie2_env || conda activate
if ! command -v bowtie2-build &> /dev/null; then
    install_bowtie2
else
    echo "Bowtie2 is already installed."
fi

# Index the reference with Bowtie2
echo "Indexing reference..."
if bowtie2-build "$input_fasta" "$index_prefix"; then
    echo "Indexing complete. Results stored in $output_dir"
else
    echo "Error: Bowtie2 indexing failed!"
    exit 1
fi
