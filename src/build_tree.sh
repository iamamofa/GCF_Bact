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

# Function to install IQ-TREE using Conda
install_iqtree() {
    echo "Installing IQ-TREE..."
    conda create -n iqtree_env -c bioconda iqtree -y
    conda activate iqtree_env
    echo "IQ-TREE installation complete."
}

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    install_conda
else
    echo "Conda is already installed."
fi

# Activate the Conda environment and install IQ-TREE
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate iqtree_env || conda activate
if ! command -v iqtree &> /dev/null; then
    install_iqtree
else
    echo "IQ-TREE is already installed."
fi

# Build tree using IQ-TREE
echo "Building phylogenetic tree..."
if iqtree -s "$input_fasta" -o "${output_dir}/tree_output"; then
    echo "Tree building complete. Results stored in $output_dir"
else
    echo "Error: IQ-TREE tree building failed!"
    exit 1
fi
