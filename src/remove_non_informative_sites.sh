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

# Function to install snp-sites using Conda
install_snp_sites() {
    echo "Installing snp-sites..."
    conda create -n snp_sites_env -c bioconda snp-sites -y
    conda activate snp_sites_env
    echo "snp-sites installation complete."
}

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    install_conda
else
    echo "Conda is already installed."
fi

# Activate the Conda environment and install snp-sites
source "$(conda info --base)/etc/profile.d/conda.sh"
if ! conda env list | grep -q 'snp_sites_env'; then
    install_snp_sites
else
    echo "snp-sites environment already exists."
fi

conda activate snp_sites_env

if ! command -v snp-sites &> /dev/null; then
    echo "Error: snp-sites installation failed!"
    exit 1
else
    echo "snp-sites is installed."
fi

# Remove non-informative sites using snp-sites
echo "Removing non-informative sites..."
if snp-sites -o "${output_dir}/informative_sites.fasta" "$input_fasta"; then
    echo "Non-informative sites removed. Result stored in $output_dir"
else
    echo "Error: snp-sites processing failed!"
    exit 1
fi
