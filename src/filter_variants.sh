#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_vcf> <output_dir>"
    exit 1
fi

input_vcf=$1
output_dir=$2

# Check if input file exists
if [ ! -f "$input_vcf" ]; then
    echo "Error: Input VCF file '$input_vcf' not found!"
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

# Function to install bcftools using Conda
install_bcftools() {
    echo "Installing bcftools..."
    conda create -n bcftools_env -c bioconda bcftools -y
    conda activate bcftools_env
    echo "bcftools installation complete."
}

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    install_conda
else
    echo "Conda is already installed."
fi

# Activate the Conda environment and install bcftools
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bcftools_env || conda activate
if ! command -v bcftools &> /dev/null; then
    install_bcftools
else
    echo "bcftools is already installed."
fi

# Filter variants using bcftools
echo "Filtering variants..."
if bcftools filter -s LOWQUAL -e '%QUAL<20 || DP<10' "$input_vcf" > "${output_dir}/filtered_variants.vcf"; then
    echo "Filtering complete. Filtered VCF stored in $output_dir"
else
    echo "Error: bcftools filtering failed!"
    exit 1
fi
