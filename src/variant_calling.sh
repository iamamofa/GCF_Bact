#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_bam> <reference> <output_dir>"
    exit 1
fi

input_bam=$1
reference=$2
output_dir=$3

# Check if input file and reference file exist
if [ ! -f "$input_bam" ]; then
    echo "Error: Input BAM file '$input_bam' not found!"
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

# Function to install GATKv4 and bcftools using Conda
install_tools() {
    echo "Installing GATKv4 and bcftools..."
    conda create -n gatk_bcftools_env -c bioconda gatk4 bcftools -y
    conda activate gatk_bcftools_env
    echo "GATKv4 and bcftools installation complete."
}

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    install_conda
else
    echo "Conda is already installed."
fi

# Activate the Conda environment and install GATKv4 and bcftools
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gatk_bcftools_env || conda activate
if ! command -v gatk &> /dev/null || ! command -v bcftools &> /dev/null; then
    install_tools
else
    echo "GATKv4 and/or bcftools are already installed."
fi

# Call variants using GATKv4
echo "Calling variants using GATKv4..."
if gatk HaplotypeCaller -R "$reference" -I "$input_bam" -O "${output_dir}/variants_gatk.vcf"; then
    echo "GATKv4 variant calling complete."
else
    echo "Error: GATKv4 variant calling failed!"
    exit 1
fi

# Call variants using bcftools mpileup
echo "Calling variants using bcftools mpileup..."
if bcftools mpileup -Ou -f "$reference" "$input_bam" | bcftools call -mv -Oz -o "${output_dir}/variants_bcftools.vcf.gz"; then
    echo "bcftools variant calling complete."
else
    echo "Error: bcftools variant calling failed!"
    exit 1
fi

echo "Variant calling complete. VCF files stored in $output_dir"
s