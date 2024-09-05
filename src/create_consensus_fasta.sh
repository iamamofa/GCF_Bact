#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_vcf> <reference> <output_dir>"
    exit 1
fi

input_vcf=$1
reference=$2
output_dir=$3

# Check if input file and reference file exist
if [ ! -f "$input_vcf" ]; then
    echo "Error: Input VCF file '$input_vcf' not found!"
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

# Function to install Python and dependencies using Conda
install_python_and_dependencies() {
    echo "Installing Python and dependencies..."
    conda create -n vcf2pseudogenome_env python=3.8 -y
    conda activate vcf2pseudogenome_env
    pip install -r requirements.txt
    echo "Python and dependencies installation complete."
}

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    install_conda
else
    echo "Conda is already installed."
fi

# Activate the Conda environment and install Python and dependencies
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate vcf2pseudogenome_env || conda activate
if ! command -v python &> /dev/null; then
    install_python_and_dependencies
else
    echo "Python is already installed."
fi

# Create consensus FASTA using vcf2pseudogenome
echo "Creating consensus FASTA..."
if python vcf2pseudogenome.py "$input_vcf" "$reference" > "${output_dir}/consensus.fasta"; then
    echo "Consensus FASTA creation complete. Result stored in $output_dir"
else
    echo "Error: vcf2pseudogenome.py execution failed!"
    exit 1
fi
