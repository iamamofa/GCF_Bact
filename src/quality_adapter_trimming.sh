#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fastq> <output_dir>"
    exit 1
fi

input_fastq=$1
output_dir=$2

# Check if input file exists
if [ ! -f "$input_fastq" ]; then
    echo "Error: Input FASTQ file '$input_fastq' not found!"
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

# Function to install Trimmomatic using Conda
install_trimmomatic() {
    echo "Installing Trimmomatic..."
    conda create -n trimmomatic_env -c bioconda trimmomatic -y
    conda activate trimmomatic_env
    echo "Trimmomatic installation complete."
}

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    install_conda
else
    echo "Conda is already installed."
fi

# Activate the Conda environment and install Trimmomatic
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate trimmomatic_env || conda activate
if ! command -v trimmomatic &> /dev/null; then
    install_trimmomatic
else
    echo "Trimmomatic is already installed."
fi

# Define Trimmomatic options and perform trimming
echo "Trimming adapters and performing quality control with Trimmomatic..."
if trimmomatic PE -phred33 "$input_fastq" "${output_dir}/trimmed_output_forward.fastq.gz" "${output_dir}/trimmed_output_reverse.fastq.gz" \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36; then
    echo "Trimming complete. Results stored in $output_dir"
else
    echo "Error: Trimmomatic trimming failed!"
    exit 1
fi
