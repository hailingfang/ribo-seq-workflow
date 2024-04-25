# Ribo-seq Analysis Workflow

This repository contain codes and information to re-present the ribo-seq analysis of the paper.

## Requirements

- Python3 >= 3.11
- biobrary == 0.1.4
    A python3 module, can be installed by `pip install biobrary==0.1.4`
- annoread == v0.1
    The code of this project has been copied to corresponding directories, you do not need to install it sepratly. The code is also avaiable at
    https://github.com/benjaminfang/annoread.
- matplotlib
- numpy
- pandas
- scipy
- cutadapt == 4.6
- bowtie2 == 2.5.2

## Installation

All scripts in this workflow is programmed using Python or Bash. So you can run those script directly after you have install the requirement.

## Usage

1. Download the Ribo-seq and RNA-seq raw data from NCBI as well as the reference data.

2. Place the raw data under corresponding directories.

3. Extract mRNA and non-coding RNA infromation from GCF_000002035.6_GRCz11_rna_from_genomic.fna. And place it under ${working_dire}/ref_data.

4. Make index of reference transcriptome using bowtie2 for all-RNA, mRNA and
non-coding RNA. And place it under ${working_dire}/ref_data.

5. Go to low-dose directory, and change direcory into 1-*/scripts direcotory in order. And run the Bash script.

6. Do same operation for high-dose data.
