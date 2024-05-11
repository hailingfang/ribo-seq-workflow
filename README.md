# The Ribo-seq Analysis Workflow

This repository contains code and information to reproduce
the Ribo-seq analysis in our *paper* currently
under submission.

For any inquiries or assistance needed, please don't hesitate
to contact us at shaoming@sdu.edu.cn.

## Requirements

- Python3 >= 3.11

- biobrary == 0.1.5

    A Python3 module, can be installed using `pip install biobrary==0.1.5`.

- annoread == 0.1.5

    The code for this utility has already been copied to the corresponding
    directories of this workflow, so you do not need to install it separately. 
    The code is also available at https://github.com/benjaminfang/annoread.

- matplotlib

- numpy

- pandas

- cutadapt == 4.6

- bowtie2 == 2.5.2

- samtools == 1.18

## Installation

All scripts in this workflow are programmed using Python or Bash,
so you can run those scripts directly after you have
installed the requirements.
The version number listed in Requirements is the software
version used in our workflow, but other versions may work as well.

## Usage

1. Download the Ribo-seq and RNA-seq raw data from NCBI,
as well as the transcriptomic reference data and the GTF data.

2. Place the raw data under the corresponding 'raw_data' 
directories in the 'low-dose' or 'high-dose' directories, respectively.

3. Extract mRNA and non-coding RNA sequences from the transcriptomic
reference data and place them under "${working_dir}/ref_data".
Then, use Bowtie2 to index the reference data. Replace the
existing directories/files under "${working_dir}/ref_data" 
with the newly indexed directories/files.

4. Navigate to the 'low-dose' directory, then change into each 
subdirectory (1-qc-mapping-ribo-seq, 2-qc-mapping-rna-seq, 
3-reads-annotation-and-counting, 4-read-distribution-along-mRNA, 
5-ribosome-density, and 6-reads-coverage-of-genes) in order, 
and run the Bash script under the 'scripts' directory. 
You may need to replace the variable value of 'working_dire' 
with your actual working directory path in the Bash scripts 
before running them.

5. Perform the same operation for the high-dose data.

## Citations

Tong Lu. et al. Double-stranded RNA triggers a distinct innate immune
response in the early embryo. (under submission)
