working_dire=/home/fanghl/work/sdu/lutong/ribo-seq/run-20240411/working_dire
script_dire=${working_dire}/high-dose/1-qc-mapping-ribo-seq/scripts
raw_data=${working_dire}/high-dose/raw-data/ribo-seq
clean_data=${working_dire}/high-dose/1-qc-mapping-ribo-seq/1-clean-data
mapped_data=${working_dire}/high-dose/1-qc-mapping-ribo-seq/2-mapped-data
ref_dire=${working_dire}/ref-data

#cutadapt
for ((i=1; i<3; i++))
do
    cutadapt \
        -a ACGACGCTCTTCCGATCT...CTGTAGGCACCATCAATAGATCGG \
        -A ATTGATGGTGCCTACAG...AGATCGGAAGAGCGTCG \
        --minimum-length 20 \
        --maximum-length 40 \
        -o ${clean_data}/dsR_${i}-clean.R1.fastq.gz \
        -p ${clean_data}/dsR_${i}-clean.R2.fastq.gz \
        ${raw_data}/dsR_${i}.raw.1.fastq.gz \
        ${raw_data}/dsR_${i}.raw.2.fastq.gz

    cutadapt \
        -a ACGACGCTCTTCCGATCT...CTGTAGGCACCATCAATAGATCGG \
        -A ATTGATGGTGCCTACAG...AGATCGGAAGAGCGTCG \
        --minimum-length 20 \
        --maximum-length 40 \
        -o ${clean_data}/un_${i}-clean.R1.fastq.gz \
        -p ${clean_data}/un_${i}-clean.R2.fastq.gz \
        ${raw_data}/un_${i}.raw.1.fastq.gz \
        ${raw_data}/un_${i}.raw.2.fastq.gz
done

#mapping
#map-to-rRNA and statistic
for ((i=1; i<3; i++))
do
    bowtie2 \
        -x ${ref_dire}/rRNA-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic_rRNA \
        -U ${clean_data}/dsR_${i}-clean.R1.fastq.gz \
        -S ${mapped_data}/dsR_${i}-map-to-rRNA.sam \
        --un ${mapped_data}/dsR_${i}-unmapped-to-rRNA.fastq \
        -k 3 \
        --threads 8

    bowtie2 \
        -x ${ref_dire}/rRNA-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic_rRNA \
        -U ${clean_data}/un_${i}-clean.R1.fastq.gz \
        -S ${mapped_data}/un_${i}-map-to-rRNA.sam \
        --un ${mapped_data}/un_${i}-unmapped-to-rRNA.fastq \
        -k 3 \
        --threads 8
done

python3 ${script_dire}/stats-fastq.py \
    ${mapped_data}/dsR_1-unmapped-to-rRNA.fastq \
    ${mapped_data}/dsR_2-unmapped-to-rRNA.fastq \
    ${mapped_data}/un_1-unmapped-to-rRNA.fastq \
    ${mapped_data}/un_2-unmapped-to-rRNA.fastq

#mapped-to-non-mRNA
for ((i=1; i<3; i++))
do
    bowtie2 \
        -x ${ref_dire}/non-mRNA-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic_non_mrna \
        -U ${clean_data}/dsR_${i}-clean.R1.fastq.gz \
        -S ${mapped_data}/dsR_${i}-map-to-non-mRNA.sam \
        --un ${mapped_data}/dsR_${i}-unmapped-to-non-mRNA.fastq \
        -k 3 \
        --threads 8

    bowtie2 \
        -x ${ref_dire}/non-mRNA-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic_non_mrna \
        -U ${clean_data}/un_${i}-clean.R1.fastq.gz \
        -S ${mapped_data}/un_${i}-map-to-non-mRNA.sam \
        --un ${mapped_data}/un_${i}-unmapped-to-non-mRNA.fastq \
        -k 3 \
        --threads 8

done

python3 ${script_dire}/stats-fastq.py \
    ${mapped_data}/dsR_1-unmapped-to-non-mRNA.fastq \
    ${mapped_data}/dsR_2-unmapped-to-non-mRNA.fastq \
    ${mapped_data}/un_1-unmapped-to-non-mRNA.fastq \
    ${mapped_data}/un_2-unmapped-to-non-mRNA.fastq


#map to mRNA
for ((i=1; i<3; i++))
do
    bowtie2 \
        -x ${ref_dire}/mRNA-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic_mrna \
        -U ${clean_dire}/dsR_${i}-clean.R1.fastq.gz \
        -S ${mapped_dire}/dsR_${i}-map-to-mRNA.sam \
        --un ${mapped_dire}/dsR_${i}-unmapped-to-mRNA.fastq \
        -k 3 \
        --threads 8

    bowtie2 \
        -x ${ref_dire}/mRNA-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic_mrna \
        -U ${clean_data}/un_${i}-clean.R1.fastq.gz \
        -S ${mapped_data}/un_${i}-map-to-mRNA.sam \
        --un ${mapped_data}/un_${i}-unmapped-to-mRNA.fastq \
        -k 3 \
        --threads 8
done

python3 ${script_dire}/stats-fastq.py \
    ${mapped_data}/dsR_1-unmapped-to-mRNA.fastq \
    ${mapped_data}/dsR_2-unmapped-to-mRNA.fastq \
    ${mapped_data}/un_1-unmapped-to-mRNA.fastq \
    ${mapped_data}/un_2-unmapped-to-mRNA.fastq

#map to all RNA
for ((i=1; i<3; i++))
do
    bowtie2 \
        -x ${ref_data}/rna-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic \
        -U ${clean_data}/dsR_${i}-clean.R1.fastq.gz \
        -S ${mapped_data}/dsR_${i}-to-all-RNA.sam \
        --un ${mapped_data}/dsR_${i}-unmapped-to-all-RNA.fastq \
        -k 3 \
        --threads 8

    bowtie2 \
        -x ${ref_data}/rna-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic \
        -U ${clean_data}/un_${i}-clean.R1.fastq.gz \
        -S ${mapped_data}/un_${i}-to-all-RNA.sam \
        --un ${mapped_data}/un_${i}-unmapped-to-all-RNA.fastq \
        -k 3 \
        --threads 8
done

python3 ${script_dire}/stats-fastq.py \
    ${mapped_data}/dsR_1-unmapped-to-all-RNA.fastq \
    ${mapped_data}/dsR_2-unmapped-to-all-RNA.fastq \
    ${mapped_data}/un_1-unmapped-to-all-RNA.fastq \
    ${mapped_data}/un_2-unmapped-to-all-RNA.fastq

#map unmapped to rRNA to mRNA
for ((i=1; i<3; i++))
do
    bowtie2 \
        -x ${ref_data}/mRNA-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic_mrna \
        -U ${mapped_data}/dsR_${i}-unmapped-to-rRNA.fastq \
        -S ${mapped_data}/dsR_${i}-unmapped-to-rRNA-to-mRNA.sam \
        --un ${mapped_data}/dsR_${i}-unmapped-to-rRNA-and-mRNA.fastq \
        -k 3 \
        --threads 8

    bowtie2 \
        -x ${ref_data}/mRNA-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic_mrna \
        -U ${mapped_data}/un_${i}-unmapped-to-rRNA.fastq \
        -S ${mapped_data}/un_${i}-unmapped-to-rRNA-to-mRNA.sam \
        --un ${mapped_data}/un_${i}-unmapped-to-rRNA-and-mRNA.fastq \
        -k 3 \
        --threads 8

done

python3 ${script_dire}/stats-fastq.py \
    ${mapped_data}/dsR_1-unmapped-to-rRNA-and-mRNA.fastq \
    ${mapped_data}/dsR_2-unmapped-to-rRNA-and-mRNA.fastq \
    ${mapped_data}/un_1-unmapped-to-rRNA-and-mRNA.fastq \
    ${mapped_data}/un_2-unmapped-to-rRNA-and-mRNA.fastq

#map unmapped to non-mRNA to mRNA

for ((i=1; i<3; i++))
do
    bowtie2 \
        -x ${ref_data}/mRNA-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic_mrna \
        -U ${mapped_data}/dsR_${i}-unmapped-to-non-mRNA.fastq \
        -S ${mapped_data}/dsR_${i}-unmapped-to-non-mRNA-to-mRNA.sam \
        --un ${mapped_data}/dsR_${i}-unmapped-to-non-mRNA-and-mRNA.fastq \
        -k 3 \
        --threads 8

    bowtie2 \
        -x ${ref_data}/mRNA-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic_mrna \
        -U ${mapped_data}/un_${i}-unmapped-to-non-mRNA.fastq \
        -S ${mapped_data}/un_${i}-unmapped-to-non-mRNA-to-mRNA.sam \
        --un ${mapped_data}/un_${i}-unmapped-to-non-mRNA-and-mRNA.fastq \
        -k 3 \
        --threads 8
done

python3 ${script_dire}/stats-fastq.py \
    ${mapped_data}/dsR_1-unmapped-to-non-mRNA-and-mRNA.fastq \
    ${mapped_data}/dsR_2-unmapped-to-non-mRNA-and-mRNA.fastq \
    ${mapped_data}/un_1-unmapped-to-non-mRNA-and-mRNA.fastq \
    ${mapped_data}/un_2-unmapped-to-non-mRNA-and-mRNA.fastq

#analyse the length of read which mapped to mRNA uniquely.
python3 ${script_dire}/stats-map-to-mRNA-len.py \
    ${mapped_data}/dsR_1-unmapped-to-non-mRNA-to-mRNA.sam,${mapped_data}/dsR_1-unmapped-to-non-mRNA.fastq \
    ${mapped_data}/dsR_2-unmapped-to-non-mRNA-to-mRNA.sam,${mapped_data}/dsR_2-unmapped-to-non-mRNA.fastq \
    ${mapped_data}/un_1-unmapped-to-non-mRNA-to-mRNA.sam,${mapped_data}/un1_1-unmapped-to-non-mRNA.fastq \
    ${mapped_data}/un_2-unmapped-to-non-mRNA-to-mRNA.sam,${mapped_data}/un_2-unmapped-to-non-mRNA.fastq


