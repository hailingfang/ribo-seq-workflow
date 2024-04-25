working_dire=/home/fanghl/work/sdu/lutong/ribo-seq/run-20240411/working_dire
script_dire=${working_dire}/high-dose/2-qc-mapping-rna-seq/scripts
raw_data=${working_dire}/high-dose/raw-data/rna-seq
clean_data=${working_dire}/high-dose/2-qc-mapping-rna-seq/1-clean-data
mapped_data=${working_dire}/high-dose/2-qc-mapping-rna-seq/2-mapped-data
ref_dire=${working_dire}/ref-data

#remove adapters
for((i=1; i<4; i++))
do
    cutadapt \
        -a AGATCGGAAGAGCACACG \
        -A AGATCGGAAGAGCGTCGT \
        --minimum-length 50 \
        -q 15,10 \
        --trim-n \
        -o ${clean_data}/dsR_${i}-clean-R1.fq.gz \
        -p ${clean_data}/dsR_${i}-clean-R2.fq.gz \
        ${raw_data}/dsR${i}_L1_1.fq.gz \
        ${raw_data}/dsR${i}_L1_2.fq.gz

    cutadapt \
        -a AGATCGGAAGAGCACACG \
        -A AGATCGGAAGAGCGTCGT \
        --minimum-length 50 \
        -q 15,10 \
        --trim-n \
        -o ${clean_data}/un_${i}-clean-R1.fq.gz \
        -p ${clean_data}/un_${i}-clean-R2.fq.gz \
        ${raw_dire}/un${i}_L1_1.fq.gz \
        ${raw_dire}/un${i}_L1_2.fq.gz
done

#map
for ((i=1; i<4; i++))
do
    bowtie2 \
        -x ${ref_data}/rna-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic \
        -1 ${clean_data}/dsR_${i}-clean-R1.fq.gz \
        -2 ${clean_data}/dsR_${i}-clean-R2.fq.gz \
        -k 3 \
        --threads 8 \
        -S ${mapped_data}/dsR-${i}.sam

    bowtie2 \
        -x ${ref_data}/rna-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic \
        -1 ${clean_data}/un_${i}-clean-R1.fq.gz \
        -2 ${clean_data}/un_${i}-clean-R2.fq.gz \
        -k 3 \
        --threads 8 \
        -S ${mapped_data}/un-${i}.sam
done

for ((i=1; i<4; i++))
do
    samtools fixmate \
        --threads 8 \
        -m \
        -O bam \
        ${mapped_data}/dsR-${i}.sam \
        ${mapped_data}/dsR-${i}-fixmate.bam

    samtools sort \
        --threads 8 \
        -O bam \
        -o ${mapped_data}/dsR-${i}-fixmate-sorted.bam \
        ${mapped_data}/dsR-${i}-fixmate.bam

    samtools markdup \
        -O bam \
        --threads 8 \
        -r \
        ${mapped_data}/dsR-${i}-fixmate-sorted.bam \
        ${mapped_data}/dsR-${i}-fixmate-sorted-rmdup.bam

    samtools sort \
    --threads 8 \
    -n \
    -O sam \
    -o ${mapped_data}/dsR-${i}-fixmate-sorted-rmdup-sorted.sam \
    ${mapped_data}/dsR-${i}-fixmate-sorted-rmdup.bam

    samtools fixmate \
        --threads 8 \
        -m \
        -O bam \
        ${mapped_data}/un-${i}.sam \
        ${mapped_data}/un-${i}-fixmate.bam

    samtools sort \
        --threads 8 \
        -O bam \
        -o ${mapped_data}/un-${i}-fixmate-sorted.bam \
        ${mapped_data}/un-${i}-fixmate.bam

    samtools markdup \
        -O bam \
        --threads 8 \
        -r \
        ${mapped_data}/un-${i}-fixmate-sorted.bam \
        ${mapped_data}/un-${i}-fixmate-sorted-rmdup.bam

    samtools sort \
    --threads 8 \
    -n \
    -O sam \
    -o ${mapped_data}/un-${i}-fixmate-sorted-rmdup-sorted.sam \
    ${mapped_data}/un-${i}-fixmate-sorted-rmdup.bam
done


