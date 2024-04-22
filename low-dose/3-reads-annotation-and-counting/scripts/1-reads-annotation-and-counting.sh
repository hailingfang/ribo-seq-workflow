working_dire="/home/fanghl/work/sdu/lutong/ribo-seq/run-20240411/working_dire"
script_dire=${working_dire}/low-dose/3-reads-annotation-and-counting/scripts
ref_data=${working_dire}/ref_data
ribo_mapped_data=${working_dire}/low-dose/1-qc-mapping-ribo-seq/2-mapped-data
rna_mapped_data=${working_dire}/low-dose/2-qc-mapping-rna-seq/2-mapped-data
anno_res=${working_dire}/3-reads-annotation-and-counting/1-read-anno
assign_res=${working_dire}/3-reads-annotation-and-counting/2-reads-assigned
count_res=${working_dire}/3-reads-annotation-and-counting/3-count-mtx

#read annotation
for ((i=1; i<3; i++))
do
    python3 ${script_dire}/annoread.py \
        --fasta ${ref_data}/GCF_000002035.6_GRCz11_rna_from_genomic.fna \
        --gtf ${ref_data}/GCF_000002035.6_GRCz11_genomic.gtf \
        --read_type se \
        --out ${anno_res}/ribo-seq/dsR_${i} \
        ${ribo_mapped_data}/dsR_${i}-unmapped-to-non-mRNA-to-mRNA.sam

    python3 ${script_dire}/annoread.py \
        --fasta ${ref_dire}/GCF_000002035.6_GRCz11_rna_from_genomic.fna \
        --gtf ${ref_dire}/GCF_000002035.6_GRCz11_genomic.gtf \
        --read_type se \
        --out ${anno_res}/ribo-seq/un_${i} \
        ${ribo_mapped_data}/un_${i}-unmapped-to-non-mRNA-to-mRNA.sam
done

for ((i=1; i<4; i++))
do
    python3 ${script_dire}/annoread.py \
        --fasta_file ${ref_data}/GCF_000002035.6_GRCz11_rna_from_genomic.fna \
        --gtf_file ${ref_data}/GCF_000002035.6_GRCz11_genomic.gtf \
        --read_type pe \
        --out ${anno_res}/rna-seq/dsR-${i} \
        ${rna_mapppd_data}/dsR-${i}-fixmate-sorted-rmdup-sorted.sam \

    python3 ${script_dire}/annoread.py \
        --fasta_file ${ref_data}/GCF_000002035.6_GRCz11_rna_from_genomic.fna \
        --gtf_file ${ref_data}/GCF_000002035.6_GRCz11_genomic.gtf \
        --read_type pe \
        --out ${anno_res}/rna-seq/un-${i} \
        ${rna_mapped_data}/un-${i}-fixmate-sorted-rmdup-sorted.sam \
done

#assignn read
for ((i=1; i<3; i++))
do
    python3 ${script_dire}/assignread.py \
        --gene_name proportion \
        --gene_id largest \
        --transcript_id largest \
        --unique_transcript_read \
        --out ${assign_res}/ribo-seq/dsR-${i}-assigned \
        ${anno_res}/ribo-seq/dsR_${i}.rdpos

    python3 ${script_dire}/assignread.py \
        --gene_name proportion \
        --gene_id largest \
        --transcript_id largest \
        --unique_transcript_read \
        --out ${assign_res}/ribo-seq/dsR-${i}-assigned \
        ${anno_res}/ribo-seq/dsR_${i}.rdpos
done

for ((i=1; i<4; i++))
do
    python3 ${script_dire}/assignread.py \
        --gene_name proportion \
        --gene_id largest \
        --transcript_id largest \
        --unique_transcript_read \
        --out ${assign_res}/rna-seq/dsR-${i}-assigned \
        ${anno_res}/rna-seq/dsR-${i}.rdpos

    python3 ${script_dire}/assignread.py \
        --gene_name proportion \
        --gene_id largest \
        --transcript_id largest \
        --unique_transcript_read \
        --out ${assign_res}/rna-seq/un-${i}-assigned \
        ${anno_res}/rna-seq/un-${i}.rdpos
done

#counting
python3 ${script_dire}/countmtx.py \
    --count_level gene_name \
    --out ${count_res}/ribo-seq/ribo-count-mtx \
    ${assign_res}/ribo-seq/dsR-{1,2}-assigned.rdpos \
    ${assign_res}/ribo-seq/un-{1,2}-assigned.rdpos


python3 ${script_dire}/countmtx.py \
    --count_level gene_name \
    --out ${count_res}/rna-seq/rna-count-mtx \
    ${assign_res}/rna-seq/dsR-{1,2,3}-assigned.rdpos \
    ${assign_res}/rna-seq/un-{1,2,3}-assigned.rdpos




