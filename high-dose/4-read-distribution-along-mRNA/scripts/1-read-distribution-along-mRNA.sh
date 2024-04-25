working_dire="/home/fanghl/work/sdu/lutong/ribo-seq/run-20240411/working_dire"
assigned_read=${working_dire}/high-dose/3-reads-annotation-and-counting/2-reads-assigned
read_pos_res=${working_dire}/high-dose/4-read-distribution-along-mRNA/1-read-pos
script_dire=${working_dire}/high-dose/4-read-distribution-along-mRNA/scripts

for ((i=1; i<3; i++))
do
    python3 ${script_dire}/position-relative-to-start-stop-codon.py \
        --out ${read_pos_res}/ribo-seq/dsR-${i}-codon \
        ${assigned_read}/ribo-seq/dsR-${i}-assigned.rdpos

    python3 ${script_dire}/position-relative-to-start-stop-codon.py \
        --out ${read_pos_res}/ribo-seq/un-${i}-codon \
        ${assigned_read}/ribo-seq/un-${i}-assigned.rdpos
done

for ((i=1; i<4; i++))
do
    python3 ${script_dire}/position-relative-to-start-stop-codon.py \
        --out ${read_pos_res}/rna-seq/dsR-${i}-codon \
        ${assigned_read}/rna-seq/dsR-${i}-assigned.rdpos

    python3 ${script_dire}/position-relative-to-start-stop-codon.py \
        --out ${read_pos_res}/rna-seq/un-${i}-codon \
        ${assigned_read}/rna-seq/un-${i}-assigned.rdpos
done

python3 ${script_dire}/plot-read-postion-and-distribution.py \
    --out high-dose \
    --ribo_data_tr ${read_pos_res}/ribo-seq/dsR-1-codon.rdpos \
        ${read_pos_res}/ribo-seq/dsR-2-codon.rdpos \
    --ribo_data_ut ${read_pos_res}/ribo-seq/un-1-codon.rdpos \
        ${read_pos_res}/ribo-seq/un-2-codon.rdpos \
    --rna_data_tr ${read_pos_res}/rna-seq/dsR-1-codon.rdpos \
        ${read_pos_res}/rna-seq/dsR-2-codon.rdpos \
        ${read_pos_res}/rna-seq/dsR-3-codon.rdpos \
    --rna_data_ut ${read_pos_res}/rna-seq/un-1-codon.rdpos \
        ${read_pos_res}/rna-seq/un-2-codon.rdpos \
        ${read_pos_res}/rna-seq/un-3-codon.rdpos 
