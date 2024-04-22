working_dire="/home/fanghl/work/sdu/lutong/ribo-seq/run-20240411/working_dire"
script_dire=${working_dire}/high-dose/6-reads-coverage-of-genes/scripts
read_anno_data=${working_dire}/high-dose/3-reads-annotation-and-counting/2-reads-assigned
result_dire=${working_dire}/high-dose/6-reads-coverage-of-genes/results

for ((i=1; i<3; i++))
do
    python3 ${script_dire}/plot-read-coverage.py \
        --out ${result_dire}/low-dose-dsR-replicate-${i}-${1} \
        ${read_anno_data}/ribo-seq/dsR-${i}-assigned.rdpos \
        ${read_anno_data}/rna-seq/dsR-${i}-assigned.rdpos \
        ${1}

    python3 ${script_dire}/plot-read-coverage.py \
        --out ${result_dire}/low-dose-un-replicate-${i}-${1} \
        ${read_anno_data}/ribo-seq/un-${i}-assigned.rdpos \
        ${read_anno_data}/rna-seq/un-${i}-assigned.rdpos \
        ${1}
done
