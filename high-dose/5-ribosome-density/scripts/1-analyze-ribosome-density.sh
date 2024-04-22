working_dire="/home/fanghl/work/sdu/lutong/ribo-seq/run-20240411/working_dire"
count_data=${working_dire}/high-dose/3-reads-annotation-and-counting/3-count-mtx
result_data=${working_dire}/high-dose/5-ribosome-density/results
script_dire=${woking_dire}/high-dose/5-ribosome-density/scripts

python3 ${script_dire}/ribosome-density.py \
    --out ${result_data}/low-dose-ribosome-density.tsv \
    ${count_data}/ribo-seq/ribo-count-mtx.tsv \
    ${count_data}/rna-seq/rna-count-mtx.tsv

python3 ${script_dire}/plot-ribosome-density-vs-mRNA-change.py \
    --min_count 64 \
    --trans_factor .5 \
    --marker_gene ${script_dire}/plot_gene_items \
    --out ${result_data}/low-dose-min-count-64-trans-factor-0.5 \
    ${result_data}/low-dose-ribosome-density.tsv
