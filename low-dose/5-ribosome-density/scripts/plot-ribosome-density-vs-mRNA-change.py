import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("ribodensity_data", help="ribosome density data file")
    parser.add_argument("--min_count", default=100, type=int, help="the minium count of ribo-seq or rna-seq data.")
    parser.add_argument("--trans_factor", default=1, type=float, help="a factor scale total translation efficency")
    parser.add_argument("--marker_gene", default=None, help="a list of gene which would be marked on plot")
    parser.add_argument("--out", default="out", help="out file name")

    args = parser.parse_args()
    return (args.ribodensity_data, args.min_count, args.trans_factor,
            args.marker_gene, args.out)


def read_data(filename):
    data = pd.read_csv(filename, header=0, index_col=0, sep="\t")
    header = np.asarray(data.columns)
    gene_id_s = np.asarray(data.index)
    data = data.to_numpy()
    return header, gene_id_s, data


def filter_data(data, min_count):
    keep_bool = []
    for idx, row in enumerate(data):
        ribo_tr = row[1:3]
        ribo_ut = row[3:5]
        rna_tr = row[5:8]
        rna_ut = row[8:11]
        if np.mean(ribo_tr) > min_count or np.mean(ribo_ut) > min_count or \
            np.mean(rna_tr) > min_count or np.mean(rna_ut) > min_count:
            keep_bool.append(idx)
    return keep_bool


def check_zero_means(data, gene_id_s, out):
    """
    check 0 in range(20, 24), if present, using 1 to replace it.
    """
    fout = open(out + '-zero-means.tsv', "w")
    for ii, row in enumerate(data):
        flg = 0
        for idx in range(20, 24):
            if row[idx] == 0:
                flg = 1
                row[idx] = 1
        if flg:
            print("\t".join([gene_id_s[ii]] + [str(ele) for ele in row]), file=fout)
        row[24] = row[20] / row[21]
        row[25] = row[22] / row[23]
        row[26] = row[20] / row[22]
        row[27] = row[21] / row[23]
        row[28] = row[26] / row[27]
    fout.close()


def multi_trans_factor(data, trans_factor):
    new_ribo_density_change_tr_ut = []
    for row in data:
        ribo_change_factored = row[24] * trans_factor
        ribo_density_change_factored = row[-1] * trans_factor
        log2_fc_ribo = np.log2(ribo_density_change_factored)
        rna_change = row[25]
        log2_fc_rna = np.log2(rna_change)

        new_ribo_density_change_tr_ut.append([ribo_change_factored,
                                              ribo_density_change_factored,
                                              log2_fc_ribo, log2_fc_rna])
    new_ribo_density_change_tr_ut = np.asarray(new_ribo_density_change_tr_ut)
    data = np.hstack((data, new_ribo_density_change_tr_ut))
    return data


def output_filter_data(header, gene_id_s, data, out):
    fout = open(out + "-filter.tsv", "w")
    print('\t'.join(["gene_name"] + list(header) + ["ribo_change_tr_ut_factored", "ribo_density_change_tr_ut_factored", \
          "log2fc_ribo_density", "log2fc_rna_change"]), file=fout)
    for idx, gene_id in enumerate(gene_id_s):
        print("\t".join([gene_id] + [str(ele) for ele in data[idx]]), file=fout)
    fout.close()


def plot_scatter(gene_id_s, data, marker_gene_file, out):
    fig, ax = plt.subplots(layout="constrained", figsize=(5, 5))

    if marker_gene_file:
        marker_gene = [line.rstrip() for line in open(marker_gene_file)]
    else:
        marker_gene = []

    x_1 = []
    y_1 = []
    x_2 = []
    y_2 = []
    x_3 = []
    y_3 = []

    marker_gene_dic = {}
    ax.axvline(0, color="black", alpha=.5, ls="--", lw=1)
    ax.axhline(0, color="black", alpha=.5, ls="--", lw=1)
    for idx, row in enumerate(data):
        gene_id = gene_id_s[idx]
        ribo_mean_change = row[-4]
        rna_fc = row[-1]
        ribo_fc = row[-2]
        if gene_id in marker_gene:
            x_1.append(rna_fc)
            y_1.append(ribo_fc)
            marker_gene_dic[gene_id] = (rna_fc, ribo_fc)
        elif ribo_mean_change > 1.0:
            x_2.append(rna_fc)
            y_2.append(ribo_fc)
        else:
            x_3.append(rna_fc)
            y_3.append(ribo_fc)
    print(marker_gene_dic)
    ax.scatter(x_3, y_3, s=8, color="gray", alpha=.8, linewidths=0)
    ax.scatter(x_2, y_2, s=8, color="#ffaaaa", linewidths=0)
    ax.scatter(x_1, y_1, s=8, color="red", linewidths=0)
    for gene in marker_gene_dic:
        ax.annotate(gene, np.asarray(marker_gene_dic[gene]) + .03,
                np.asarray(marker_gene_dic[gene]) + .18,
                arrowprops=dict(facecolor='black', headlength=4, headwidth=4, width=1))
    ax.set_title("mRNA abundance and translational efficiency (ds/un)")
    ax.set_xlabel("mRNA $log_2$ fold change")
    ax.set_ylabel("ribosome density $log_2$ fold change")
    #ax.set_xlim([-2, 3])

    fig.savefig(out + ".svg")


def main():
    ribodensity_data, min_count, trans_factor, marker_gene, out = getargs()
    header, gene_id_s, data = read_data(ribodensity_data)
    keep_bool = filter_data(data, min_count)
    gene_id_s = gene_id_s[keep_bool]
    data = data[keep_bool]
    check_zero_means(data, gene_id_s, out)
    data = multi_trans_factor(data, trans_factor)
    output_filter_data(header, gene_id_s, data, out)

    plot_scatter(gene_id_s, data, marker_gene, out)


if __name__ == "__main__":
    main()