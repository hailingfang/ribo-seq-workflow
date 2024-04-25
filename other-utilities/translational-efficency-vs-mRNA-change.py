import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="out", help="out file name")
    parser.add_argument("ribocount", help="file name ribo count")
    parser.add_argument("rnacount", help="file name of RNA-seq count data")
    parser.add_argument("--min_count", default=100, type=int, help="filtering threthold to filter Ribosome and RNA count data")
    parser.add_argument("--trans_factor", type=float, help="a factor scale total translation efficency")
    parser.add_argument("--plot_gene_id", help="gene id list which will shown in plot.")
    args = parser.parse_args()
    return args.out, args.ribocount, args.rnacount, args.min_count, args.trans_factor, args.plot_gene_id


def filter_and_norm_data(data, min_count):
    f_bool = []
    row_len = len(data[0]) // 2
    for row in data:
        if np.mean(row[:row_len]) > min_count or np.mean(row[row_len:]) > min_count:
            f_bool.append(True)
        else:
            f_bool.append(False)

    data_sum = np.sum(data, axis=0)
    data *= 1e6 / data_sum

    return f_bool, data


def get_intersection(f_bool1, ribodt_geneid, f_bool2, rnadt_geneid):
    comm_gene = set(ribodt_geneid[f_bool1]) & set(rnadt_geneid[f_bool2]) 
    comm_gene = list(comm_gene)
    bool1 = []
    bool2 = []
    ribodt_geneid = list(ribodt_geneid)
    rnadt_geneid = list(rnadt_geneid)
    for gg in comm_gene:
        bool1.append(ribodt_geneid.index(gg))
        bool2.append(rnadt_geneid.index(gg))

    return bool1, bool2, comm_gene


def cal_density(ribodt, rnadt, trans_factor):
    ribo_mean = []
    rna_mean = []
    ribo_row_len = len(ribodt[0]) // 2
    rna_row_len = len(rnadt[0]) // 2
    
    for row in ribodt:
        ribo_mean.append([np.mean(row[: ribo_row_len]), \
                        np.mean(row[ribo_row_len:])])
    for row in rnadt:
        rna_mean.append([np.mean(row[: rna_row_len]), \
                        np.mean(row[rna_row_len: ])])

    density = []

    for idx, row in enumerate(ribo_mean):
        row2 = rna_mean[idx]
        density.append([row[0]/row2[0], row[1]/row2[1]])
    print(density)
    rna_mean = np.asarray(rna_mean)
    density = np.asarray(density)
    density[:, 0] = density[:, 0] * trans_factor
    rna_fold_change = rna_mean[:, 0] / rna_mean[:, 1]
    density_fold_change = density[:, 0] / density[:, 1]

    return rna_fold_change, density_fold_change


def plot_scatter(rna_fold_change, density_fold_change, comm_gene, gene_plot_list, out):
    fig, ax = plt.subplots(layout="constrained")
    rna_change = np.log2(rna_fold_change)
    density_change = np.log2(density_fold_change)

    rna_change_list = list(rna_change)
    density_change_list = list(density_change)
    rna_change_sort = copy.deepcopy(rna_change_list)
    density_change_sort = copy.deepcopy(density_change_list)
    rna_change_sort.sort()
    density_change_sort.sort()
    
    boundary_point = []
    boundary_point.append(rna_change_list.index(rna_change_sort[0]))
    boundary_point.append(rna_change_list.index(rna_change_sort[1]))
    boundary_point.append(rna_change_list.index(rna_change_sort[-1]))
    boundary_point.append(rna_change_list.index(rna_change_sort[-2]))
    boundary_point.append(density_change_list.index(density_change_sort[0]))
    boundary_point.append(density_change_list.index(density_change_sort[1]))
    boundary_point.append(density_change_list.index(density_change_sort[-1]))
    boundary_point.append(density_change_list.index(density_change_sort[-2]))

    direc1 = np.array([1,1])
    direc2 = np.array([-1, 1])
    direc1_v = direc1 @ np.vstack((rna_change, density_change))
    direc2_v = direc2 @ np.vstack((rna_change, density_change))
    
    direc1_v_list = list(direc1_v)
    direc1_v_list.sort()
    boundary_point.append(list(direc1_v).index(direc1_v_list[0]))
    boundary_point.append(list(direc1_v).index(direc1_v_list[1]))
    boundary_point.append(list(direc1_v).index(direc1_v_list[-1]))
    boundary_point.append(list(direc1_v).index(direc1_v_list[-2]))

    direc2_v_list = list(direc2_v)
    direc2_v_list.sort()
    boundary_point.append(list(direc2_v).index(direc2_v_list[0]))
    boundary_point.append(list(direc2_v).index(direc2_v_list[1]))
    boundary_point.append(list(direc2_v).index(direc2_v_list[-1]))
    boundary_point.append(list(direc2_v).index(direc2_v_list[-2]))

    for pp in boundary_point:
        gene_plot_list.append(comm_gene[pp])
    gene_plot_list = list(set(gene_plot_list))
    print(gene_plot_list)
    
    point1x = []
    point1y = []
    point2x = []
    point2y = []
    plot_loc_dic = {}
    for idx, gg in enumerate(comm_gene):
        if gg in gene_plot_list:
            point1x.append(rna_change[idx])
            point1y.append(density_change[idx])
            plot_loc_dic[gg] = [rna_change[idx], density_change[idx]]
        else:
            point2x.append(rna_change[idx])
            point2y.append(density_change[idx])

    ax.scatter(point2x, point2y, s=2)
    ax.scatter(point1x, point1y, s=2, c="red")
    ax.set_title("mRNA abundance and translational efficiency (ds/un)")
    ax.set_xlabel("mRNA $log_2$ fold change")
    ax.set_ylabel("ribosome dnesity $log_2$ fold change")
    ax.axvline(0, color="black", alpha=.5, ls="--")
    ax.axhline(0, color="black", alpha=.5, ls="--")
    for gg in plot_loc_dic:
        ax.annotate(gg, plot_loc_dic[gg], np.asarray(plot_loc_dic[gg]) + .1,\
            arrowprops=dict(facecolor='black', headlength=4, headwidth=6, shrink=.9,  width=2))
    fig.savefig(out + "-fold-change.svg")
    return


def print_density_data(comm_gene, ribodt_head, rnadt_head, ribodt, rnadt, out, trans_factor):
    fout = open(out + "-ribo-density.tsv", "w")
    print('\t'.join(["geneid"] + ["ribo_" + e for e in ribodt_head]+ ["rna_" + e for e in rnadt_head] + 
           ["ribo_treatment_ave", "ribo_control_ave", "rna_treatment_ave",
             "rna_control_ave", "treatment_density", "control_density", "density-log2-fold-change", "rna-log2-fold-chage"]), file=fout)
    ribo_len = len(ribodt[0]) // 2
    rna_len = len(rnadt[0]) // 2
    print(ribo_len, rna_len)
    for idx, geneid in enumerate(comm_gene):
        ribo_a, ribo_b = np.mean(ribodt[idx][:ribo_len]), np.mean(ribodt[idx][ribo_len:])
        rna_a, rna_b = np.mean(rnadt[idx][:rna_len]), np.mean(rnadt[idx][rna_len:])
        density_a = ribo_a / rna_a * trans_factor
        density_b = ribo_b / rna_b
        print("\t".join([geneid] + [str(e) for e in ribodt[idx]] + 
                        [str(e) for e in rnadt[idx]] + 
                        [format(ribo_a, ".6f"), format(ribo_b, ".6f"), format(rna_a, ".6f"), format(rna_b, ".6f")] +
                        [format(density_a, ".6f"), format(density_b, ".6f"), format(np.log2(density_a / density_b), ".6f"), format(np.log2(rna_a / rna_b), ".6f")]), file=fout)
    fout.close()


def print_not_filtered_data(ribodt_geneid, ribodt_head, ribodt, rnadt_geneid, rnadt_head, rnadt, out, trans_factor):
    ribodt = np.copy(ribodt)
    rnadt = np.copy(rnadt)
    ribo_dt_sum = np.sum(ribodt, axis=0)
    ribodt *= 1e6 / ribo_dt_sum
    rnadt_dt_sum = np.sum(rnadt, axis=0)
    rnadt *= 1e6 / rnadt_dt_sum

    common_gene = set(ribodt_geneid) & set(rnadt_geneid)
    common_gene = list(common_gene)
    ribo_bool = []
    rna_bool = []
    ribodt_geneid = list(ribodt_geneid)
    rnadt_geneid = list(rnadt_geneid)
    for gid in common_gene:
        ribo_bool.append(ribodt_geneid.index(gid))
        rna_bool.append(rnadt_geneid.index(gid))

    ribodt = ribodt[ribo_bool]
    rnadt = rnadt[rna_bool]
    print(ribodt.shape)
    print(rnadt.shape)
    ribo_len = len(ribodt[0]) // 2
    rna_len = len(rnadt[0]) // 2
    ribo_mean_treat = np.mean(ribodt[:, : ribo_len], axis=1)
    ribo_mean_contro = np.mean(ribodt[:, ribo_len:], axis=1)
    rna_mean_treat = np.mean(rnadt[:, : rna_len], axis=1)
    rna_mean_contro = np.mean(rnadt[:, rna_len:], axis=1)
    print(ribo_mean_treat.shape)

    fout = open(out + "-before-filter.tsv", "w")
    print('\t'.join(["geneid"] + ["ribo_" + e for e in ribodt_head]+ ["rna_" + e for e in rnadt_head] + 
           ["ribo_treatment_ave", "ribo_control_ave", "rna_treatment_ave",
             "rna_control_ave", "treatment_density", "control_density", "density-log2-fold-change", "rna-log2-fold-chage"]), file=fout)
    print(len(common_gene))
    for idx, gid in enumerate(common_gene):
        print(idx, gid)
        ribo_mean_t = ribo_mean_treat[idx]
        ribo_mean_c = ribo_mean_contro[idx]
        rna_mean_t = rna_mean_treat[idx]
        rna_mean_c = rna_mean_contro[idx]
        if rna_mean_t >0 and rna_mean_c > 0 and ribo_mean_t > 0 and ribo_mean_c > 0:
            density_t = ribo_mean_t / rna_mean_t * trans_factor
            density_c = ribo_mean_c / rna_mean_c
            density_lfc = np.log2(density_t / density_c)
            rna_lfc = np.log2(rna_mean_t / rna_mean_c)
            print("\t".join([gid] + [str(e) for e in ribodt[idx]] + [str(e) for e in rnadt[idx]]  + [str(ribo_mean_t), str(ribo_mean_c), str(rna_mean_t), str(rna_mean_c)] + [format(density_t, ".6f"), format(density_c, ".6f"), format(density_lfc, ".6f"), format(rna_lfc, ".6f")]), file=fout)
        else:
            density_t = "NA"
            density_c = "NA"
            density_lfc = "NA"
            rna_lfc = "NA"
            print("\t".join([gid] + [str(e) for e in ribodt[idx]] + [str(e) for e in rnadt[idx]] + [str(ribo_mean_t), str(ribo_mean_c), str(rna_mean_t), str(rna_mean_c)] + [density_t, density_c, density_lfc, rna_lfc]), file=fout)


def main():
    out, ribocount, rnacount, min_count, trans_factor, plot_gene_id = getargs()
    ribodt = pd.read_csv(ribocount, header=0, index_col="geneid", sep="\t")
    rnadt = pd.read_csv(rnacount, header=0, index_col="geneid", sep="\t")

    ribodt_head = ribodt.columns
    ribodt_geneid = np.asarray(ribodt.index)
    ribodt = ribodt.to_numpy(dtype="float32")

    rnadt_head = rnadt.columns
    rnadt_geneid = np.asarray(rnadt.index)
    rnadt = rnadt.to_numpy(dtype="float32")

    print_not_filtered_data(ribodt_geneid, ribodt_head, ribodt, rnadt_geneid, rnadt_head, rnadt, out, trans_factor)

    f_bool1, ribodt = filter_and_norm_data(ribodt, min_count)
    f_bool2, rnadt = filter_and_norm_data(rnadt, min_count)

    boo1, boo2, comm_gene = get_intersection(f_bool1, ribodt_geneid, f_bool2, rnadt_geneid)

    ribodt = ribodt[boo1]
    rnadt = rnadt[boo2]

    rna_fold_change, density_fold_change = cal_density(ribodt, rnadt, trans_factor)

    print_density_data(comm_gene, ribodt_head, rnadt_head, ribodt, rnadt, out, trans_factor)
    gene_plot_list = [l.rstrip() for l in open(plot_gene_id)]
    plot_scatter(rna_fold_change, density_fold_change, comm_gene, gene_plot_list, out)



if __name__ == "__main__":
    main()