import argparse
import numpy as np
import pandas as pd


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("ribo_count", help="ribo-seq count data")
    parser.add_argument("rna_count", help="rna-seq count data")
    parser.add_argument("--out", default="out", help="out file name")

    args = parser.parse_args()
    return args.ribo_count, args.rna_count, args.out


def read_data(filename):
    dt = pd.read_csv(filename, header=0, index_col=0, sep="\t")
    header = np.asarray(dt.columns)
    gene_id_s = np.asarray(dt.index)
    data = dt.to_numpy(dtype=np.int32)
    return header, gene_id_s, data


def norm_data(data):
    data_sum = np.sum(data, axis=0)
    data_out = data * 1e6
    data_out /= data_sum
    return data_out


def get_common_gene_s(gene_1, gene_2):
    comm_gene = set(gene_1) & set(gene_2)
    comm_gene = list(comm_gene)
    comm_gene.sort()
    gene_1_idx = []
    gene_1 = list(gene_1)
    gene_2 = list(gene_2)
    for gene in comm_gene:
        gene_1_idx.append(gene_1.index(gene))
    gene_2_idx = []
    for gene in comm_gene:
        gene_2_idx.append(gene_2.index(gene))

    return comm_gene, gene_1_idx, gene_2_idx


def cal_mean_by_condition(data, c1, c2):
    c1_mean = np.mean(data[:, :c1], axis=1)
    c2_mean = np.mean(data[:, c1: c1 + c2], axis=1)
    return np.vstack((c1_mean, c2_mean)).T


def output_res(ribo_header, rna_header, comm_gene, ribo_data, rna_data,
               ribo_data_norm, rna_data_norm,
               ribo_data_norm_mean, rna_data_norm_mean,
               ribo_change, rna_change,
               ribo_density, ribo_density_change,
               out):
    
    fout = open(out, "w")
    ribo_header = ["ribo_" + ele for ele in ribo_header]
    rna_header = ["rna_" + ele for ele in rna_header]
    ribo_header_norm = [ele + "_norm" for ele in ribo_header]
    rna_header_norm = [ele + "_norm" for ele in rna_header]
    ribo_rna_mean_header = ["ribo_tr_mean", "ribo_ut_mean", "rna_tr_mean", "rna_ut_mean"]
    ribo_rna_change_header = ["ribo_change_tr_ut", "rna_change_tr_ut"]
    ribo_density_header = ["ribo_density_tr", "ribo_densit_ut"]
    ribo_density_change_header = ["ribo_density_change_tr_ut"]
    print("\t".join(["gene"] + ribo_header + rna_header + 
                    ribo_header_norm + rna_header_norm +
                    ribo_rna_mean_header + ribo_rna_change_header +
                    ribo_density_header + ribo_density_change_header), file=fout)

    for idx in range(len(comm_gene)):
        ribo_data_row = [str(ele) for ele in ribo_data[idx]]

        rna_data_row = [str(ele) for ele in rna_data[idx]]
        ribo_data_norm_row = [str(ele) for ele in ribo_data_norm[idx]]
        rna_data_norm_row = [str(ele) for ele in rna_data_norm[idx]]

        ribo_rna_mean_row = [str(ele) for ele in ribo_data_norm_mean[idx]] + \
            [str(ele) for ele in rna_data_norm_mean[idx]]

        ribo_rna_change_row = [str(ribo_change[idx]), str(rna_change[idx])]
        ribo_density_row = [str(ele) for ele in ribo_density[idx]]
        ribo_density_change_row = [str(ribo_density_change[idx])]


        row_str = [comm_gene[idx]] + ribo_data_row + rna_data_row + \
                        ribo_data_norm_row + rna_data_norm_row + \
                        ribo_rna_mean_row + ribo_rna_change_row + \
                        ribo_density_row + ribo_density_change_row
        '''
        print(ribo_data_row)
        print(rna_data_row)
        print(ribo_rna_mean_row)
        print(ribo_rna_change_row)
        print(ribo_density_row)
        print(ribo_density_change_row)
        print(row_str)
        '''

        print("\t".join(row_str), file=fout)

    fout.close()



def main():
    ribo_count, rna_count, out = getargs()
    ribo_header, ribo_gene_s, ribo_data = read_data(ribo_count)
    rna_header, rna_gene_s, rna_data = read_data(rna_count)

    ribo_data_norm = norm_data(ribo_data)
    rna_data_norm = norm_data(rna_data)
    comm_gene, ribo_comm_idx, rna_comm_idx = \
        get_common_gene_s(ribo_gene_s, rna_gene_s)

    ribo_data = ribo_data[ribo_comm_idx]
    rna_data = rna_data[rna_comm_idx]
    ribo_data_norm = ribo_data_norm[ribo_comm_idx]
    rna_data_norm = rna_data_norm[rna_comm_idx]

    ribo_data_norm_mean = cal_mean_by_condition(ribo_data_norm, 2, 2)
    rna_data_norm_mean = cal_mean_by_condition(rna_data_norm, 3, 3)


    ribo_change = []
    rna_change = []
    for row in ribo_data_norm_mean:
        if row[1]:
            ribo_change.append(row[0] / row[1])
        else:
            ribo_change.append(np.nan)
    for row in rna_data_norm_mean:
        if row[1]:
            rna_change.append(row[0] / row[1])
        else:
            rna_change.append(np.nan)

    ribo_density = []
    ribo_density_change = []
    for idx in range(len(comm_gene)):
        ribo_mean_tr, ribo_mean_ut = ribo_data_norm_mean[idx]
        rna_mean_tr, rna_mean_ut = rna_data_norm_mean[idx]
        if rna_mean_tr:
            ribo_density_tr = ribo_mean_tr / rna_mean_tr
        else:
            ribo_density_tr = np.nan
        if rna_mean_ut:
            ribo_density_ut = ribo_mean_ut / rna_mean_ut
        else:
            ribo_density_ut = np.nan
        ribo_density.append([ribo_density_tr, ribo_density_ut])

        if ribo_density_ut and (not np.isnan(ribo_density_ut)) and (not np.isnan(ribo_density_tr)):
            ribo_density_change.append(ribo_density_tr / ribo_density_ut)
        else:
            ribo_density_change.append(np.nan)

    output_res(ribo_header, rna_header, comm_gene, ribo_data, rna_data,
               ribo_data_norm, rna_data_norm,
               ribo_data_norm_mean, rna_data_norm_mean,
               ribo_change, rna_change,
               ribo_density, ribo_density_change,
               out)

if __name__ == "__main__":
    main()