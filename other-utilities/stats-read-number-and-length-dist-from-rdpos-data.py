import sys
sys.path.append("/home/fanghl/work/git-repo/annoread")
from annodata import ANNOREAD_DATA
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="out", help="out file name")
    parser.add_argument("rdpos_files", nargs="+", help="ribos_files")
    args = parser.parse_args()
    return args.out, args.rdpos_files


def stats_file_read_len(filename):
    data = ANNOREAD_DATA(filename)
    data = data.get_data()
    len_dic = {}
    read_id_set = set()
    for gene_name in data:
        for gene_id in data[gene_name]:
            for transcript_id in data[gene_name][gene_id]["transcripts"]:
                for read in data[gene_name][gene_id]["transcripts"][transcript_id]["reads"]:
                    read_id = read[0]
                    read_len = read[2]
                    read_id_set.add(read_id)
                    if read_len in len_dic:
                        len_dic[read_len] += 1
                    else:
                        len_dic[read_len] = 1
    keys = list(len_dic.keys())
    keys.sort()
    vales = []
    for kk in keys:
        vales.append(len_dic[kk])

    return keys, vales


def main():
    out, rdpos_files = getargs()
    fig, axs = plt.subplots(nrows=1, ncols=2, layout="constrained", figsize=(8, 4))
    for ff in rdpos_files:
        ff_basename = os.path.basename(ff)
        length, len_count = stats_file_read_len(ff)
        print("read number:", ff_basename, sum(len_count))
        axs[0].plot(length, len_count, label=ff_basename)
        len_count_pct = np.asarray(len_count) / sum(len_count)
        axs[1].plot(length, len_count_pct, label=ff_basename)
    axs[0].set_xticks(length, length, rotation=90)
    axs[1].set_xticks(length, length, rotation=90)
    axs[0].legend()
    axs[1].legend()
    axs[0].grid()
    axs[1].grid()
    axs[0].set_xlabel("read length")
    axs[0].set_ylabel("count")  
    axs[1].set_xlabel("read length")
    axs[1].set_ylabel("percent")
    fig.savefig(out + ".svg")


if __name__ == "__main__":
    main()