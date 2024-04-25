import argparse
import matplotlib.pyplot as plt
import numpy as np
import os


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="out", help="out file name.")
    parser.add_argument("rbpos", nargs="+", help="ribosome position file")
    args = parser.parse_args()
    return args.out, args.rbpos


def stats_read_len(rbpos):
    readlen_dic = {}
    read_set = set()
    fin = open(rbpos, "r")
    for line in fin:
        if line[0] != ">" and line[0] != "$":
            line = line.rstrip().split()
            readid = line[0]
            read_len = int(line[2])
            if readid not in read_set:
                read_set.add(readid)
                if read_len in readlen_dic:
                    readlen_dic[read_len] += 1
                else:
                    readlen_dic[read_len] = 1
    fin.close()
    return readlen_dic


def main():
    out, rbpos_files = getargs()
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 4), layout="constrained")
    for rbpos in rbpos_files:
        readlen_dic = stats_read_len(rbpos)
        keys = list(readlen_dic.keys())
        keys.sort()
        len_list = []
        for kk in keys:
            len_list.append(readlen_dic[kk])
        axs[0].plot(keys, len_list, label=os.path.basename(rbpos))
        len_list_pct = np.asarray(len_list) / sum(len_list)
        axs[1].plot(keys, len_list_pct, label=os.path.basename(rbpos))
    axs[0].set_xticks(keys, keys, rotation=90)
    axs[1].set_xticks(keys, keys, rotation=90)
    axs[0].grid()
    axs[1].grid()
    axs[0].legend()
    axs[1].legend()
    axs[0].set_xlabel("read length")
    axs[0].set_ylabel("count")  
    axs[1].set_xlabel("read length")
    axs[1].set_ylabel("percent")
    fig.savefig(out + ".svg")


if __name__ == "__main__":
    main()