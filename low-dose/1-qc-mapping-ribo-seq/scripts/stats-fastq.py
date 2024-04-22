import argparse
import numpy as np
import matplotlib.pyplot as plt
import gzip
import os


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="out", help="file name of result.")
    parser.add_argument("fastq", nargs="+", help="fastq file, support .gz format.")
    args = parser.parse_args()
    return args.out, args.fastq


def stats_fastq_gz(filename):
    len_dic = {}
    fin = gzip.open(filename, "r")
    idx = 0
    for line in fin:
        if idx % 4 == 1:
            line_len = len(line.decode().rstrip())
            if line_len in len_dic:
                len_dic[line_len] += 1
            else:
                len_dic[line_len] = 1
        idx += 1

    fin.close()
    return len_dic


def stats_fastq(filename):
    len_dic = {}
    fin = open(filename, "r")
    idx = 0
    for line in fin:
        if idx % 4 == 1:
            line_len = len(line.rstrip())
            if line_len in len_dic:
                len_dic[line_len] += 1
            else:
                len_dic[line_len] = 1
        idx += 1
    fin.close()
    return len_dic


def main():
    out, fastq = get_args()
    len_dics = []
    for filename in fastq:
        if filename.endswith("gz"):
            len_dics.append(stats_fastq_gz(filename))
        else:
            len_dics.append(stats_fastq(filename))
    
    fig, axs = plt.subplots(nrows=1, ncols=2, layout="constrained", figsize=(8, 4))
    for idx, len_dic in enumerate(len_dics):
        keys = list(len_dic.keys())
        keys.sort()
        len_list = []
        for kk in keys:
            len_list.append(len_dic[kk])
        print("read number", os.path.basename(fastq[idx]), sum(len_list))
        axs[0].plot(keys, len_list, label=os.path.basename(fastq[idx]))
        #axs[0].plot(keys, len_list, label=os.path.basename(fastq[idx]).split(".")[0])
        len_list_pct = np.asarray(len_list) / sum(len_list)
        axs[1].plot(keys, len_list_pct, label=os.path.basename(fastq[idx]))
        #axs[1].plot(keys, len_list_pct, label=os.path.basename(fastq[idx]).split(".")[0])
    axs[0].set_xticks(keys, keys, rotation=90)
    axs[1].set_xticks(keys, keys, rotation=90)
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
