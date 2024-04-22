import argparse
import gzip
import matplotlib.pyplot as plt
import numpy as np
import os


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("sam_fastq", help="sam file and fastq file paires. Exp: sam,fastq", nargs="+")
    parser.add_argument("--out", default="out", help="out file name")
    args = parser.parse_args()
    return args.sam_fastq, args.out


def construct_len_dic(fastqfile):
    len_dic = {}
    if fastqfile.endswith("fastq"):
        fin = open(fastqfile, "r")
        idx = 0
        key = None
        for line in fin:
            if idx % 4 == 0:
                key = line.split()[0][1:]
            elif idx % 4 == 1:
                length = len(line.rstrip())
                len_dic[key] = length
            idx += 1
        fin.close()
    elif fastqfile.endswith(".gz"):
        fin = gzip.open(fastqfile, "r")
        idx = 0
        key = None
        for line in fin:
            if idx % 4 == 0:
                key = line.decode().split()[0][1:]
            elif idx % 4 == 1:
                length = len(line.decode().rstrip())
                len_dic[key] = length
            idx += 1
        fin.close()
    else:
        print("file extension not be regnized.")
    return len_dic


def choose_unique_to_mRNA_read(samfile):
    read_dic = {}
    fin = open(samfile, "r")
    
    while True:
        line = fin.readline()
        if line[0] != "@":
            break
    
    line = line.rstrip().split()
    readid, flg, refseq = line[0], line[1], line[2]
    if flg == "0" or flg == "256":
        ref_type = 1 if refseq.split("_")[2] == "mrna" else 0
        if readid not in read_dic:
            read_dic[readid] = [ref_type]
        else:
            read_dic[readid].append(ref_type)
    for line in fin:
        line = line.rstrip().split()
        readid, flg, refseq = line[0], line[1], line[2]
        if flg == "0" or flg == "256":
            ref_type = 1 if refseq.split("_")[2] == "mrna" else 0
            if readid not in read_dic:
                read_dic[readid] = [ref_type]
            else:
                read_dic[readid].append(ref_type)
    fin.close()

    read_list = []
    for readid in read_dic:
        if all(read_dic[readid]):
            read_list.append(readid)
    return read_list


def main():
    samfastqs, out = getargs()
    print(samfastqs)
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 4), layout="constrained")
    for samfastq in samfastqs:
        samfile, fastqfile = samfastq.split(",")
        read_list = choose_unique_to_mRNA_read(samfile)
        len_dic = construct_len_dic(fastqfile)
        len_count_dic = {}
        for readid in read_list:
            length = len_dic[readid]
            if length in len_count_dic:
                len_count_dic[length] += 1
            else:
                len_count_dic[length] = 1
        keys = list(len_count_dic.keys())
        keys.sort()
        len_list = []
        for kk in keys:
            len_list.append(len_count_dic[kk])
        
        axs[0].plot(keys, len_list, label=os.path.basename(samfile).split("-")[0])
        len_list = np.asarray(len_list)
        len_list_pct = len_list / np.sum(len_list)
        axs[1].plot(keys, len_list_pct, label=os.path.basename(samfile).split("-")[0])
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