#!
"""
plot read count correlation between repeats.
"""
import argparse
import matplotlib.pyplot as plt
import numpy as np


def getargs():
    parse = argparse.ArgumentParser()
    parse.add_argument("countmtx", help="the count matrix data")
    parse.add_argument("--col_idx_r1", help="column index of replication 1", type=int)
    parse.add_argument("--col_idx_r2", help="colume index of replication 2", type=int)
    parse.add_argument("--out", default="out", help="out file")
    args = parse.parse_args()
    return args.countmtx, args.col_idx_r1, args.col_idx_r2, args.out


def read_data(filename):
    dt = []
    fin = open(filename, "r")
    fin.readline()
    for line in fin:
        line = line.rstrip().split('\t')
        dt.append([int(ele) for ele in line[1:]])
    fin.close()
    dt = np.asarray(dt)
    return dt


def filter_log2_great_zero(r1, r2):
    dt = np.vstack((r1, r2))
    dt = dt.T
    dt_out = []
    for row in dt:
        if row[0] >= 1 and row[1] >= 1:
            dt_out.append(row)
    dt_out = np.asarray(dt_out)
    return dt_out[:, 0], dt_out[:, 1]


def norm_dt(dt_array):
    dt_sum = sum(dt_array)
    dt_array = dt_array * 1e6
    dt_array /= dt_sum
    dt_array = np.around(dt_array, 5)
    return dt_array


def hist_bar(dt_array):
    counter = {}
    dt_array = np.log2(dt_array)
    max_v = np.max(dt_array)
    min_v = np.min(dt_array)
    print(max_v, min_v)
    linespace = np.linspace(min_v, max_v + .00000001, 50)
    print(linespace)
    bw = linespace[1] - linespace[0]
    for i in range(len(linespace) - 1):
        counter[(linespace[i], linespace[i + 1])] = []
    
    for point in dt_array:
        flg = 1
        for key in counter:
            if point >= key[0] and point < key[1]:
                counter[key].append(point)
                flg = 0
                break
        if flg:
            print(point)
    key_s = list(counter.keys())
    key_s.sort()
    x = []
    y = []
    for key in key_s:
        x.append(key[0])
        y.append(sum(counter[key]))
    return x, y, bw


if __name__ == "__main__":
    countmtx_file, col_idx_r1, col_idx_r2, out = getargs()
    countmtx = read_data(countmtx_file)
    countmtx = countmtx.T
    r1 = countmtx[col_idx_r1]
    r2 = countmtx[col_idx_r2]
    r1 = norm_dt(r1)
    r2 = norm_dt(r2)
    r1, r2 = filter_log2_great_zero(r1, r2)
    corr = np.corrcoef(r1, r2)
    corr = np.round(corr[0, 1], 8)

    fig, axs = plt.subplots(nrows=2, ncols=2, width_ratios=[4, .4], height_ratios=[.4, 4], gridspec_kw={'hspace': 0.01, 'wspace': 0.01}, figsize=(5, 5))
    axs[1, 0].scatter(r1, r2, s=1)
    axs[1, 0].set_xscale("log")
    axs[1, 0].set_yscale("log")
    x1, y1, bw = hist_bar(r1)
    print(x1)
    x2, y2, bw = hist_bar(r2)
    axs[0, 0].bar(x1, y1, color="gray", edgecolor="white", width=bw)
    axs[1, 1].barh(x2, y2, color="gray", edgecolor="white", height=bw)

    spines = axs[0, 1].spines
    axs[0, 1].set_xticks([], [])
    axs[0, 1].set_yticks([], [])
    spines["top"].set_visible(False)
    spines["bottom"].set_visible(False)
    spines["left"].set_visible(False)
    spines["right"].set_visible(False)
    
    spines = axs[0, 0].spines
    axs[0, 0].set_xticks([], [])
    axs[0, 0].set_yticks([], [])
    spines["top"].set_visible(False)
    spines["bottom"].set_visible(False)
    spines["left"].set_visible(False)
    spines["right"].set_visible(False)
    
    spines = axs[1, 1].spines
    axs[1, 1].set_xticks([], [])
    axs[1, 1].set_yticks([], [])
    spines["top"].set_visible(False)
    spines["bottom"].set_visible(False)
    spines["left"].set_visible(False)
    spines["right"].set_visible(False)
    fig.suptitle("Correlation between repeats " + "r=" + str(corr))
    fig.savefig(out + "-correlation.svg")

