import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="out", help="out file name")
    parser.add_argument("--compare", help="column to do correlation plot. ie, 1,3")
    parser.add_argument("countfile", help="count file")
    args = parser.parse_args()
    return args.out, args.countfile, args.compare


def plot_corr(dt, x_idx, x_label, y_idx, y_label, ax):
    x = dt[:, x_idx]
    y = dt[:, y_idx]
    ax.scatter(x, y, s=2)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xscale("log")
    ax.set_yscale("log")


def plot_count_dist(dt, ax, head):
    for idx, head_label in enumerate(head):
        count = {}
        col = dt[:, idx]
        for e in col:
            if e in count:
                count[e] += 1
            else:
                count[e] = 1        
        keys = list(count.keys())
        keys.sort()
        val = []
        for k in keys:
            val.append(count[k])
        ax.plot(keys[1:], val[1:], label=head_label, lw=.7)
    ax.set_xlabel("read count")
    ax.set_xscale("log")
    ax.legend()


def main():
    out, countfile, compare = getargs()
    compare = [int(e) -1 for e in compare.split(",")]
    dt = pd.read_csv(countfile, header=0, index_col="geneid", sep="\t")
    head = dt.columns
    dt = dt.to_numpy()

    fig, axs = plt.subplots(nrows=1, ncols=2, layout="constrained", figsize=(8, 4))

    plot_corr(dt, compare[0], head[compare[0]], compare[1], head[compare[1]], axs[0])
    plot_count_dist(dt, axs[1], head)

    fig.savefig(out + "-count-stats.svg")


if __name__ == "__main__":
    main()