import argparse
import matplotlib.pyplot as plt
import numpy as np

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="out", help="out file")

    parser.add_argument("--ribo_data_tr", nargs="+", help="ribo-seq position data")
    parser.add_argument("--ribo_data_ut", nargs="+", help="ribo-seq position data")
    parser.add_argument("--rna_data_tr", nargs="+", help="rna-seq positon data")
    parser.add_argument("--rna_data_ut", nargs="+", help="rna-seq positon data")
    args = parser.parse_args()
    return args.out, args.ribo_data_tr, args.ribo_data_ut, args.rna_data_tr, args.rna_data_ut


def read_data(filename):
    dt = {}
    fin = open(filename, "r")
    key = None
    for line in fin:
        if line[0] == ">" or line[0] == "&":
            continue
        if line[0] == "$":
            line = line[1:].rstrip().split("\t")
            transcript_id, start_codon_pos, stop_codon_pos = line[0], int(line[3]), int(line[4])
            key = (transcript_id, start_codon_pos, stop_codon_pos)
            dt[key] = []
            continue
        line = line.rstrip().split('\t')
        read_len, start_dist, stop_dist = int(line[2]), int(line[3]), int(line[4])
        dt[key].append([read_len, start_dist, stop_dist])

    fin.close()
    return dt


def fetch_data(data, start_dist_range, stop_dist_range, cds_len_cutoff=0, read_len_select=0):
    dt_out_start = []
    dt_out_stop = []
    for key in data:
        transcript_id, start_codon_pos, stop_codon_pos = key
        cds_len = stop_codon_pos - start_codon_pos
        if cds_len > cds_len_cutoff:
            for read in  data[key]:
                read_len, start_dist, stop_dist = read
                if read_len_select != 0 and read_len == read_len_select:
                    if start_dist_range:
                        if start_dist >= start_dist_range[0] and start_dist <= start_dist_range[1]:
                            dt_out_start.append(start_dist)
                    else:
                        dt_out_start.append(start_dist)
                    if stop_dist_range:
                        if stop_dist >= stop_dist_range[0] and stop_dist <= stop_dist_range[1]:
                            dt_out_stop.append(stop_dist)
                    else:
                        dt_out_stop.append(stop_dist)
                else:
                    if start_dist_range:
                        if start_dist >= start_dist_range[0] and start_dist <= start_dist_range[1]:
                            dt_out_start.append(start_dist)
                    else:
                        dt_out_start.append(start_dist)
                    if stop_dist_range:
                        if stop_dist >= stop_dist_range[0] and stop_dist <= stop_dist_range[1]:
                            dt_out_stop.append(stop_dist)
                    else:
                        dt_out_stop.append(stop_dist)
    return dt_out_start, dt_out_stop


def count_data(dt):
    count_dic = {}
    for ele in dt:
        if ele in count_dic:
            count_dic[ele] += 1
        else:
            count_dic[ele] = 1
    key = list(count_dic.keys())
    key.sort()
    val = [count_dic[kk] for kk in key]
    return key, val


def transfer_to_codon(base_pos, counting):
    dt = {}
    for idx, xx in enumerate(base_pos):
        key = xx // 3
        if key in dt:
            dt[key] += counting[idx]
        else:
            dt[key] = counting[idx]

    x = list(dt.keys())
    x.sort()
    y = [dt[xx] for xx in x]

    return x, y


def read_frame(pos, count):
    dt = {}
    for pp, cc in zip(pos, count):
        key = pp % 3
        if key in dt:
            dt[key] +=cc
        else:
            dt[key] = cc
    x = list(dt.keys())
    x.sort()
    y = [dt[xx] for xx in x]
    return x, y


if __name__ == "__main__":
    out, ribo_data_tr, ribo_data_ut, rna_data_tr, rna_data_ut = getargs()

    ribo_data_tr_dt = []
    for ff in ribo_data_tr:
        ribo_data_tr_dt.append(read_data(ff))
    
    ribo_data_ut_dt = []
    for ff in ribo_data_ut:
        ribo_data_ut_dt.append(read_data(ff))
    
    rna_data_tr_dt = []
    for ff in rna_data_tr:
        rna_data_tr_dt.append(read_data(ff))

    rna_data_ut_dt = []
    for ff in rna_data_ut:
        rna_data_ut_dt.append(read_data(ff))

    # the periodicity of number of 5' of read.
    #fig, 3 periodicity
    fig1, axs = plt.subplots(nrows=2, ncols=2, layout="constrained", figsize=(12, 4), sharey=True)

    color = ["blue", "red"]
    alpha = [1, .5]
    start_range = [-24, 24]
    stop_range = [-36, 12]
    for idx, dt in enumerate(ribo_data_ut_dt):
        start_dist, stop_dist = fetch_data(dt, start_range, stop_range, cds_len_cutoff=0, read_len_select=0)
        start_x, start_y = count_data(start_dist)
        stop_x, stop_y = count_data(stop_dist)
        axs[0, 0].bar(start_x, start_y, color=color[idx], alpha=alpha[idx])
        axs[0, 1].bar(stop_x, stop_y, color=color[idx], alpha=alpha[idx])

    for idx, dt in enumerate(ribo_data_tr_dt):
        start_dist, stop_dist = fetch_data(dt, start_range, stop_range, cds_len_cutoff=0, read_len_select=0)
        start_x, start_y = count_data(start_dist)
        stop_x, stop_y = count_data(stop_dist)
        axs[1, 0].bar(start_x, start_y, color=color[idx], alpha=alpha[idx])
        axs[1, 1].bar(stop_x, stop_y, color=color[idx], alpha=alpha[idx])

    for idx, dt in enumerate(rna_data_ut_dt[:2]):
        start_dist, stop_dist = fetch_data(dt, start_range, stop_range, cds_len_cutoff=0, read_len_select=0)
        start_x, start_y = count_data(start_dist)
        stop_x, stop_y = count_data(stop_dist)
        axs[0, 0].plot(start_x, start_y, color=color[idx], alpha=alpha[idx])
        axs[0, 1].plot(stop_x, stop_y, color=color[idx], alpha=alpha[idx])

    for idx, dt in enumerate(rna_data_tr_dt[:2]):
        start_dist, stop_dist = fetch_data(dt, start_range, stop_range, cds_len_cutoff=0, read_len_select=0)
        start_x, start_y = count_data(start_dist)
        stop_x, stop_y = count_data(stop_dist)
        axs[1, 0].plot(start_x, start_y, color=color[idx], alpha=alpha[idx])
        axs[1, 1].plot(stop_x, stop_y, color=color[idx], alpha=alpha[idx])

    axs[0, 0].set_xticks(list(range(-24, 25))[::3], [])
    axs[0, 0].grid(axis="x")
    axs[0, 0].set_ylabel("count")
    axs[0, 1].set_xticks(list(range(-36, 13))[::3], [])
    axs[0, 1].grid(axis="x")
    
    axs[1, 0].set_xticks(list(range(-24, 25))[::3], list(range(-24, 25))[::3])
    axs[1, 0].grid(axis="x")
    axs[1, 0].set_xlabel("distance from start codon /base")
    axs[1, 0].set_ylabel("count")
    axs[1, 1].set_xticks(list(range(-36, 13))[::3], list(range(-36, 13))[::3])
    axs[1, 1].grid(axis="x")
    axs[1, 1].set_xlabel("distance from stop codon /base")

    fig1.suptitle("RPFs/Read positon near start and stop codon\n dsR: first row, un: second row")
    fig1.savefig(out + "-periodicity.svg")


    # frame
    fig2, ax = plt.subplots(layout="constrained", figsize=(4, 4))
    start_range = [-12, 72]
    ribo_read_frame = []
    rna_read_frame = []
    for idx, dt in enumerate(ribo_data_ut_dt):
        start_dist, stop_dist = fetch_data(dt, start_range, None, cds_len_cutoff=0, read_len_select=31)
        start_x, start_y = count_data(start_dist)
        start_x, start_y = read_frame(start_x, start_y)
        ribo_read_frame.append(start_y)
    
    for idx, dt in enumerate(ribo_data_tr_dt):
        start_dist, stop_dist = fetch_data(dt, start_range, None, cds_len_cutoff=0, read_len_select=31)
        start_x, start_y = count_data(start_dist)
        start_x, start_y = read_frame(start_x, start_y)
        ribo_read_frame.append(start_y)

    for idx, dt in enumerate(rna_data_ut_dt):
        start_dist, stop_dist = fetch_data(dt, None, None, cds_len_cutoff=0, read_len_select=0)
        start_x, start_y = count_data(start_dist)
        start_x, start_y = read_frame(start_x, start_y)
        rna_read_frame.append(start_y)

    for idx, dt in enumerate(rna_data_tr_dt):
        start_dist, stop_dist = fetch_data(dt, None, None, cds_len_cutoff=0, read_len_select=0)
        start_x, start_y = count_data(start_dist)
        start_x, start_y = read_frame(start_x, start_y)
        rna_read_frame.append(start_y)
    
    ribo_read_frame_mean = np.mean(np.asarray(ribo_read_frame), axis=0)
    rna_read_frame_mean = np.mean(np.asarray(rna_read_frame), axis=0)
    ribo_read_frame_norm = ribo_read_frame_mean / np.sum(ribo_read_frame_mean)
    rna_read_frame_norm = rna_read_frame_mean / np.sum(rna_read_frame_mean)
    x = np.array([0, 1, 2])
    ax.bar(x - .16, ribo_read_frame_norm, width=.3, color="blue")
    ax.bar(x + .16, rna_read_frame_norm, width=.3, color="red")
    print("distance of 5' of reads from start codon. The distance [-12, 72] is analyzed.")
    print("fraction of ribo-seq data", ribo_read_frame_norm)
    print("fraction of rna-seq data", rna_read_frame_norm)
    ax.set_title("RPFs/Reads frames")
    ax.set_ylabel("Fraction of 31mer reads")
    ax.set_xlabel("Frame")
    ax.set_xticks([0, 1, 2], [0, 1, 2])
    fig2.savefig(out + "-read-frame.svg")

    #read distribution normalised.
    fig3 = plt.figure(figsize=(12, 4), layout="constrained")
    fig3_up, fig3_bot = fig3.subfigures(nrows=2, ncols=1, height_ratios=[4, 2])
    axs_up = fig3_up.subplots(nrows=1, ncols=2, sharey=True)
    axs_bot = fig3_bot.subplots(nrows=1, ncols=2, sharey=True)

    start_range = [-100, 300]
    stop_range = [-300, 100]
    alpha = .25
    for idx, dt in enumerate(ribo_data_ut_dt):
        start_dist, stop_dist = fetch_data(dt, start_range, stop_range, cds_len_cutoff=600, read_len_select=0)
        start_x, start_y = count_data(start_dist)
        stop_x, stop_y = count_data(stop_dist)
        #start_x, start_y = transfer_to_codon(start_x, start_y)
        #stop_x, stop_y = transfer_to_codon(stop_x, stop_y)
        denominator = np.sum(start_y) + np.sum(stop_y)
        start_y = np.asarray(start_y) / denominator
        stop_y = np.asarray(stop_y) / denominator
        axs_up[0].plot(start_x, start_y, color="blue", alpha=alpha, lw=1, label="Ribo-seq uninj")
        axs_up[1].plot(stop_x, stop_y, color="blue", alpha=alpha, lw=1, label="Ribo-seq uninj")

    for idx, dt in enumerate(ribo_data_tr_dt):
        start_dist, stop_dist = fetch_data(dt, start_range, stop_range, cds_len_cutoff=600, read_len_select=0)
        start_x, start_y = count_data(start_dist)
        stop_x, stop_y = count_data(stop_dist)
        #start_x, start_y = transfer_to_codon(start_x, start_y)
        #stop_x, stop_y = transfer_to_codon(stop_x, stop_y)
        denominator = np.sum(start_y) + np.sum(stop_y)
        start_y = np.asarray(start_y) / denominator
        stop_y = np.asarray(stop_y) / denominator
        axs_up[0].plot(start_x, start_y, color="red", alpha=alpha, lw=1, label="Ribo-seq dsRNA")
        axs_up[1].plot(stop_x, stop_y, color="red", alpha=alpha, lw=1, label="Ribo-seq dsRNA")

    for idx, dt in enumerate(rna_data_ut_dt[:2]):
        start_dist, stop_dist = fetch_data(dt, start_range, stop_range, cds_len_cutoff=600, read_len_select=0)
        start_x, start_y = count_data(start_dist)
        stop_x, stop_y = count_data(stop_dist)
        #start_x, start_y = transfer_to_codon(start_x, start_y)
        #stop_x, stop_y = transfer_to_codon(stop_x, stop_y)
        denominator = np.sum(start_y) + np.sum(stop_y)
        start_y = np.asarray(start_y) / denominator
        stop_y = np.asarray(stop_y) / denominator
        axs_bot[0].plot(start_x, start_y, color="blue", alpha=alpha, lw=1, label="RNA-seq uninj")
        axs_bot[1].plot(stop_x, stop_y, color="blue", alpha=alpha, lw=1, label="RNA-seq uninj")

    for idx, dt in enumerate(rna_data_tr_dt[:2]):
        start_dist, stop_dist = fetch_data(dt, start_range, stop_range, cds_len_cutoff=600, read_len_select=0)
        start_x, start_y = count_data(start_dist)
        stop_x, stop_y = count_data(stop_dist)
        #start_x, start_y = transfer_to_codon(start_x, start_y)
        #stop_x, stop_y = transfer_to_codon(stop_x, stop_y)
        denominator = np.sum(start_y) + np.sum(stop_y)
        start_y = np.asarray(start_y) / denominator
        stop_y = np.asarray(stop_y) / denominator
        axs_bot[0].plot(start_x, start_y, color="red", alpha=alpha, lw=1, label="RNA-seq dsRNA")
        axs_bot[1].plot(stop_x, stop_y, color="red", alpha=alpha, lw=1, label="RNA-seq dsRNA")
      
    
    fig3.suptitle("Read density around start/stop codon\n (the lengths of the CDSs are greater than 600 bps)")
    axs_up[0].set_ylabel("normed count")
    axs_bot[0].set_ylabel("normed count")
    #axs_up[0].set_xticks([0, 100, 200, 300, 400, 500], [])
    #axs_up[1].set_xticks([-500, -400, -300, -200, -100, 0], [])
    axs_bot[0].set_xlabel("distance from start codon /nt")
    axs_bot[1].set_xlabel("distance from stop codon /nt")
    handles, labels = axs_up[0].get_legend_handles_labels()
    axs_up[0].legend([handles[0], handles[2]], [labels[0], labels[2]])
    
    handles, labels = axs_up[1].get_legend_handles_labels()
    axs_up[1].legend([handles[0], handles[2]], [labels[0], labels[2]])
    
    handles, labels = axs_bot[0].get_legend_handles_labels()
    axs_bot[0].legend([handles[0], handles[3]], [labels[0], labels[3]])
    
    handles, labels = axs_bot[1].get_legend_handles_labels()
    axs_bot[1].legend([handles[0], handles[3]], [labels[0], labels[3]])
    
    fig3.savefig(out + "-distribution.svg")


