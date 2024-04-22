#!
"""
"""
import argparse
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/home/fanghl/work/git-repo/annoread")
from annodata import ANNOREAD_DATA


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("ribo_pos", help="read position file of ribo-seq")
    parser.add_argument("rna_pos", help="read position file of RNA-seq")
    parser.add_argument("gene_id", help="the gene_id which to plot")
    parser.add_argument("--transcript_id", default=None, 
                        help="the transcript should be plotted of this gene")
    parser.add_argument("--plot_meta_gene", action="store_true", 
                        help="Where to plot meta gene")
    parser.add_argument("--out", default="out", help="out file name")
    args = parser.parse_args()
    return (args.ribo_pos, args.rna_pos, args.gene_id, args.transcript_id,
            args.plot_meta_gene, args.out)


def get_gene_data(annodata, gene_id_query):
    data = annodata.get_data()
    gene_data = None
    for gene_name in data:
        for gene_id in data[gene_name]:
            if gene_id == gene_id_query:
                gene_data = data[gene_name][gene_id]
                break
    return gene_data


def check_gene_consistance(ribo_gene_data, rna_gene_data):
    check_flg = True
    gene_info = [None, None, None]
    if ribo_gene_data["ref_id"] != rna_gene_data["ref_id"]:
        check_flg = False
    else:
        gene_info[0] = ribo_gene_data["ref_id"]
    if ribo_gene_data["gene_range"] != rna_gene_data["gene_range"]:
        check_flg = False
    else:
        gene_info[1] = ribo_gene_data["gene_range"]
    if ribo_gene_data["ori"] != rna_gene_data["ori"]:
        check_flg = False
    else:
        gene_info[2] = ribo_gene_data["ori"]
    return check_flg, gene_info


def select_transcript(ribo_gene_data, rna_gene_data, transcript_id):
    ribo_trans = ribo_gene_data["transcripts"]
    rna_trans = rna_gene_data["transcripts"]
    if transcript_id and transcript_id in ribo_trans and transcript_id in rna_trans:
        return transcript_id
    else:
        ribo_trans_read_num = {}
        rna_trans_read_num = {}
        for transcript_id in ribo_trans:
            ribo_trans_read_num[transcript_id] = len(ribo_trans[transcript_id]["reads"])
        for transcript_id in rna_trans:
            rna_trans_read_num[transcript_id] = len(rna_trans[transcript_id]["reads"])

        comm_trans_id = set(ribo_trans_read_num.keys()) & set(rna_trans_read_num.keys())
        trans_read_num = []
        for trans_id in comm_trans_id:
            trans_read_num.append([ribo_trans_read_num[trans_id] + 
                                   rna_trans_read_num[trans_id], trans_id])
        trans_read_num.sort()
        if trans_read_num:
            return trans_read_num[-1][1]
        else:
            return None


def check_transcript_consistance(ribo_gene_data, rna_gene_data, transcript_id):
    ribo_trans_info = ribo_gene_data["transcripts"][transcript_id]
    rna_trans_info = rna_gene_data["transcripts"][transcript_id]
    trans_info = [None, None, None, None]
    check_flg = True
    if ribo_trans_info["gbkey"] != rna_trans_info["gbkey"]:
        check_flg = False
    else:
        trans_info[0] = ribo_trans_info["gbkey"]
    if ribo_trans_info["transcript_range"] != rna_trans_info["transcript_range"]:
        check_flg = False
    else:
        trans_info[1] = ribo_trans_info["transcript_range"]
    if ribo_trans_info["start_codon"] != rna_trans_info["start_codon"]:
        check_flg = False
    else:
        trans_info[2] = ribo_trans_info["start_codon"]
    if ribo_trans_info["stop_codon"] != rna_trans_info["stop_codon"]:
        check_flg = False
    else:
        trans_info[3] = ribo_trans_info["stop_codon"]

    return check_flg, trans_info


def get_read_data(gene_data, transcript_id):
    return gene_data["transcripts"][transcript_id]["reads"]


def get_transcript_joint_info(trans_info):
    trans_range = trans_info[1]
    start_codon = trans_info[2]
    stop_codon = trans_info[3]
    block_len = []
    for ele in trans_range:
        block_len.append(ele[1] - ele[0] + 1)
    block_len_cmu = np.cumsum(block_len)

    trans_len = block_len_cmu[-1]
    junction = block_len_cmu[:-1]

    if start_codon:
        start_codon = start_codon[0][0]
        flg = False
        dist = 0
        for idx, ele in enumerate(trans_range):
            if start_codon >= ele[0] and start_codon <= ele[1]:
                flg = True
                dist = start_codon - ele[0] + 1
                break
        if flg:
            if idx > 0:
                start_codon = block_len_cmu[idx - 1] + dist
            else:
                start_codon = dist
        else:
            start_codon = None


    if stop_codon:
        stop_codon = stop_codon[0][0]
        flg = False
        dist = 0
        for idx, ele in enumerate(trans_range):
            if stop_codon >= ele[0] and stop_codon <= ele[1]:
                flg = True
                dist = stop_codon - ele[0] + 1
                break
        if flg:
            if idx > 0:
                stop_codon = block_len_cmu[idx - 1] + dist
            else:
                stop_codon = dist
        else:
            stop_codon = None

    return trans_len, junction, start_codon, stop_codon


def count_coverage(count_dic, reads):
    for read in reads:
        read_pos, read_len = read[1], read[2]
        for pos in range(read_pos, read_pos + read_len):
            count_dic[pos] += 1


def plot(gene_id, gene_info, selected_trans_id, trans_info, ribo_reads, rna_reads, out):
    print(gene_id)
    print(gene_info)
    print(selected_trans_id)
    print(trans_info)

    trans_len, junction, start_codon, stop_codon = get_transcript_joint_info(trans_info)
    ribo_count = {pos: 0 for pos in range(1, trans_len + 1)}
    rna_count = {pos: 0 for pos in range(1, trans_len + 1)}
    junction_dic = {pos: 0 for pos in range(1, trans_len + 1)}
    for ele in junction:
        junction_dic[ele] = 1

    count_coverage(ribo_count, ribo_reads)
    count_coverage(rna_count, rna_reads)
    
    x = list(range(1, trans_len + 1))
    ribo_y = [ribo_count[ele] for ele in x]
    rna_y = [rna_count[ele] for ele in x]
    junction_y = [junction_dic[ele] for ele in x]
    
    fig, axs = plt.subplots(nrows=4, ncols=1, height_ratios=[4, 5, 1, 5],
                            figsize=(.015 * trans_len, 5),
                            sharex=False,
                            layout="constrained", gridspec_kw={'hspace': 0.01})

    gene_info[1] = ",".join([str(ele) for ele in gene_info[1]])
    gene_str = " ".join([gene_id] + gene_info)
    gbkey = trans_info[0]
    trans_range = ";".join([",".join([str(bb) for bb in ele]) for ele in trans_info[1]])
    if trans_info[2]:
        start_str = ";".join([",".join([str(bb) for bb in ele]) for ele in trans_info[2]])
    else:
        start_str = "*"
    if trans_info[3]:
        stop_str = ";".join([",".join([str(bb) for bb in ele]) for ele in trans_info[3]])
    else:
        stop_str = "*"
    trans_info_str = "  ".join([selected_trans_id] + [gbkey, trans_range, start_str, stop_str])
    ribo_reads_num_str = "RPFs number: " + str(len(ribo_reads))
    rna_reads_num_str = "RNA read number: " + str(len(rna_reads))

    axs[0].set_ylim(0, 5)
    axs[0].text(1, 4, gene_str)
    axs[0].text(1, 3, trans_info_str)
    axs[0].text(1, 2, ribo_reads_num_str)
    axs[0].text(1, 1, rna_reads_num_str)
    axs[0].set_xlim([0, trans_len])
    spines = axs[0].spines
    for bod in spines:
        spines[bod].set_visible(False)
    axs[0].set_yticks([], [])
    axs[0].set_xticks([], [])

    axs[1].bar(x, ribo_y)
    axs[1].set_ylabel("RPFs")
    axs[1].set_xlim([0, trans_len])
    axs[1].set_xticks([], [])

    spines = axs[2].spines
    for bod in spines:
        spines[bod].set_visible(False)
    axs[2].set_yticks([], [])
    axs[2].set_xticks([], [])
    axs[2].bar(x, junction_y, color="green", width=3)
    if start_codon:
        start_y = [0] * trans_len
        start_y[start_codon - 1] = 1
        axs[2].bar(x, start_y, color="blue", width=3)
    if stop_codon:
        stop_y = [0] * trans_len
        stop_y[stop_codon - 1] = 1
        axs[2].bar(x, stop_y, color="red", width=3)
    axs[2].set_xlim([0, trans_len])

    axs[3].bar(x, rna_y)
    axs[3].set_ylabel("RNA")
    axs[3].invert_yaxis()
    axs[3].set_xlim([0, trans_len])

    fig.savefig(out + ".svg")


def main():
    ribo_pos_file, rna_pos_file, gene_id, transcript_id, plot_meta_gene, out = getargs()
    annoread_ribo = ANNOREAD_DATA(ribo_pos_file)
    annoread_rna = ANNOREAD_DATA(rna_pos_file)
    ribo_gene_data = get_gene_data(annoread_ribo, gene_id)
    rna_gene_data = get_gene_data(annoread_rna, gene_id)
    check_flg, gene_info = check_gene_consistance(ribo_gene_data, rna_gene_data)
    assert check_flg

    if plot_meta_gene:
        print("Not support yet")
    else:
        selected_trans_id = select_transcript(ribo_gene_data, rna_gene_data, transcript_id)
        if not selected_trans_id:
            print("not found transcipt id")
            return
        check_flg, trans_info = check_transcript_consistance(ribo_gene_data,
                                    rna_gene_data, selected_trans_id)
        assert check_flg
    
        ribo_reads = get_read_data(ribo_gene_data, selected_trans_id)
        rna_reads = get_read_data(rna_gene_data, selected_trans_id)
        plot(gene_id, gene_info, selected_trans_id, trans_info, ribo_reads, rna_reads, out)


if __name__ == "__main__":
    main()