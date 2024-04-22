#!
import argparse
import sys
sys.path.append("/home/fanghl/work/git-repo/annoread")
from countmtx import ANNOREAD_DATA

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("annofile", help="read annotation file.")
    parser.add_argument("--out", help="out file", default="out")
    args = parser.parse_args()
    return args.annofile, args.out


def get_start_stop_codon_pos(transcript_range_rel, start_codon_range, stop_codon_range):
    transcript_range_block_len = [ele[1] - ele[0] + 1 for ele in transcript_range_rel]
    start_most_left = start_codon_range[0][0]
    stop_most_left = stop_codon_range[0][0]
    find = 0
    for idx_start, block in enumerate(transcript_range_rel):
        if start_most_left <= block[1] and start_most_left >= block[0]:
            find = 1
            break
    assert find == 1

    find = 0
    for idx_stop, block in enumerate(transcript_range_rel):
        if stop_most_left <= block[1] and stop_most_left >= block[0]:
            find = 1
            break
    assert find == 1

    start_pos = sum(transcript_range_block_len[: idx_start]) + \
        start_most_left - transcript_range_rel[idx_start][0] + 1
    stop_pos = sum(transcript_range_block_len[: idx_stop]) + \
        stop_most_left - transcript_range_rel[idx_stop][0] + 1
    trans_len = sum(transcript_range_block_len)
    return trans_len, start_pos, stop_pos


def calculate_relative_position(annodata, out):
    fout = open(out + ".rdpos", "w")
    data = annodata.get_data()
    for gene_name in data:
        print(">" + gene_name, file=fout)
        for gene_id in data[gene_name]:
            print("&" + gene_id, file=fout)
            for transcript_id in data[gene_name][gene_id]["transcripts"]:
                gbkey = data[gene_name][gene_id]["transcripts"][transcript_id]["gbkey"]
                transcript_range = data[gene_name][gene_id]["transcripts"][transcript_id]["transcript_range"]
                start_codon_range = data[gene_name][gene_id]["transcripts"][transcript_id]["start_codon"]
                stop_codon_range = data[gene_name][gene_id]["transcripts"][transcript_id]["stop_codon"]
                reads = data[gene_name][gene_id]["transcripts"][transcript_id]["reads"]
                if gbkey == "mRNA" and start_codon_range and stop_codon_range:
                    trans_len, start_pos, stop_pos = get_start_stop_codon_pos(transcript_range, 
                                                            start_codon_range, stop_codon_range)
                    print("\t".join(["$" + transcript_id, gbkey, str(trans_len), str(start_pos), str(stop_pos)]), file=fout)
                    for read in reads:
                        read_id, read_pos, read_len, read_align = read
                        rel_start = read_pos - start_pos
                        rel_stop = read_pos - stop_pos
                        print('\t'.join([read_id, str(read_pos), str(read_len), str(rel_start), str(rel_stop)]), file=fout)
    
    fout.close()


if __name__ == "__main__":
    """
    Read annofile, and calculate the distance of bwtween first base
    of start or stop codon and first base of read(5' side).
    All transcripts need be mRNA.
    
    output:
    >gene_name
    &gene_id
    $transcript_id transcript_len, start_codon_pos, stop_codon_pos
    read_id, read_len, pos_rel_to_start, pos_rel_to_stop
    """
    annofile, out = getargs()
    print(annofile)
    print("reading...")
    annodata = ANNOREAD_DATA(annofile)
    print("calculating and output...")
    calculate_relative_position(annodata, out)

