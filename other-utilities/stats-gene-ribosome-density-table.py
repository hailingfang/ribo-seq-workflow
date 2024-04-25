import pandas as pd
import argparse

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("table", help="ribosome density file filtered")
    args = parser.parse_args()
    return args.table


if __name__ == "__main__":
    data_file = getargs()
    data = pd.read_csv(data_file, header=0, index_col=0, sep="\t")
    data = data.to_numpy()

    total = 0
    red = 0
    below_zero = 0
    for row in data:
        ribo_dens_fc = row[-2]
        ribo_change = row[-4]

        if ribo_change > 1.0:
            red += 1
        if ribo_dens_fc < 0:
            below_zero += 1
        total += 1
    print("\t".join(["total_gene", "light_red_point_ribo_change_facted", "ribo_density_change_factored_below_zero"]))
    print("\t".join([str(total), str(red), str(below_zero)]))
