import argparse


def getargs():
    parser = argparse.ArgumentParser(description="statstic read map to positive oritation")
    parser.add_argument("samfile", help="same file")
    args = parser.parse_args()
    return args.samfile


def main():
    samfile = getargs()
    readdic = {}
    
    fin = open(samfile, "r")
    while True:
        line = fin.readline()
        if line[0] != "@":
            break
    
    line = line.rstrip().split()
    if line[1] == "99" and (line[0] not in readdic):
        readdic[line[0]] = 1
    
    for line in fin:
        line = line.rstrip().split()
        if line[1] == "99" and (line[0] not in readdic):
            readdic[line[0]] = 1

    fin.close()
    print(len(readdic))


if __name__ == "__main__":
    main()
