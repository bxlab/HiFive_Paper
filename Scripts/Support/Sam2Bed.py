#!/usr/bin/env python

import sys


def main():
    in_fname, out_fname = sys.argv[1:3]
    data = {}
    for line in open(in_fname, 'r'):
        if line[0] == "@":
            continue
        temp = line[:-1].split('\t')
        if temp[2] not in data:
            data[temp[2]] = []
            length = int(temp[5][:-1])
        data[temp[2]].append(int(temp[3]))
    chroms = data.keys()
    chroms.sort()
    output = open(out_fname, 'w')
    for chrom in chroms:
        if chrom.count("_") > 0:
            continue
        data[chrom].sort()
        for i in range(len(data[chrom])):
            print >> output, "%s\t%i\t%i\t.\t.\t." % (chrom, data[chrom][i] - 1, data[chrom][i] + length - 1)
    output.close()


if __name__ == "__main__":
    main()
