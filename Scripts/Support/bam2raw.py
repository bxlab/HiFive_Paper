#!/usr/bin/env python

import sys

import pysam


def main():
    in_fname1, in_fname2, out_fname = sys.argv[1:4]
    unpaired = {}
    input = pysam.Samfile(in_fname1, 'rb')
    for read in input.fetch(until_eof=True):
        # Only consider reads with an alignment
        if read.is_unmapped:
            continue
        chrom = input.getrname(read.tid).strip('chr')
        # skip multiply-aligned reads
        for tag in read.tags:
            if tag[0] == 'XS':
                break
        else:
            if read.is_reverse:
                strand = "-"
                end = read.pos + len(read.seq)
            else:
                strand = "+"
                end = read.pos
            unpaired[read.qname] = "%s\t%i\t%s" % (chrom, end, strand)
    input.close()
    output = open(out_fname, 'w')
    input = pysam.Samfile(in_fname2, 'rb')
    for read in input.fetch(until_eof=True):
        # Only consider reads with an alignment
        if read.is_unmapped:
            continue
        chrom = input.getrname(read.tid).strip('chr')
        # skip multiply-aligned reads
        for tag in read.tags:
            if tag[0] == 'XS':
                break
        else:
            if read.is_reverse:
                strand = "-"
                end = read.pos + len(read.seq)
            else:
                strand = "+"
                end = read.pos
            if read.qname in unpaired:
                print >> output, "%s\t%s\t%i\t%s" % (unpaired[read.qname], chrom, end, strand)
    input.close()
    output.close()


if __name__ == "__main__":
    main()