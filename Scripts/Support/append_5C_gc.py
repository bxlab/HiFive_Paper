#!/usr/bin/env python

import sys


def main():
    bed_fname, fa_fname, out_fname = sys.argv[1:4]
    bed_data = []
    for line in open(bed_fname):
        bed_data.append(line[:-1].split('\t'))
    # if header in bed_data, remove
    try:
        int(bed_data[0][1])
    except:
        bed_data = bed_data[1:]
    name = ''
    fa_data = {}
    for line in open(fa_fname):
        if line[0] == '>':
            name = line[1:-1]
        else:
            seq = line[:-1].upper()
            seq.replace('TAATACGACTCACTATAGCC', '')
            seq.replace('TCCCTTTAGTGAGGGTTAATA', '')
            gc = float(seq.count('G') + seq.count('C')) / len(seq)
            fa_data[name] = gc
    for i in range(len(bed_data)):
        if bed_data[i][3] in fa_data:
            bed_data[i].append(str(fa_data[bed_data[i][3]]))
        else:
            bed_data[i].append('0.0')
    output = open(out_fname, 'w')
    print >> output, "chr\tstart\tstop\tname\tscore\tstrand\tgc"
    for line in bed_data:
        print >> output, '\t'.join(line)
    output.close()

if __name__ == "__main__":
    main()
