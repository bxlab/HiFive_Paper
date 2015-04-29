#!/usr/bin/env python

import sys

import h5py


def main():
    infname, outfname = sys.argv[1:3]
    hcd = h5py.File(infname, 'r')
    output = open(outfname, 'w')
    print >> output, "fend1\tfend2\tcount"
    data = hcd['cis_data'][...]
    for i in range(data.shape[0]):
        print >> output, "%i\t%i\t%i" % (data[i, 0] + 1, data[i, 1] + 1, data[i, 2])
    data = hcd['trans_data'][...]
    hcd.close()
    for i in range(data.shape[0]):
        print >> output, "%i\t%i\t%i" % (data[i, 0] + 1, data[i, 1] + 1, data[i, 2])
    output.close()


if __name__ == "__main__":
    main()
