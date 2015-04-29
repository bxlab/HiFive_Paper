#!/usr/bin/env python
#(c) 2014 Emory University. All Rights Reserved
# Code written by: Michael Sauria (mgehrin@emory.edu)

import sys
import glob


def main():
    stats_prefix, out_fname = sys.argv[1:3]
    stats_fnames = glob.glob("%s*" % stats_prefix)
    ordered_fnames = {}
    for fname in stats_fnames:
        ordered_fnames[int(fname.split('.')[-1])] = fname
    stats = {'processed':0, 'alignment':0, '-m':0}
    keys = ordered_fnames.keys()
    keys.sort()
    output = open(out_fname, 'w')
    for key in keys:
        print >> output, ordered_fnames[key].split('/')[-1]
        for line in open(ordered_fnames[key], 'r'):
            if line.startswith('#'):
                temp = line.strip('\n\r')
                print >> output, temp
                name = line.split(':')[0].split(' ')[-1]
                if name == 'processed':
                   stats[name] = max(stats[name], int(line.split(':')[1].split('(')[0].strip(' ')))
                elif name == 'alignment' or name == '-m':
                    stats[name] += int(line.split(':')[1].split('(')[0].strip(' '))
    print >> output, "Overall Stats"
    print >> output, "# reads processed: %i" % (stats['processed'])
    print >> output, "# reads with at least one reported alignment: %i (%0.2f%%)" % (stats['alignment'],
                     100.0 * stats['alignment'] / float(stats['processed']))
    print >> output, "# reads that failed to align: %i (%0.2f%%)" % (stats['processed'] - stats['-m'] -
                     stats['alignment'], 100.0 * (stats['processed'] - stats['-m'] - stats['alignment']) /
                     float(stats['processed']))
    print >> output, "# reads with alignments suppressed due to -m: %i (%0.2f%%)" % (stats['-m'],
                     100.0 * stats['-m'] / float(stats['processed']))
    output.close()
    return None


if __name__ == '__main__':
    main()
