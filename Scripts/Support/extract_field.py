#!/usr/bin/python

import sys


def main():
    fname = sys.argv[1]
    if len(sys.argv) > 2:
        names = sys.argv[2:]
    else:
        names = None
    values = {}
    input = open(fname, 'r')
    for line in input:
        if line[0] == '#':
            continue
        fields = line[:-1].split(' ')
        if fields[0] not in values:
            values[fields[0]] = []
        values[fields[0]] += fields[2].split(',')
    input.close()
    if names is None:
        print ' '.join(values.keys())
    else:
        all_values = []
        for name in names:
            for key in values:
                if key.count(name) > 0:
                    all_values += values[key]
        print ' '.join(all_values)


if __name__ == "__main__":
    main()