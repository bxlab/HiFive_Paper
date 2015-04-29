#!/usr/bin/env python

import sys


def main():
	in_fname, out_fname = sys.argv[1:3]
	if len(sys.argv) > 3:
		names = sys.argv[3:]
		count = 0
	else:
		names = None
		count = 1
	input = open(in_fname, 'r')
	output = open(out_fname, 'w')
	for line in input:
		if line[0] == '>':
			if not names is None:
				count = 0
				for name in names:
					count += line.count(name)
			if count > 0:
				print >> output, "%s_%04i" % ('_'.join(line.split('_')[:-1]), int(line.split('_')[-1].split('|')[0]))
		else:
			if count > 0:
				print >> output, line.strip('\n')
	input.close()
	output.close()
	return None


if __name__ == '__main__':
	main()
