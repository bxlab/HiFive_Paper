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
	line_num = 0
	primers = []
	input = open(in_fname, 'r')
	for line in input:
		if line_num % 2 == 0:
			primer = "%s_%04i" % ('_'.join(line.strip('>').split('_')[:-1]), int(line.split('_')[-1].split('|')[0]))
			if not names is None:
				count = 0
				for name in names:
					count += primer.count(name)
			if count > 0:
				chr = line.split('|')[-1].split(':')[0]
				start = str(int(line.split(':')[1].split('-')[0]) - 2)
				stop = str(int(line.strip('\n').split(':')[1].split('-')[1]) - 2)
				if primer.split('_')[3] == 'FOR':
					primers.append([primer, chr, start, stop, '+'])
				else:
					primers.append([primer, chr, start, stop, '-'])
		line_num += 1
	input.close()
	primers.sort()
	output = open(out_fname, 'w')
	for line in primers:
		print >> output, '\t'.join(line[1:4] + [line[0], '0', line[-1]])
	output.close()
	return None


if __name__ == '__main__':
	main()
