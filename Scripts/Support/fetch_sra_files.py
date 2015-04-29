#!/usr/bin/env python
#(c) 2014 Emory University. All Rights Reserved
# Code written by: Michael Sauria (mgehrin@emory.edu)

import sys
import os
import subprocess


file_list, tmp_dir, out_dir, fastq_dump = sys.argv[1:5]
files = []
for line in open(file_list, 'r'):
    line = line.strip()
    if not line or line.startswith('#'):
        continue
    fields = line.split()
    srx = fields[1]
    for srr in fields[2].split(','):
        files.append([srr, srx])
for file in files:
	srr, srx = file
        if (not os.path.exists("%s/%s_1.fastq" % (out_dir, srr)) or
            not os.path.exists("%s/%s_2.fastq" % (out_dir, srr))):
            if not os.path.exists("%s/%s.sra" % (tmp_dir, srr)):
                subprocess.call('wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/%s/%s/%s/%s/%s.sra -O %s' % (srx[:3], srx[:6], srx, srr, srr, "%s/%s.sra" % (tmp_dir, srr)), shell=True)
for file in files:
	srr, srx = file
        if (not os.path.exists("%s/%s_1.fastq" % (out_dir, srr)) or
            not os.path.exists("%s/%s_2.fastq" % (out_dir, srr))):
            subprocess.call('cd %s; %s %s.sra --split-3' % (tmp_dir, fastq_dump, srr), shell=True)
            subprocess.call('mv %s/%s_1.fastq %s/' % (tmp_dir, srr, out_dir), shell=True)
            subprocess.call('mv %s/%s_2.fastq %s/' % (tmp_dir, srr, out_dir), shell=True)
            subprocess.call('rm %s/%s.sra' % (tmp_dir, srr), shell=True)
