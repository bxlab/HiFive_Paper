#!/usr/bin/env python
#(c) 2014 Emory University. All Rights Reserved
# Code written by: Michael Sauria (mgehrin@emory.edu)

import subprocess
import os
import sys
import glob


def main():
    filelist, fastq_dir, tmp_dir, out_dir, bowtie, index_dir, samtools, compile_script, threads = sys.argv[1:12]
    threads = int(threads)
    files = []
    for line in open(filelist, 'r'):
        if line[0] == '#':
            continue
        temp = line[:-1].split(' ')
        for srr in temp[2].split(','):
            if os.path.exists("%s/%s_1.fastq" % (fastq_dir, srr)):
                files.append(["%s/%s_1.fastq" % (fastq_dir, srr), temp[4]])
            if os.path.exists("%s/%s_2.fastq" % (fastq_dir, srr)):
                files.append(["%s/%s_2.fastq" % (fastq_dir, srr), temp[4]])
            if os.path.exists("%s/%s.fastq" % (fastq_dir, srr)):
                files.append(["%s/%s.fastq" % (fastq_dir, srr), temp[4]])
    for fastq_fname, index in files:
        bowtie_index = "%s/%s/%s" % (index_dir, index, index)
        trim_5 = 0
        fastq_prefix = os.path.basename(fastq_fname).split('.')[0]
        process = subprocess.Popen("head -n 2 %s | tail -n 1" % fastq_fname, shell=True, stdout=subprocess.PIPE)
        read_len = len(process.communicate()[0].strip('\n\r'))
        trim_3 = max(0, read_len - trim_5 - 50)
        if read_len >= 40:
            step = 5
        else:
            step = 4
        if not os.path.exists("%s/%s.bam" % (out_dir, fastq_prefix)):
            subprocess.call("rm -f %s/%s.*" % (tmp_dir, fastq_prefix), shell=True)
            subprocess.call("ln -s %s %s/%s.fastq" % (fastq_fname, tmp_dir, fastq_prefix), shell=True)
            first_round = True
            while True:
                if first_round:
                    header = ""
                    first_round = False
                else:
                    header = "--sam-nohead --sam-nosq"
                print >> sys.stderr, ("Mapping:%s\tlength:%i\ttrimmed:%i\n") % (fastq_prefix, read_len, trim_3),
                subprocess.call("%s -3 %i -5 %i --threads %i -m 1 -v 2 --tryhard --sam %s --phred33-quals --un %s/%s.unaligned.%i --max %s/%s.multiple %s %s/%s.fastq >> %s/%s.sam 2> %s/%s.stats.%i" % \
                    (bowtie, trim_3, trim_5, threads, header, tmp_dir, fastq_prefix, trim_3, tmp_dir, fastq_prefix,
                     bowtie_index, tmp_dir, fastq_prefix, tmp_dir, fastq_prefix, tmp_dir, fastq_prefix, trim_3),
                    shell=True)
                subprocess.call("rm %s/%s.fastq" % (tmp_dir, fastq_prefix), shell=True)
                subprocess.call("rm %s/%s.multiple" % (tmp_dir, fastq_prefix), shell=True)
                subprocess.call("ln -s %s/%s.unaligned.%i %s/%s.fastq" %
                    (tmp_dir, fastq_prefix, trim_3, tmp_dir, fastq_prefix), shell=True)
                trim_3 += step
                if read_len - trim_5 - trim_3 < 21:
                    break
            subprocess.call("%s view -bS -o %s/%s.bam %s/%s.sam" % (samtools, out_dir,
                fastq_prefix, tmp_dir, fastq_prefix), shell=True)
            subprocess.call("python %s %s/%s.stats %s/%s.stats" %
                (compile_script, tmp_dir, fastq_prefix, out_dir, fastq_prefix), shell=True)
            subprocess.call("rm %s/%s.*" % (tmp_dir, fastq_prefix), shell=True)
    return None


if __name__ == "__main__":
    main()