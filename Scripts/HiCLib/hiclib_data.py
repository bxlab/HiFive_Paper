#!/usr/bin/env python

import sys
import os

from hiclib import mapping, fragmentHiC
from mirnylib import h5dict, genome

fasta_dir, re_name, out_fname, in_dir = sys.argv[1:5]
in_prefices = sys.argv[5:]
basedir = os.path.split(os.path.abspath(out_fname))[0]

mapped_reads = []
for prefix in in_prefices:
    mapped_reads.append(h5dict.h5dict('%s/%s.hdf5' % (basedir, prefix)))
genome_db = genome.Genome(fasta_dir, readChrms=['#', 'X'], chrmFileTemplate="%s.fa")

for i, name in enumerate(mapped_reads):
    mapping.parse_sam(
        sam_basename1="%s/%s_1.bam" % (in_dir, in_prefices[i]),
        sam_basename2="%s/%s_2.bam" % (in_dir, in_prefices[i]),
        out_dict=name,
        genome_db=genome_db, 
        enzyme_name=re_name)

for i, name in enumerate(mapped_reads):
    fragments = fragmentHiC.HiCdataset(
        filename='temp',
        genome=genome_db,
        maximumMoleculeLength=500,
        mode='w',
        enzymeName=re_name,
        inMemory=True)
    fragments.parseInputData(dictLike="%s/%s.hdf5" % (basedir, prefix))
    if i != len(mapped_reads) - 1:
        fragments.save("%s/%s_data.hdf5" % (basedir, prefix))
    else:
        frag_files = []
        for prefix in in_prefices[:-1]:
            frag_files.append("%s/%s_data.hdf5" % (basedir, prefix))
        fragments.merge(frag_files)
        fragments.filterDuplicates()
        fragments.filterExtreme(cutH=0, cutL=0.005)        
        fragments.save(out_fname)
