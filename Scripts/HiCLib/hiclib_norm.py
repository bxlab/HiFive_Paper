#!/usr/bin/env python

import sys

from hiclib import fragmentHiC
from mirnylib import genome

fasta_dir, re_name, data_fname, out_fname = sys.argv[1:5]
genome_db = genome.Genome(fasta_dir, readChrms=['#','X'], chrmFileTemplate="%s.fa")
fragments = fragmentHiC.HiCdataset(
    filename='temp',
    genome=genome_db,
    maximumMoleculeLength=500,
    mode='w',
    enzymeName=re_name,
    inMemory=True)
fragments.load(data_fname)
weights = fragments.iterativeCorrectionFromMax(minimumCount=10)
fragments.save(out_fname)
