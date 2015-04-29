#!/usr/bin/env python

import sys

from hiclib import mapping, fragmentHiC
from mirnylib import h5dict, genome

basedir = sys.argv[1]
genome_db    = genome.Genome('%s/Data/Genome/mm9_fasta' % basedir, readChrms=['1'], chrmFileTemplate="%s.fa")
fragments = fragmentHiC.HiCdataset(
    filename='temp1',
    genome=genome_db,
    maximumMoleculeLength=500,
    mode='w',
    enzymeName="NcoI",
    inMemory=True)
fragments.load('%s/Data/Timing/hiclib_data.hdf5' % basedir)
fragments.filterDuplicates()
fragments.filterExtreme(cutH=0, cutL=0.005)
fragments.save('%s/Data/Timing/hiclib_data_filt.hdf5' % basedir)
