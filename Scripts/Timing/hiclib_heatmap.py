#!/usr/bin/env python

import sys

from hiclib import mapping, fragmentHiC
from mirnylib import h5dict, genome
import h5py

basedir = sys.argv[1]
genome_db    = genome.Genome('%s/Data/Genome/mm9_fasta' % basedir, readChrms=['1'], chrmFileTemplate="%s.fa")
temp = h5py.File('%s/Data/Timing/hiclib_data_norm.hdf5' % basedir, 'r')
weights = temp['weights'][...]
temp.close()
fragments = fragmentHiC.HiCdataset(
    filename='temp',
    genome=genome_db,
    maximumMoleculeLength=500,
    mode='a',
    enzymeName="NcoI",
    inMemory=True)
fragments.load('%s/Data/Timing/hiclib_data_norm.hdf5' % basedir)
fragments.weights = weights
fragments.fragmentWeights = weights
fragments.vectors['weights'] = 'float32'

fragments.saveHeatmap('%s/Data/Timing/hiclib_heatmap.hdf5' % basedir, resolution=10000, useWeights=True)