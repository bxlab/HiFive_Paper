#!/usr/bin/env python

import sys

from hiclib import fragmentHiC
from mirnylib import genome
import h5py
import numpy


fasta_dir, re_name, data_fname, out_fname, binsize, trans = sys.argv[1:7]
binsize = int(binsize)
if trans in ['1', 'true', 'True', 'TRUE']:
	trans = True
else:
	trans = False
genome_db    = genome.Genome(fasta_dir, readChrms=['#', 'X'], chrmFileTemplate="%s.fa")
temp = h5py.File(data_fname, 'r')
if 'weights' in temp:
    weights = temp['weights'][...]
else:
    weights = temp['fragmentWeights'][...]
temp.close()
fragments = fragmentHiC.HiCdataset(
    filename='temp',
    genome=genome_db,
    maximumMoleculeLength=500,
    mode='a',
    enzymeName=re_name,
    inMemory=True)
fragments.load(data_fname)
fragments.weights = weights
fragments.fragmentWeights = weights
fragments.vectors['weights'] = 'float32'

if trans:
    heatmap = fragments.buildAllHeatmap(resolution=binsize, useWeights=True)
    fragments.genome.setResolution(binsize)
    chr_indices = numpy.r_[fragments.genome.chrmStartsBinCont(numpy.arange()), heatmap.shape[0]]
    output = h5py.File(out_fname, 'w')
    for i in range(chr_indices.shape[0] - 1):
        positions = numpy.zeros((chr_indices[i + 1] - chr_indices[i], 2), dtype=numpy.int32)
        positions[:, 0] = numpy.arange(positions.shape[0]) * binsize
        positions[:, 1] = positions[:, 0] + binsize
        output.create_dataset(name='%s.positions' % (int2chr[i]), data=positions)
        indices = numpy.triu_indices(chr_indices[i + 1] - chr_indices[1], 1)
        output.create_dataset(name='%s.enrichments' % (int2chr[i]),
                data=heatmap[indices[0] + chr_indices[i], indices[1] + chr_indices[i]])
        for j in range(i + 1, chr_indices.shape[0] - 1):
            output.create_dataset(name='%s_by_%s.enrichments' % (int2chr[i], int2chr[j]),
                    data=heatmap[chr_indices[i]:chr_indices[i + 1], chr_indices[j]:chr_indices[j + 1]])
    output.close()
else:
    