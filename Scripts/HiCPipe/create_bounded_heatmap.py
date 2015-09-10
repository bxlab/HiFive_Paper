#!/usr/bin/env python

import sys
import subprocess
import os
from math import floor, ceil

import numpy
import h5py

import hifive
import heatmap


def main():
    datafname, modelprefix, outfname, fivecfname, regions = sys.argv[1:6]
    regions = regions.split(',')
    fivec = hifive.FiveC(fivecfname, 'r')
    data = h5py.File(datafname, 'r')
    # get fend range from data file
    data_path = os.path.dirname(datafname).split('/')
    fendfilename = data['/'].attrs['fendfilename']
    if fendfilename.count('../') > 0:
        data_path = data_path[:-fendfilename.count('../')]
    fendfilename = '/'.join(data_path + fendfilename.split('/')[-1:])
    fends = h5py.File(fendfilename, 'r')
    output = h5py.File(outfname, 'w')
    for region in regions:
        region = int(region)
        # get region bounds
        index, start_frag, stop_frag, chrom, start, stop = fivec.frags['regions'][region]
        bounds = numpy.zeros((stop_frag - start_frag, 2), dtype=numpy.int32)
        bounds[:, 0] = fivec.frags['fragments']['start'][start_frag:stop_frag]
        bounds[:, 1] = fivec.frags['fragments']['stop'][start_frag:stop_frag] + 1
        chr_indices = fends['chr_indices'][...]
        chromosomes = fends['chromosomes'][...]
        chr2int = {}
        mids = fends['fends']['mid'][...]
        for i in range(chromosomes.shape[0]):
            chr2int[chromosomes[i]] = i
        chrint = chr2int[chrom]
        start_fend = chr_indices[chrint]
        stop_fend = chr_indices[chrint + 1]
        while start_fend < stop_fend and mids[start_fend] < bounds[0, 0]:
            start_fend += 1
        while stop_fend > start_fend and mids[stop_fend - 1] > bounds[-1, 1]:
            stop_fend -= 1
        num_fends = stop_fend - start_fend
        # get data and indices
        start_index = data['cis_indices'][start_fend]
        stop_index = data['cis_indices'][stop_fend]
        cis_data = data['cis_data'][start_index:stop_index, :]
        cis_indices = data['cis_indices'][start_fend:(stop_fend + 1)]
        cis_indices -= cis_indices[0]
        print >> sys.stderr, ("Creating binned heatmap..."),
        # load model parameters
        gc_matrix = load_matrix("%s_frag_gc_bin.f" % modelprefix)
        len_matrix = load_matrix("%s_frag_len_bin.f" % modelprefix)
        map_matrix = load_matrix("%s_map_bin.f" % modelprefix)
        fend_bins = load_fend_bins("%s.binned" % modelprefix, mids.shape[0])
        # allocate results arrays
        num_bins = bounds.shape[0]
        observed = numpy.zeros(num_bins * (num_bins - 1) / 2, dtype=numpy.int32)
        expected = numpy.zeros(num_bins * (num_bins - 1) / 2, dtype=numpy.float32)
        unbinned_observed = numpy.zeros(num_fends * (num_fends - 1) / 2, dtype=numpy.int32)
        unbinned_expected = numpy.zeros(num_fends * (num_fends - 1) / 2, dtype=numpy.float32)
        # determine bins for fends
        mids =  fends['fends']['mid'][start_fend:stop_fend]
        start_indices = numpy.searchsorted(bounds[:, 1], mids, side='right')
        stop_indices = numpy.searchsorted(bounds[:, 0], mids, side='right')
        where = numpy.where(start_indices < stop_indices)[0]
        bins = numpy.zeros(mids.shape[0], dtype=numpy.int32) - 1
        bins[where] = start_indices[where]
        # Find cis heatmap
        heatmap.bin_bounds_cis(cis_data, cis_indices, bins, fend_bins, map_matrix, gc_matrix, len_matrix,
                                observed, expected, start_fend, stop_fend, bounds.shape[0])
        heatmap.bin_bounds_cis(cis_data, cis_indices, numpy.arange(num_fends).astype(numpy.int32), fend_bins, map_matrix, gc_matrix, len_matrix,
                                unbinned_observed, unbinned_expected, start_fend, stop_fend, num_fends)
        print >> sys.stderr, ("Done\n"),
        # Write results to file
        print >> sys.stderr, ("Writing to file..."),
        output.create_dataset(name='%i.counts' % region, data=observed)
        output.create_dataset(name='%i.expected' % region, data=expected)
        output.create_dataset(name='%i.unbinned_counts' % region, data=unbinned_observed)
        output.create_dataset(name='%i.unbinned_expected' % region, data=unbinned_expected)
        output.create_dataset(name='%i.bounds' % region, data=bounds)
        output.create_dataset(name='%i.mids' % region, data=mids)
        output.attrs['%i.chromosome' % region] = chrom
    output.close()
    print >> sys.stderr, ("Done\n"),
    return None


def load_matrix(fname):
    input = open(fname, 'r')
    data = []
    for line in input:
        data.append(line.strip('\n').split('\t'))
    input.close()
    data = data[1:]
    for i in range(len(data)):
        for j in range(len(data[i])):
            data[i][j] = float(data[i][j])
    data = numpy.array(data, dtype=numpy.float32)
    num_bins = int(numpy.amax(data[:,0]))
    matrix = numpy.zeros((num_bins, num_bins), dtype=numpy.float32)
    for i in range(data.shape[0]):
        matrix[int(data[i, 0]) - 1, int(data[i, 1]) - 1] = data[i, 2]
    matrix += matrix.T
    matrix[numpy.arange(num_bins), numpy.arange(num_bins)] /= 2.0
    return matrix


def load_fend_bins(fname, num_fends):
    input = open(fname, 'r')
    bins = numpy.zeros((num_fends, 4), dtype=numpy.int32)
    for line in input:
        temp = line.strip('\n').split('\t')
        if temp[0] == 'fend':
            continue
        fend = int(temp[0]) - 1
        map_bin = int(temp[8]) - 1
        gc_bin = int(temp[9]) - 1
        len_bin = int(temp[10]) - 1
        bins[fend, 0] = 1
        bins[fend, 1] = map_bin
        bins[fend, 2] = gc_bin
        bins[fend, 3] = len_bin
    input.close()
    return bins


if __name__ == "__main__":
    main()
