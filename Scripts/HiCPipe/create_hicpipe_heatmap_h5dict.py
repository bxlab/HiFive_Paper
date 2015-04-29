#!/usr/bin/env python

import sys
import subprocess
import os
from math import floor, ceil

import numpy
try:
    from mpi4py import MPI
except:
    pass
import h5py

import hifive
import heatmap


def write_heatmap_dict(datafname, modelprefix, outfname, binsize=1000000, includetrans=False, chroms=[]):
    # check if MPI is available
    if 'mpi4py' in sys.modules.keys():
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        num_procs = comm.Get_size()
    else:
        comm = None
        rank = 0
        num_procs = 1
    # get chr_indices from data file
    data = h5py.File(datafname, 'r')
    data_path = os.path.dirname(datafname).split('/')
    fendfilename = data['/'].attrs['fendfilename']
    if fendfilename.count('../') > 0:
        data_path = data_path[:-fendfilename.count('../')]
    fendfilename = '/'.join(data_path + fendfilename.split('/')[-1:])
    fends = h5py.File(fendfilename, 'r')
    chr_indices = fends['chr_indices'][...]
    chromosomes = fends['chromosomes'][...]
    chr2int = {}
    mids = fends['fends']['mid'][...]
    for i in range(chromosomes.shape[0]):
        chr2int[chromosomes[i]] = i
    # Check if filename already exists, and remove if it does    
    if rank == 0:
        if os.path.exists(outfname):
            print >> sys.stderr, ("%s already exists, overwriting.") % outfname
            #subprocess.call('rm %s' % outfname, shell=True)
        print >> sys.stderr, ("Creating binned heatmap...\n"),
        output = h5py.File(outfname, 'a')
        output.attrs['resolution'] = binsize
        # If chromosomes not specified, fill list
        if len(chroms) == 0:
            chroms = list(fends['chromosomes'][...])
        # Assemble list of requested arrays
        needed = []
        for chrom in chroms:
            needed.append((chrom,))
        if includetrans:
            for i in range(len(chroms)-1):
                for j in range(i + 1, len(chroms)):
                    if chr2int[chroms[i]] < chr2int[chroms[j]]:
                        needed.append((chroms[i],chroms[j]))
                    else:
                        needed.append((chroms[j],chroms[i]))
        if num_procs == 1:
            node_needed = needed
        else:
            worker_size = int(ceil(len(needed) / float(num_procs)))
            for i in range(1, num_procs):
                comm.send(needed[(worker_size * (i - 1)):(worker_size * i)], dest=i, tag=11)
            node_needed = needed[(worker_size * (num_procs - 1)):]
    else:
        node_needed = comm.recv(source=0, tag=11)
    # load model parameters
    gc_matrix = load_matrix("%s_frag_gc_bin.f" % modelprefix)
    len_matrix = load_matrix("%s_frag_len_bin.f" % modelprefix)
    map_matrix = load_matrix("%s_map_bin.f" % modelprefix)
    fend_bins = load_fend_bins("%s.binned" % modelprefix, mids.shape[0])
    # Determine the requested data type
    heatmaps = {}
    # Find heatmaps
    for chrom in node_needed:
        if len(chrom) == 1:
            # Find cis heatmap
            start_fend = chr_indices[chr2int[chrom[0]]]
            stop_fend = chr_indices[chr2int[chrom[0]] + 1] 
            temp_indices = data['cis_indices'][start_fend:(stop_fend + 1)]
            temp_data = data['cis_data'][temp_indices[0]:temp_indices[-1], :]
            temp_indices -= temp_indices[0]
            num_bins = mids[stop_fend - 1] / binsize + 1
            observed = numpy.zeros(num_bins * (num_bins - 1) / 2, dtype=numpy.int32)
            expected = numpy.zeros(num_bins * (num_bins - 1) / 2, dtype=numpy.float32)
            heatmap.bin_cis(temp_data, temp_indices, mids, fend_bins, map_matrix, gc_matrix, len_matrix,
                            observed, expected, start_fend, num_bins, binsize)
            heatmaps[chrom] = [num_bins, observed, expected]
            print >> sys.stderr, ("%s Done\n") % (chrom[0]),
        else:
            # Find trans heatmap
            start_fend1 = chr_indices[chr2int[chrom[0]]]
            stop_fend1 = chr_indices[chr2int[chrom[0]] + 1]
            start_fend2 = chr_indices[chr2int[chrom[1]]]
            stop_fend2 = chr_indices[chr2int[chrom[1]] + 1]
            temp_indices = data['trans_indices'][start_fend1:(stop_fend1 + 1)]
            temp_data = data['trans_data'][temp_indices[0]:temp_indices[-1], :]
            temp_indices -= temp_indices[0]
            num_bins1 = mids[stop_fend1 - 1] / binsize + 1
            num_bins2 = mids[stop_fend2 - 1] / binsize + 1
            observed = numpy.zeros((num_bins1, num_bins2), dtype=numpy.int32)
            expected = numpy.zeros((num_bins1, num_bins2), dtype=numpy.float32)
            heatmap.bin_trans(temp_data, temp_indices, mids, fend_bins, map_matrix, gc_matrix, len_matrix,
                              observed, expected, start_fend1, start_fend2, stop_fend2, binsize)
            heatmaps[chrom] = [observed, expected]
            print >> sys.stderr, ("%s by %s Done\n") % (chrom[0], chrom[1]),
    # Collect heatmaps at node 0 and write to h5dict
    if rank == 0:
        if num_procs > 1:
            for i in range(1, num_procs):
                temp = comm.recv(source=i, tag=11)
                heatmaps.update(temp)
            del temp
        print >> sys.stderr, ("Writing to file..."),
        for chrom in heatmaps.keys():
            if len(chrom) == 1:
                positions = numpy.zeros((heatmaps[chrom][0], 2), dtype=numpy.int32)
                positions[:, 0] = numpy.arange(heatmaps[chrom][0]) * binsize
                positions[:, 1] = positions[:, 0] + binsize
                if '%s.counts' % chrom[0] in output.keys():
                    output['%s.counts' % chrom[0]][:] = heatmaps[chrom][1][:]
                else:
                    output.create_dataset('%s.counts' % chrom[0], data=heatmaps[chrom][1])
                if '%s.expected' % chrom[0] in output.keys():
                    output['%s.expected' % chrom[0]][:] = heatmaps[chrom][2][:]
                else:
                    output.create_dataset('%s.expected' % chrom[0], data=heatmaps[chrom][2])
                if '%s.positions' % chrom[0] in output.keys():
                    output['%s.positions' % chrom[0]][:, :] = positions
                else:
                    output.create_dataset('%s.positions' % chrom[0], data=positions)
            else:
                if '%s_by_%s.counts' % (chrom[0], chrom[1]) in output.keys():
                    output['%s_by_%s.counts' % (chrom[0], chrom[1])][:, :] = heatmaps[chrom][0][:, :]
                else:
                    output.create_dataset('%s_by_%s.counts' % (chrom[0], chrom[1]), data=heatmaps[chrom][0])
                if '%s_by_%s.expected' % (chrom[0], chrom[1]) in output.keys():
                    output['%s_by_%s.expected' % (chrom[0], chrom[1])][:, :] = heatmaps[chrom][1][:, :]
                else:
                    output.create_dataset('%s_by_%s.expected' % (chrom[0], chrom[1]), data=heatmaps[chrom][1])
        if 'chromosomes' in output.keys():
            output['chromosomes'][:] = chroms[:]
        else:
            output.create_dataset('chromosomes', data=numpy.array(chroms))
        output.close()
        print >> sys.stderr, ("Done\n"),
    else:
        comm.send(heatmaps, dest=0, tag=11)
        del heatmaps
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
    datafname, modelprefix, outfname, binsize, includetrans = sys.argv[1:6]
    binsize = int(binsize)
    if includetrans in ['true','True','1']:
        includetrans = True
    else:
        includetrans = False
    if len(sys.argv) > 6:
        chroms = sys.argv[6].split(',')
    else:
        chroms = []
    write_heatmap_dict(datafname, modelprefix, outfname, binsize, includetrans, chroms)
