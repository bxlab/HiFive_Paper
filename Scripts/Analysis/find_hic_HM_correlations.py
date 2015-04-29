#!/usr/bin/env python

import sys
import os

import numpy
import glob
try:
    from mpi4py import MPI
except:
    pass

import h5py


def main():
    infname1, infname2, chroms, num_steps, trans, max_size = sys.argv[1:7]
    num_steps = int(num_steps)
    max_size = int(max_size)
    chroms = chroms.split(',')
    if trans.lower() in ['true', '1']:
        trans = True
        needed = []
    else:
        trans = False
        needed = list(chroms)
    if trans:
        for i in range(len(chroms) - 1):
            for j in range(i + 1, len(chroms)):
                needed.append((chroms[i], chroms[j]))
    res_str = infname1.split('_')[-1].split('.')[0]
    res = int(res_str.replace('M','000000').replace('K','000'))
    infile1 = h5py.File(infname1, 'r')
    infile2 = h5py.File(infname2, 'r')
    if rank == 0:
        positions = {}
        for chrom in chroms:
            pos1 = infile1['%s.positions' % chrom][...]
            pos2 = infile2['%s.positions' % chrom][...]
            where = numpy.searchsorted(pos2[:-1, 0], pos1[:, 0])
            positions[chrom] = [pos1[numpy.where(pos1[:, 0] == pos2[where, 0])[0], 0]]
            positions[chrom].append(numpy.searchsorted(pos1[:, 0], positions[chrom][0]))
            positions[chrom].append(numpy.searchsorted(pos2[:, 0], positions[chrom][0]))
            positions[chrom].append(numpy.zeros(pos1.shape[0], dtype=numpy.int32) - 1)
            positions[chrom][-1][positions[chrom][1]] = numpy.arange(positions[chrom][1].shape[0])
            positions[chrom].append(numpy.zeros(pos2.shape[0], dtype=numpy.int32) - 1)
            positions[chrom][-1][positions[chrom][2]] = numpy.arange(positions[chrom][2].shape[0])
        for i in range(1, num_procs):
            comm.send(positions, dest=i, tag=11)
    else:
        positions = comm.recv(source=0, tag=11)
    size_cutoffs = numpy.ceil(numpy.exp(numpy.linspace(numpy.log(5.0 * res), numpy.log(max_size), num_steps))).astype(numpy.int64)
    data = {}
    size_bins = {}
    for chrom in needed:
        if isinstance(chrom, str):
            indices = numpy.triu_indices(positions[chrom][0].shape[0], 1)
            node_ranges = numpy.round(numpy.linspace(0, indices[0].shape[0], num_procs + 1)).astype(numpy.int32)
            distances = (positions[chrom][0][indices[1][node_ranges[rank]:node_ranges[rank + 1]]] -
                         positions[chrom][0][indices[0][node_ranges[rank]:node_ranges[rank + 1]]])
            fend1 = positions[chrom][1][indices[0][node_ranges[rank]]]
            fend2 = positions[chrom][1][indices[1][node_ranges[rank]]]
            start1 = fend1 * (positions[chrom][3].shape[0] - 1) - (fend1 * (fend1 + 1)) / 2 - 1 + fend2
            fend1 = positions[chrom][1][indices[0][node_ranges[rank + 1] - 1]]
            fend2 = positions[chrom][1][indices[1][node_ranges[rank + 1] - 1]]
            stop1 = fend1 * (positions[chrom][3].shape[0] - 1) - (fend1 * (fend1 + 1)) / 2 + fend2
            fend1 = positions[chrom][2][indices[0][node_ranges[rank]]]
            fend2 = positions[chrom][2][indices[1][node_ranges[rank]]]
            start2 = fend1 * (positions[chrom][4].shape[0] - 1) - (fend1 * (fend1 + 1)) / 2 - 1 + fend2
            fend1 = positions[chrom][2][indices[0][node_ranges[rank + 1] - 1]]
            fend2 = positions[chrom][2][indices[1][node_ranges[rank + 1] - 1]]
            stop2 = fend1 * (positions[chrom][4].shape[0] - 1) - (fend1 * (fend1 + 1)) / 2 + fend2
            indices = numpy.triu_indices(positions[chrom][3].shape[0], 1)
            valid = numpy.where((positions[chrom][3][indices[0][start1:stop1]] >= 0) *
                                (positions[chrom][3][indices[1][start1:stop1]] >= 0))
            if "%s.counts" % chrom in infile1:
                temp1_1 = infile1["%s.counts" % chrom][start1:stop1][valid]
                temp1 = infile1["%s.expected" % chrom][start1:stop1][valid].astype(numpy.float64)
                nonzero = numpy.where(temp1 > 0.0)
                temp1[nonzero] = temp1_1[nonzero] / temp1[nonzero]
                del temp1_1
                del nonzero
            else:
                temp1 = infile1["%s.enrichments" % chrom][start1:stop1][valid].astype(numpy.float64)
            indices = numpy.triu_indices(positions[chrom][4].shape[0], 1)
            valid = numpy.where((positions[chrom][4][indices[0][start2:stop2]] >= 0) *
                                (positions[chrom][4][indices[1][start2:stop2]] >= 0))
            if "%s.counts" % chrom in infile2:
                temp2_1 = infile2["%s.counts" % chrom][start2:stop2][valid]
                temp2 = infile2["%s.expected" % chrom][start2:stop2][valid].astype(numpy.float64)
                nonzero = numpy.where(temp2 > 0.0)
                temp2[nonzero] = temp2_1[nonzero] / temp2[nonzero]
                del temp2_1
                del nonzero
            else:
                temp2 = infile2["%s.enrichments" % chrom][start2:stop2][valid].astype(numpy.float64)
            valid = numpy.where((temp1 > 0.0) * (temp2 > 0.0))
            data[chrom] = numpy.empty((valid[0].shape[0], 2), dtype=numpy.float64)
            data[chrom][:, 0] = numpy.log(temp1[valid])
            data[chrom][:, 1] = numpy.log(temp2[valid])
            if not trans:
                size_bins[chrom] = numpy.searchsorted(size_cutoffs, distances[valid])
        else:
            chrom1, chrom2 = chrom
            if "%s_by_%s.counts" % (chrom2, chrom1) in infile1 or "%s_by_%s.enrichments" % (chrom2, chrom1) in infile1:
                chrom2, chrom1 = chrom
            node_ranges = numpy.round(numpy.linspace(0, positions[chrom1][0].shape[0], num_procs + 1)).astype(numpy.int32)
            if "%s_by_%s.counts" % (chrom1, chrom2) in infile1:
                temp1 = infile1["%s_by_%s.counts" % (chrom1, chrom2)][positions[chrom1][1][node_ranges[rank]]:(positions[chrom1][1][node_ranges[rank + 1] - 1] + 1), :].astype(numpy.float64)
                temp1_1 = infile1["%s_by_%s.expected" % (chrom1, chrom2)][positions[chrom1][1][node_ranges[rank]]:(positions[chrom1][1][node_ranges[rank + 1] - 1] + 1), :]
                nonzero = numpy.where(temp1_1 > 0.0)
                temp1[nonzero] /= temp1_1[nonzero]
                del temp1_1
                del nonzero
            else:
                temp1 = infile1["%s_by_%s.enrichments" % (chrom1, chrom2)][positions[chrom1][1][node_ranges[rank]]:(positions[chrom1][1][node_ranges[rank + 1] - 1] + 1), :].astype(numpy.float64)
            temp1 = temp1[positions[chrom1][1][node_ranges[rank]:node_ranges[rank + 1]] - positions[chrom1][1][node_ranges[rank]], :]
            temp1 = temp1[:, positions[chrom2][1]].ravel()
            node_ranges = numpy.round(numpy.linspace(0, positions[chrom1][0].shape[0], num_procs + 1)).astype(numpy.int32)
            if "%s_by_%s.counts" % (chrom1, chrom2) in infile2:
                temp2 = infile2["%s_by_%s.counts" % (chrom1, chrom2)][positions[chrom1][2][node_ranges[rank]]:(positions[chrom1][2][node_ranges[rank + 1] - 1] + 1), :].astype(numpy.float64)
                temp2_1 = infile2["%s_by_%s.expected" % (chrom1, chrom2)][positions[chrom1][2][node_ranges[rank]]:(positions[chrom1][2][node_ranges[rank + 1] - 1] + 1), :]
                nonzero = numpy.where(temp2_1 > 0.0)
                temp2[nonzero] /= temp2_1[nonzero]
                del temp2_1
                del nonzero
            else:
                temp2 = infile2["%s_by_%s.enrichments" % (chrom1, chrom2)][positions[chrom1][2][node_ranges[rank]]:(positions[chrom1][2][node_ranges[rank + 1] - 1] + 1), :].astype(numpy.float64)
            temp2 = temp2[positions[chrom1][2][node_ranges[rank]:node_ranges[rank + 1]] - positions[chrom1][2][node_ranges[rank]], :]
            temp2 = temp2[:, positions[chrom2][2]].ravel()
            valid = numpy.where((temp1 > 0.0) * (temp2 > 0.0))
            data[chrom] = numpy.empty((valid[0].shape[0], 2), dtype=numpy.float64)
            data[chrom][:, 0] = numpy.log(temp1[valid])
            data[chrom][:, 1] = numpy.log(temp2[valid])
    if trans:
        means = numpy.zeros(2, dtype=numpy.float64)
        counts = 0
        for chrom in data:
            counts += data[chrom].shape[0]
            means += numpy.sum(data[chrom], axis=0)
        if rank == 0:
            for i in range(1, num_procs):
                means += comm.recv(source=i, tag=11)
                counts += comm.recv(source=i, tag=11)
            means /= counts
            for i in range(1, num_procs):
                comm.send(means, dest=i, tag=11)
        else:
            comm.send(means, dest=0, tag=11)
            comm.send(counts, dest=0, tag=11)
            means = comm.recv(source=0, tag=11)
        stdevs = numpy.zeros(2, dtype=numpy.float64)
        for chrom in data:
            stdevs += numpy.sum((data[chrom] - means.reshape(1, -1)) ** 2.0, axis=0)
        if rank == 0:
            for i in range(1, num_procs):
                stdevs += comm.recv(source=i, tag=11)
            stdevs = stdevs ** 0.5
            for i in range(1, num_procs):
                comm.send(stdevs, dest=i, tag=11)
        else:
            comm.send(stdevs, dest=0, tag=11)
            stdevs = comm.recv(source=0, tag=11)
        corrs = 0.0
        for chrom in data:
            corrs += numpy.sum(numpy.product(data[chrom] - means.reshape(1, -1), axis=1)) / numpy.product(stdevs)
        if rank == 0:
            for i in range(1, num_procs):
                corrs += comm.recv(source=i, tag=11)
            print "%i\tNA\ttrans\t%i\t%f" % (res, counts, corrs)
        else:
            comm.send(corrs, dest=0, tag=11)      
    else:
        means = numpy.zeros((num_steps + 1, 2), dtype=numpy.float64)
        counts = numpy.zeros(num_steps + 1, dtype=numpy.int64)
        for chrom in data:
            means[:-1, 0] += numpy.bincount(size_bins[chrom], weights=data[chrom][:, 0], minlength=num_steps)
            means[:-1, 1] += numpy.bincount(size_bins[chrom], weights=data[chrom][:, 1], minlength=num_steps)
            counts[:-1] += numpy.bincount(size_bins[chrom], minlength=num_steps)
        if rank == 0:
            for i in range(1, num_procs):
                means += comm.recv(source=i, tag=11)
                counts += comm.recv(source=i, tag=11)
            means[-1, :] = numpy.sum(means[:-1, :], axis=0)
            counts[-1] = numpy.sum(counts[:-1])
            means /= counts.reshape(-1, 1)
            for i in range(1, num_procs):
                comm.send(means, dest=i, tag=11)
        else:
            comm.send(means, dest=0, tag=11)
            comm.send(counts, dest=0, tag=11)
            means = comm.recv(source=0, tag=11)
        stdevs = numpy.zeros((num_steps + 1, 2), dtype=numpy.float64)
        for chrom in data:
            stdevs[:-1, 0] += numpy.bincount(size_bins[chrom], weights=(data[chrom][:, 0] - means[size_bins[chrom], 0]) ** 2.0, minlength=num_steps)
            stdevs[:-1, 1] += numpy.bincount(size_bins[chrom], weights=(data[chrom][:, 1] - means[size_bins[chrom], 1]) ** 2.0, minlength=num_steps)
            stdevs[-1, :] += numpy.sum((data[chrom] - means[-1, :].reshape(1, -1)) ** 2.0, axis=0)
        if rank == 0:
            for i in range(1, num_procs):
                stdevs += comm.recv(source=i, tag=11)
            stdevs = stdevs ** 0.5
            for i in range(1, num_procs):
                comm.send(stdevs, dest=i, tag=11)
        else:
            comm.send(stdevs, dest=0, tag=11)
            stdevs = comm.recv(source=0, tag=11)
        corrs = numpy.zeros(num_steps + 1, dtype=numpy.float64)
        for chrom in data:
            corrs[:-1] += numpy.bincount(size_bins[chrom], weights=numpy.product(data[chrom] - means[size_bins[chrom], :], axis=1), minlength=num_steps) / numpy.product(stdevs[:-1, :], axis=1)
            corrs[-1] += numpy.sum(numpy.product(data[chrom] - means[-1, :].reshape(1, -1), axis=1)) / numpy.product(stdevs[-1, :])
        if rank == 0:
            for i in range(1, num_procs):
                corrs += comm.recv(source=i, tag=11)
            for i in range(num_steps):
                print "%i\t%i\tcis\t%i\t%f" % (res, size_cutoffs[i], counts[i], corrs[i])
            print "%i\tNA\tcis\t%i\t%f" % (res, counts[-1], corrs[-1])
        else:
            comm.send(corrs, dest=0, tag=11)


if "mpi4py" in sys.modules.keys():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    num_procs = comm.Get_size()
else:
    comm = None
    rank = 0
    num_procs = 1

if __name__ == "__main__":
    main()
