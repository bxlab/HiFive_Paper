#!/usr/bin/env python

import sys
import optparse
from math import ceil, floor
from random import shuffle

import h5py
import numpy
try:
    from mpi4py import MPI
except:
    pass

import hifive
from Matrix_Balancing import BNEWT, BNEWT_sparse, BNEWT_sparse_binary


def main():
    usage = "usage: [mpirun -np N_PROCS] %prog [options] <project_file> <out_file>\n\nArguments:"
    usage += "\n<project_file>  HiFive HiC project file"
    usage += "\n<out_file>      destination for new project file"
    help = {
        "-t":"include trans interactions in calculating weights [default: %default]",
        "-c":"comma-separated list of chromosomes to include in learning [default: all chromosomes]",
        "-q":"silence output messages [default: %default]",
        "-b":"Use binary indicator instead of counts [default %default]",
        }
    if rank == 0:
        parser = optparse.OptionParser(usage=usage)
    else:
        parser = optparse.OptionParser(usage=optparse.SUPPRESS_USAGE, add_help_option=False)
        for key in help:
            help[key] = optparse.SUPPRESS_HELP
        parser.add_option("-h", "--help", help=optparse.SUPPRESS_HELP, dest="help", action="store_true")
    parser.add_option("-t", "--trans", dest="trans", default=False,
                      help=help["-t"], action="store_true")
    parser.add_option("-b", "--binary", dest="binary", default=False,
                      help=help["-b"], action="store_true")
    parser.add_option("-c", "--chromosomes", dest="chroms", default="", metavar="CHROMS", type="string",
                      help=help["-c"], action="store")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help=help["-q"])
    options, args = parser.parse_args()
    if rank != 0:
        options.silent = True
    if len(args) < 2:
        if rank == 0:
            parser.error('incorrect number of arguments')
        else:
            sys.exit(1)
    hic = hifive.HiC(args[0], 'r', silent=options.silent)
    if len(options.chroms) == 0:
        chroms = list(hic.fends['chromosomes'][...])
    else:
        chroms = options.chroms.split(',')
    chr_indices = hic.fends['chr_indices'][...]
    if options.trans:
        if rank > 0:
            return None
        filt = numpy.copy(hic.filter)
        for i, chrom in enumerate(hic.fends['chromosomes'][...]):
            if chrom not in options.chroms:
                filt[chr_indices[i]:chr_indices[i + 1]] = 0
        data = hic.data['cis_data'][...]
        if options.trans:
            data = numpy.vstack((data, hic.data['trans_data'][...]))
        valid = numpy.where(filt[data[:, 0]] * filt[data[:, 1]])[0]
        data = data[valid, :]
        interactions = numpy.bincount(data[:, 0], minlength=filt.shape[0])
        interactions += numpy.bincount(data[:, 1], minlength=filt.shape[0])
        rev_mapping = numpy.where(interactions > 0)[0]
        mapping = numpy.zeros(filt.shape[0], dtype=numpy.int32) - 1
        mapping[rev_mapping] = numpy.arange(rev_mapping.shape[0])
        data[:, 0] = mapping[data[:, 0]]
        data[:, 1] = mapping[data[:, 1]]
        if options.binary:
            temp = BNEWT_sparse_binary(data, fl=(not options.silent))
        else:
            temp = BNEWT_sparse(data, fl=(not options.silent))
        corrections = 1.0 / temp[0].astype(numpy.float64)
        corr_mean = numpy.sum(corrections)
        corr_mean = corr_mean ** 2.0 - numpy.sum(corrections ** 2.0)
        corr_mean /= corrections.shape[0] * (corrections.shape[0] - 1)
        corrections /= corr_mean ** 0.5
        hic.corrections = numpy.ones(hic.filter.shape[0], dtype=numpy.float32)
        hic.corrections[rev_mapping] = corrections
        hic.normalization = 'express'
        hic.save(args[1])
    else:
        node_ranges = numpy.round(numpy.linspace(0, len(chroms), num_procs + 1)).astype(numpy.int32)
        if rank == 0:
            needed = chroms[:node_ranges[1]]
            for i in range(1, num_procs):
                comm.send(chroms[node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
        else:
            needed = comm.recv(source=0, tag=11)
        corrections = numpy.ones(hic.filter.shape[0], dtype=numpy.float32)
        for chrom in needed:
            chrint = hic.chr2int[chrom]
            start_fend = chr_indices[chrint]
            stop_fend = chr_indices[chrint + 1]
            start_index = hic.data['cis_indices'][start_fend]
            stop_index = hic.data['cis_indices'][stop_fend]
            data = hic.data['cis_data'][start_index:stop_index, :]
            data = data[numpy.where(hic.filter[data[:, 0]] * hic.filter[data[:, 1]])[0], :]
            rev_mapping = numpy.where(hic.filter[chr_indices[chrint]:chr_indices[chrint + 1]])[0]
            mapping = numpy.zeros(chr_indices[chrint + 1] - chr_indices[chrint], dtype=numpy.int32) - 1
            mapping[rev_mapping] = numpy.arange(rev_mapping.shape[0])
            data[:, 0] = mapping[data[:, 0] - chr_indices[chrint]]
            data[:, 1] = mapping[data[:, 1] - chr_indices[chrint]]
            if options.binary:
                temp = BNEWT_sparse_binary(data, fl=(not options.silent))
            else:
                temp = BNEWT_sparse(data, fl=(not options.silent))
            new_corrections = 1.0 / temp[0].astype(numpy.float64)
            corr_mean = numpy.sum(new_corrections)
            corr_mean = corr_mean ** 2.0 - numpy.sum(new_corrections ** 2.0)
            corr_mean /= new_corrections.shape[0] * (new_corrections.shape[0] - 1)
            new_corrections /= corr_mean ** 0.5
            corrections[chr_indices[chrint] + rev_mapping] = new_corrections
        if rank == 0:
            for i in range(1, num_procs):
                for j in range(node_ranges[i], node_ranges[i + 1]):
                    chrom = chroms[j]
                    chrint = hic.chr2int[chrom]
                    comm.Recv(corrections[chr_indices[chrint]:chr_indices[chrint + 1]], source=i, tag=13)
            hic.corrections = corrections
            hic.normalization = 'express'
            hic.save(args[1])
        else:
            for chrom in needed:
                chrint = hic.chr2int[chrom]
                comm.Send(corrections[chr_indices[chrint]:chr_indices[chrint + 1]], dest=0, tag=13)
    return None


if "mpi4py" in sys.modules.keys():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    num_procs = comm.Get_size()

if __name__ == "__main__":
    main()
