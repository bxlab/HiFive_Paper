#!/usr/bin/env python

import sys

from hiclib import mapping, fragmentHiC
from mirnylib import h5dict, genome

basedir = sys.argv[1]

mapped_reads1 = h5dict.h5dict('%s/Data/Timing/mapped_reads1.hdf5' % basedir)
mapped_reads2 = h5dict.h5dict('%s/Data/Timing/mapped_reads2.hdf5' % basedir)
mapped_reads3 = h5dict.h5dict('%s/Data/Timing/mapped_reads3.hdf5' % basedir)
genome_db    = genome.Genome('%s/Data/Genome/mm9_fasta' % basedir, readChrms=['1'], chrmFileTemplate="%s.fa")

mapping.parse_sam(
    sam_basename1='%s/Data/Timing/SRR443886_sub_1.bam' % basedir,
    sam_basename2='%s/Data/Timing/SRR443886_sub_2.bam' % basedir,
    out_dict=mapped_reads1,
    genome_db=genome_db, 
    enzyme_name='NcoI')

mapping.parse_sam(
    sam_basename1='%s/Data/Timing/SRR443887_sub_1.bam' % basedir,
    sam_basename2='%s/Data/Timing/SRR443887_sub_2.bam' % basedir,
    out_dict=mapped_reads2,
    genome_db=genome_db, 
    enzyme_name='NcoI')

mapping.parse_sam(
    sam_basename1='%s/Data/Timing/SRR443888_sub_1.bam' % basedir,
    sam_basename2='%s/Data/Timing/SRR443888_sub_2.bam' % basedir,
    out_dict=mapped_reads3,
    genome_db=genome_db, 
    enzyme_name='NcoI')

fragments2 = fragmentHiC.HiCdataset(
    filename='temp2',
    genome=genome_db,
    maximumMoleculeLength=500,
    mode='w',
    enzymeName="NcoI",
    inMemory=True)
fragments2.parseInputData(dictLike="%s/Data/Timing/mapped_reads2.hdf5" % basedir)
fragments2.save('%s/Data/Timing/hiclib_data2.hdf5' % basedir)
del fragments2

fragments3 = fragmentHiC.HiCdataset(
    filename='temp3',
    genome=genome_db,
    maximumMoleculeLength=500,
    mode='w',
    enzymeName="NcoI",
    inMemory=True)
fragments3.parseInputData(dictLike="%s/Data/Timing/mapped_reads3.hdf5" % basedir)
fragments3.save('%s/Data/Timing/hiclib_data3.hdf5' % basedir)
del fragments3

fragments = fragmentHiC.HiCdataset(
    filename='temp1',
    genome=genome_db,
    maximumMoleculeLength=500,
    mode='w',
    enzymeName="NcoI",
    inMemory=True)
fragments.parseInputData(dictLike="%s/Data/Timing/mapped_reads1.hdf5" % basedir)
fragments.merge(['%s/Data/Timing/hiclib_data2.hdf5' % basedir, '%s/Data/Timing/hiclib_data3.hdf5' % basedir])
fragments.save('%s/Data/Timing/hiclib_data.hdf5' % basedir)
