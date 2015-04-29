#!/usr/bin/python

import sys
import os

import h5py
import numpy as np
from math import ceil, floor
#import rpy2.robjects as robjects
#from rpy2.robjects.packages import importr
import statsmodels.api as sm

import hifive
from normalization_functions import *


#stats = importr('stats')
#base = importr('base')


def HiCNorm_Analysis():
    # get user input
    fend_fname, data_fname, out_fname, max_insert, bin_size, trans = sys.argv[1:7]
    data = h5py.File(data_fname, 'r')
    fendfilename = data['/'].attrs['fendfilename']
    data.close()
    if fendfilename[:2] == './':
        fendfilename = fendfilename[2:]
    parent_count = fendfilename.count('../')
    fendfilename = '/'.join(os.path.abspath(data_fname).split('/')[:-(1 + parent_count)] +
                                fendfilename.lstrip('/').split('/')[parent_count:])
    
    fends = h5py.File(fendfilename, 'r')
    all_chroms = list(fends['chromosomes'])
    fends.close()
    chr2int = {}
    for i, chrom in enumerate(all_chroms):
        chr2int[chrom] = i
    if len(sys.argv) >= 8:
        chroms = sys.argv[7].split(',')
    else:
        chroms = all_chroms
    if trans in ['1', 'true', 'True', 'TRUE']:
        trans = True
    else:
        trans = False
    # get data from hifive data hdf5 file
    cis_data, cis_indices, trans_data, trans_indices = Load_Data(fendfilename, data_fname, chroms)
    valid_chr = []
    if bin_size.count('M') > 0 or bin_size.count('m') > 0:
        bin_size = int(bin_size.strip('M').strip('m')) * 1000000
    elif bin_size.count('K') > 0 or bin_size.count('k') > 0:
        bin_size = int(bin_size.strip('K').strip('k')) * 1000
    else:
        bin_size = int(bin_size)
    # get features from fend file
    features = Get_Features(fend_fname, cis_data.keys(), int(max_insert), bin_size)
    # open hdf5 file for output heatmap
    output = h5py.File(out_fname, 'w')
    output.attrs['resolution'] = bin_size
    output.create_dataset(name="feature_names", data=np.array(["effective_len","gc","map","count"]))
    starts = {}
    mapping = {}
    binned_features = {}
    positions = {}
    valid_chrints = []
    valid_chroms = []
    for chrom in chroms:
        # bin features and see if there is enough data, given the bin size
        if chrom not in features:
    	    continue
        temp = Bin_Features(features, chrom, bin_size, starts, mapping)
        if not temp is None:
            binned_features[chrom] = temp[0]
            positions[chrom] = temp[1]
            output.create_dataset(name=("%s.features" % chrom), data=binned_features[chrom])
            output.create_dataset(name=("%s.positions" % chrom), data=positions[chrom])
            valid_chrints.append(chr2int[chrom])
            valid_chroms.append(chrom)
    valid_chrints.sort()
    output.create_dataset(name="coefficient_names", data=np.array(["intercept",'effective_len','gc']))
    # bin counts and find normalization
    for chrint in valid_chrints:
        chrom = all_chroms[chrint]
        binned_counts = Bin_Counts(cis_data, cis_indices, features, mapping, starts, chrom, bin_size)
        coefficients, expected_counts = Find_Expected_Counts2(binned_counts, binned_features[chrom], chrom)
        output.create_dataset(name=("%s.counts" % chrom), data=binned_counts)
        output.create_dataset(name=("%s.expected" % chrom), data=expected_counts)
        output.create_dataset(name=("%s.coefficients" % chrom), data=coefficients)
    # if needed, bin trans counts and find trans normalizations
    if trans:
        for i in range(len(valid_chrints)-1):
            chr1 = all_chroms[valid_chrints[i]]
            for j in range(i+1,len(valid_chrints)):
                chr2 = all_chroms[valid_chrints[j]]
                binned_counts = Bin_Trans_Counts(trans_data, trans_indices, features, mapping, starts, chr1, chr2,
                                                 bin_size)
                coefficients, expected_counts = Find_Expected_Trans_Counts2(binned_counts, binned_features[chr1],
                                                                           binned_features[chr2], chr1, chr2)
                output.create_dataset(name=("%s_by_%s.counts" % (chr1, chr2)), data=binned_counts)
                output.create_dataset(name=("%s_by_%s.expected" % (chr1, chr2)), data=expected_counts)
                output.create_dataset(name=("%s_by_%s.coefficients" % (chr1, chr2)), data=coefficients)
    output.create_dataset(name="chromosomes", data=np.array(valid_chroms))
    output.close()
    print >> sys.stderr, ("Done\n"),
    return None


def Load_Data(fend_fname, fname, chroms):
    print >> sys.stderr, ("Loading data..."),
    data = h5py.File(fname, 'r')
    fends = h5py.File(fend_fname, 'r')
    chromosomes = list(fends['chromosomes'][...])
    if len(chroms) == 0:
        chroms = chromosomes
    chr2int = {}
    chrints = []
    for i in range(len(chroms)):
        chr2int[chroms[i]] = chromosomes.index(chroms[i])
        chrints.append(chromosomes.index(chroms[i]))
    chrints.sort()
    chr_indices = fends['chr_indices'][...]
    cis_indices = data['cis_indices'][...]
    trans_indices = data['trans_indices'][...]
    all_data = data['cis_data'][...]
    cis_data = {}
    cis_data_indices = {}
    for i in range(len(chrints)):
        chrint = chrints[i]
        chrom = chromosomes[chrint]
        start_index = cis_indices[chr_indices[chrint]]
        stop_index = cis_indices[chr_indices[chrint + 1]]
        cis_data[chrom] = np.copy(all_data[start_index:stop_index, :])
        cis_data[chrom][:, :2] -= chr_indices[chrint]
        cis_data_indices[chrom] = np.copy(cis_indices[chr_indices[chrint]:(chr_indices[chrint + 1] + 1)])
        cis_data_indices[chrom] -= cis_data_indices[chrom][0]
    del all_data
    trans_data = {}
    trans_data_indices = {}
    all_data = data['trans_data'][...]
    for i in range(len(chrints) - 1):
        chrint1 = chrints[i]
        chrom1 = chromosomes[chrint1]
        start_index = trans_indices[chr_indices[chrint1]]
        stop_index = trans_indices[chr_indices[chrint1 + 1]]
        trans_chroms = fends['fends']['chr'][all_data[start_index:stop_index, 1]]
        for j in range(i + 1, len(chrints)):
            chrint2 = chrints[j]
            chrom2 = chromosomes[chrint2]
            trans_data[(chrom1, chrom2)] = np.copy(all_data[(np.where(trans_chroms == chrint2)[0] + start_index), :])
            trans_data[(chrom1, chrom2)][:, 0] -= chr_indices[chrint1]
            trans_data[(chrom1, chrom2)][:, 1] -= chr_indices[chrint2]
            trans_data_indices[(chrom1, chrom2)] = np.r_[0, np.bincount(trans_data[(chrom1, chrom2)][:, 0], minlength=(
                                                         chr_indices[chrint1 + 1] - 
                                                         chr_indices[chrint1]))].astype(np.int64)
            for k in range(1, trans_data_indices[(chrom1, chrom2)].shape[0]):
                trans_data_indices[(chrom1, chrom2)][k] += trans_data_indices[(chrom1, chrom2)][k - 1]
    del all_data
    data.close()
    fends.close()
    print >> sys.stderr, ("Done\n"),
    return [cis_data, cis_data_indices, trans_data, trans_data_indices]


def Get_Features(fname, chroms, max_insert, bin_size):
    print >> sys.stderr, ("Loading features..."),
    input = open(fname, 'r')
    input.readline()
    temp = input.readline()
    features = {}
    for chrom in chroms:
        features[chrom] = []
    while len(temp) > 0:
        fend_num, frag_num, chrom, start, valid, length, gc, mappability, map_bin = temp.strip('\n').split('\t')
        chrom = chrom.strip('chr')
        if chrom in features:
            fend_num = int(fend_num) - 1
            if fend_num % 2 == 0:
                cut = int(start)
            else:
                cut = int(start) + int(length)
            length = min(int(length), max_insert)
            features[chrom].append([cut, length, float(gc), float(mappability)])
        temp = input.readline()
    input.close()
    for chrom in features:
        features[chrom] = np.asarray( features[chrom], dtype=np.float32 )
    print >> sys.stderr, ("Done\n"),
    return features


def Bin_Features(features, chrom, bin_size, starts, mapping):
    start = (int(np.amin(features[chrom][:,0])) / bin_size) * bin_size
    starts[chrom] = start
    stop = (int(np.amax(features[chrom][:,0])) / bin_size) * bin_size
    num_bins = (stop - start) / bin_size + 1
    if num_bins < 5:
        return  None
    binned_features = np.zeros((num_bins,4), dtype=np.float32)
    indices = (features[chrom][:, 0].astype(np.int32) - start) / bin_size
    for i in range(1,4):
        binned_features[:, i - 1] = np.bincount(indices, weights=features[chrom][:, i], minlength=num_bins)
    binned_features[:, 3] = np.bincount(indices, minlength=num_bins)
    bin_positions = np.zeros((num_bins,2), dtype=np.int32)
    bin_positions[:,0] += np.arange(num_bins, dtype=np.int32) * bin_size + start
    bin_positions[:,1] += bin_positions[:,0] + bin_size
    mapping[chrom] = np.zeros(num_bins, dtype=np.int32)-1
    where = np.where(numpy.amin(binned_features, axis=1) > 0)[0]
    if len(where) < 5:
        return  None
    mapping[chrom][where] = np.arange(len(where), dtype=np.int32)
    binned_features = binned_features[where,:]
    binned_features[:,1:3] /= binned_features[:,3:]
    bin_positions = bin_positions[where,:]
    return [binned_features, bin_positions]


def Bin_Counts(data, indices, features, mapping, starts, chrom, binsize):
    print >> sys.stderr, ("Binning %s counts...") % (chrom),
    num_bins = np.amax(mapping[chrom]) + 1
    binned_counts = np.zeros(num_bins * (num_bins - 1) / 2, dtype=np.int32)
    BinCounts(
        binned_counts,
        data[chrom],
        indices[chrom],
        features[chrom],
        mapping[chrom],
        starts[chrom],
        binsize,
        num_bins)
    print >> sys.stderr, ("Done\n"),
    return binned_counts
    

def Bin_Trans_Counts(data, indices, features, mapping, starts, chrom1, chrom2, bin_size):
    print >> sys.stderr, ("Binning %s by %s counts...") % (chrom1, chrom2),
    num_bins1 = np.amax(mapping[chrom1]) + 1
    num_bins2 = np.amax(mapping[chrom2]) + 1
    binned_counts = np.zeros((num_bins1, num_bins2), dtype=np.int32)
    BinTransCounts(
        binned_counts,
        data[(chrom1, chrom2)],
        indices[(chrom1, chrom2)],
        features[chrom1],
        features[chrom2],
        mapping[chrom1],
        mapping[chrom2],
        starts[chrom1],
        starts[chrom2],
        bin_size)
    print >> sys.stderr, ("Done\n"),
    return binned_counts


def Find_Expected_Counts(counts, features, chrom):
    print >> sys.stderr, ("Arranging %s features...") % (chrom),
    count_vec = robjects.FloatVector( counts )
    np_len_vec = np.zeros(counts.shape[0], dtype=np.float32)
    np_gc_vec = np.zeros(counts.shape[0], dtype=np.float32)
    np_map_vec = np.zeros(counts.shape[0], dtype=np.float32)
    FindUpperTriFeatures(np_len_vec, features[:,0])
    np_len_vec -= np.mean(np_len_vec)
    np_len_vec /= np.std(np_len_vec)
    FindUpperTriFeatures(np_gc_vec, features[:,1])
    np_gc_vec -= np.mean(np_gc_vec)
    np_gc_vec /= np.std(np_gc_vec)
    FindUpperTriFeatures(np_map_vec, features[:,2])
    print >> sys.stderr, ("Done\nModeling %s counts...") % (chrom),
    len_vec = robjects.FloatVector(np_len_vec)
    gc_vec = robjects.FloatVector(np_gc_vec)
    map_vec = robjects.FloatVector(np_map_vec)
    robjects.globalenv['count'] = count_vec
    robjects.globalenv['len'] = len_vec
    robjects.globalenv['GCC'] = gc_vec
    model = stats.glm("count ~ len + GCC", offset=map_vec, family=robjects.r('poisson(link="log")'))
    coefficients = np.zeros(3, dtype=np.float32)
    for i in range(3):
        coefficients[i] = model.rx2('coefficients')[i]
    expected_counts = np.zeros(counts.shape[0], dtype=np.float32)
    print >> sys.stderr, ("Done\nFinding %s expected counts...") % (chrom),
    FindExpectedCounts(
        np_len_vec,
        np_gc_vec,
        np_map_vec,
        expected_counts,
        coefficients)
    print >> sys.stderr, ("Done\n"),
    return [coefficients, expected_counts]


def Find_Expected_Counts2(counts, features, chrom):
    print >> sys.stderr, ("Arranging %s features...") % (chrom),
    np_len_vec = np.zeros(counts.shape[0], dtype=np.float32)
    np_gc_vec = np.zeros(counts.shape[0], dtype=np.float32)
    np_map_vec = np.zeros(counts.shape[0], dtype=np.float32)
    FindUpperTriFeatures(np_len_vec, features[:,0])
    np_len_vec -= np.mean(np_len_vec)
    np_len_vec /= np.std(np_len_vec)
    FindUpperTriFeatures(np_gc_vec, features[:,1])
    np_gc_vec -= np.mean(np_gc_vec)
    np_gc_vec /= np.std(np_gc_vec)
    FindUpperTriFeatures(np_map_vec, features[:,2])
    independent = numpy.hstack((numpy.ones((np_gc_vec.shape[0], 1), dtype=numpy.float32), np_len_vec.reshape(-1, 1),
                                np_gc_vec.reshape(-1, 1)))
    print >> sys.stderr, ("Done\nModeling %s counts...") % (chrom),
    glm_model = sm.GLM(counts, independent, offset=np_map_vec, family=sm.families.Poisson(sm.families.links.log))
    coefficients = glm_model.fit().params.astype(numpy.float32)
    expected_counts = np.zeros(counts.shape[0], dtype=np.float32)
    print >> sys.stderr, ("Done\nFinding %s expected counts...") % (chrom),
    FindExpectedCounts(
        np_len_vec,
        np_gc_vec,
        np_map_vec,
        expected_counts,
        coefficients)
    print >> sys.stderr, ("Done\n"),
    return [coefficients, expected_counts]


def Find_Expected_Trans_Counts(counts, features1, features2, chrom1, chrom2):
    print >> sys.stderr, ("Arranging %s by %s features...") % (chrom1, chrom2),
    count_vec = robjects.FloatVector( counts.ravel( order='F'))
    len_mat = np.log(features1[:, 0].reshape(-1, 1) * features2[:, 0].reshape(1, -1))
    len_mat -= np.mean(len_mat)
    len_mat /= np.std(len_mat)
    gc_mat = np.log(features1[:,1].reshape(-1, 1) * features2[:,1].reshape(1, -1))
    gc_mat -= np.mean(gc_mat)
    gc_mat /= np.std(gc_mat)
    map_mat = np.log(features1[:,2].reshape(-1, 1) * features2[:,2].reshape(1, -1))
    print >> sys.stderr, ("Done\nModeling %s by %s counts...") % (chrom1, chrom2),
    len_vec = robjects.FloatVector(len_mat.ravel( order='F'))
    gc_vec = robjects.FloatVector(gc_mat.ravel( order='F'))
    map_vec = robjects.FloatVector( map_mat.ravel( order='F'))
    robjects.globalenv['count'] = count_vec
    robjects.globalenv['len'] = len_vec
    robjects.globalenv['GCC'] = gc_vec
    model = stats.glm(formula="count ~ len + GCC", offset=map_vec, family=robjects.r('poisson(link="log")'))
    coefficients = np.zeros(3, dtype=np.float32)
    for i in range(3):
        coefficients[i] = model.rx2('coefficients')[i]
    print >> sys.stderr, ("Done\nFinding %s by %s expected counts...") % (chrom1, chrom2),
    expected_counts = np.zeros(len_mat.shape, dtype=np.float32)
    FindTransExpectedCounts(
        len_mat,
        gc_mat,
        map_mat,
        expected_counts,
        coefficients)
    print >> sys.stderr, ("Done\n"),
    return [coefficients, expected_counts]


def Find_Expected_Trans_Counts2(counts, features1, features2, chrom1, chrom2):
    print >> sys.stderr, ("Arranging %s by %s features...") % (chrom1, chrom2),
    len_mat = np.log(features1[:, 0].reshape(-1, 1) * features2[:, 0].reshape(1, -1))
    len_mat -= np.mean(len_mat)
    len_mat /= np.std(len_mat)
    gc_mat = np.log(features1[:,1].reshape(-1, 1) * features2[:,1].reshape(1, -1))
    gc_mat -= np.mean(gc_mat)
    gc_mat /= np.std(gc_mat)
    map_mat = np.log(features1[:,2].reshape(-1, 1) * features2[:,2].reshape(1, -1))
    independent = numpy.hstack((numpy.ones((gc_mat.shape[0] * gc_mat.shape[1], 1), dtype=numpy.float32),
                                len_mat.ravel().reshape(-1, 1), gc_mat.ravel().reshape(-1, 1)))
    print >> sys.stderr, ("Done\nModeling %s by %s counts...") % (chrom1, chrom2),
    glm_model = sm.GLM(counts.ravel(), independent, offset=map_mat.ravel(),
                       family=sm.families.Poisson(sm.families.links.log))
    coefficients = glm_model.fit().params.astype(numpy.float32)
    expected_counts = np.zeros(len_mat.shape, dtype=np.float32)
    print >> sys.stderr, ("Done\nFinding %s by %s expected counts...") % (chrom1, chrom2),
    FindTransExpectedCounts(
        len_mat,
        gc_mat,
        map_mat,
        expected_counts,
        coefficients)
    print >> sys.stderr, ("Done\n"),
    return [coefficients, expected_counts]


if __name__ == '__main__':
    usage =    """Usage: normalization.py <fend_fname> <data.hdf5> <results.hdf5> <max_insert> <bin_size> <trans_flag> [chr_1,chr_2,...chr_n]"""
    if len(sys.argv) < 7:
        print usage
    else:
        HiCNorm_Analysis()
