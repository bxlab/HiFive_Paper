#!/usr/bin/env python

"""
This is a module contains scripts for generating compact, upper-triangle and full matrices of HiC interaction data.

Concepts
--------

These functions rely on the :class:`HiC` class in conjunction with the :class:`Fend` and :class:`HiCData` classes.

Data can either be arranged in compact, complete, or flattened (row-major) upper-triangle arrays. Compact arrays are N x M, where N is the number of fends or bins, and M is the maximum distance between fends or bins. This is useful for working with sets of short interactions. Data can be raw, fend-corrected, distance-dependence removed, or enrichment values. Arrays are 3-dimensional with observed values in the first layer of d3, expected values in the second layer of d3. The exception to this is upper-triangle arrays, which are 2d, divinding observed and expected along the second axis.

API Documentation
-----------------
"""

import os
import sys
import subprocess

import numpy
import h5py
try:
    from mpi4py import MPI
except:
    pass

from libraries._hic_interactions import find_max_fend
import libraries._hic_binning as _hic_binning


def find_cis_signal(hic, chrom, binsize=10000, binbounds=None, start=None, stop=None, startfend=None, stopfend=None,
                    datatype='enrichment', arraytype='compact', maxdistance=0, skipfiltered=False, returnmapping=False,
                    proportional=False, **kwargs):
    """
    Create an array of format 'arraytype' and fill with data requested in 'datatype'.

    :param hic: A :class:`HiC <hifive.hic.HiC>` class object containing fend and count data.
    :type hic: :class:`HiC <hifive.hic.HiC>`
    :param chrom: The name of a chromosome contained in 'hic'.
    :type chrom: str.
    :param binsize: This is the coordinate width of each bin. A value of zero indicates unbinned. If binbounds is not None, this value is ignored.
    :type binsize: int.
    :param binbounds: An array containing start and stop coordinates for a set of user-defined bins. Any fend not falling in a bin is ignored.
    :type binbounds: numpy array
    :param start: The smallest coordinate to include in the array, measured from fend midpoints or the start of the first bin. If 'binbounds' is given, this value is ignored. If both 'start' and 'startfend' are given, 'start' will override 'startfend'. If unspecified, this will be set to the midpoint of the first fend for 'chrom', adjusted to the first multiple of 'binsize' if not zero. Optional.
    :type start: int.
    :param stop: The largest coordinate to include in the array, measured from fend midpoints or the end of the last bin. If 'binbounds' is given, this value is ignored. If both 'stop' and 'stopfend' are given, 'stop' will override 'stopfend'. If unspecified, this will be set to the midpoint of the last fend plus one for 'chrom', adjusted to the last multiple of 'start' + 'binsize' if not zero. Optional.
    :type stop: int.
    :param startfend: The first fend to include in the array. If 'binbounds' is given, this value is ignored. If unspecified and 'start' is not given, this is set to the first valid fend in 'chrom'. In cases where 'start' is specified and conflicts with 'startfend', 'start' is given preference. Optional
    :type startfend: int.
    :param stopfend: The first fend not to include in the array. If 'binbounds' is given, this value is ignored. If unspecified and 'stop' is not given, this is set to the last valid fend in 'chrom' plus one. In cases where 'stop' is specified and conflicts with 'stopfend', 'stop' is given preference. Optional.
    :type stopfend: str.
    :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend', 'enrichment', and 'expected'. Observed values are always in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, unfiltered fends return value of one. Expected values are returned for 'distance', 'fend', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fend' uses only fend correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
    :type datatype: str.
    :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'compact', 'full', and 'upper'. 'compact' means data are arranged in a N x M x 2 array where N is the number of bins, M is the maximum number of steps between included bin pairs, and data are stored such that bin n,m contains the interaction values between n and n + m + 1. 'full' returns a square, symmetric array of size N x N x 2. 'upper' returns only the flattened upper triangle of a full array, excluding the diagonal of size (N * (N - 1) / 2) x 2.
    :type arraytype: str.
    :param maxdistance: This specifies the maximum coordinate distance between bins that will be included in the array. If set to zero, all distances are included.
    :type maxdistance: str.
    :param skipfiltered: If 'True', all interaction bins for filtered out fends are removed and a reduced-size array is returned.
    :type skipfiltered: bool.
    :param returnmapping: If 'True', a list containing the data array and a 2d array containing first coordinate included and excluded from each bin, and the first fend included and excluded from each bin is returned. Otherwise only the data array is returned.
    :type returnmapping: bool.
    :param proportional: Indicates whether interactions should proportionally contribute to bins based on the amount of overlap instead of being attributed solely based on midpoint. Only valid for binned heatmaps.
    :type proportional: bool.
    :returns: Array in format requested with 'arraytype' containing data requested with 'datatype'.
    """
    if 'silent' in kwargs:
        silent = kwargs['silent']
    else:
        silent = False
    # check that all values are acceptable
    if datatype not in ['raw', 'fend', 'distance', 'enrichment', 'expected']:
        if not silent:
            print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
        return None
    elif datatype in ['fend', 'enrichment'] and hic.normalization == 'none':
        if not silent:
            print >> sys.stderr, ("Normalization has not been performed yet on this project. Select either 'raw' or 'distance' for datatype. No data returned\n"),
        return None
    elif datatype in ['distance', 'enrichment'] and hic.distance_parameters is None:
        if not silent:
            print >> sys.stderr, ("No distance-dependence relationship has been calculated for this project yet. Select either 'raw' or 'fend' for datatype. No data returned\n"),
        return None
    if arraytype not in ['full', 'compact', 'upper']:
        if not silent:
            print >> sys.stderr, ("Unrecognized array type. No data returned.\n"),
        return None
    # Determine start, stop, startfend, and stopfend
    chrint = hic.chr2int[chrom.strip('chr')]
    if not binbounds is None:
        start = binbounds[0, 0]
        stop = binbounds[-1, 1]
        startfend = _find_fend_from_coord(hic, chrint, start)
        stopfend = _find_fend_from_coord(hic, chrint, stop) + 1
    else:
        if start is None and startfend is None:
            startfend = hic.fends['chr_indices'][chrint]
            while startfend < hic.fends['chr_indices'][chrint + 1] and hic.filter[startfend] == 0:
                startfend += 1
            if startfend == hic.fends['chr_indices'][chrint + 1]:
                if not silent:
                    print >> sys.stderr, ("Insufficient data.\n"),
                return None
            start = hic.fends['fends']['mid'][startfend]
            if binsize > 0:
                start = (start / binsize) * binsize
        elif start is None:
            start = hic.fends['fends']['mid'][startfend]
            if binsize > 0:
                start = (start / binsize) * binsize
        else:
            startfend = _find_fend_from_coord(hic, chrint, start)
        if (stop is None or stop == 0) and stopfend is None:
            stopfend = hic.fends['chr_indices'][chrint + 1]
            while stopfend > hic.fends['chr_indices'][chrint] and hic.filter[stopfend - 1] == 0:
                stopfend -= 1
            stop = hic.fends['fends']['mid'][stopfend - 1]
            if binsize > 0:
                stop = ((stop - 1 - start) / binsize + 1) * binsize + start
        elif stop is None or stop == 0:
            stop = hic.fends['fends']['mid'][stopfend - 1]
            if binsize > 0:
                stop = ((stop - 1 - start) / binsize + 1) * binsize + start
        else:
            if binsize > 0:
                stop = ((stop - 1 - start) / binsize + 1) * binsize + start
            stopfend = _find_fend_from_coord(hic, chrint, stop) + 1
    if not silent:
        print >> sys.stderr, ("Finding %s %s array for %s:%i-%i...") % (datatype, arraytype, chrom, start, stop),
    # If datatype is not 'expected', pull the needed slice of data
    if datatype != 'expected':
        start_index = hic.data['cis_indices'][startfend]
        stop_index = hic.data['cis_indices'][stopfend]
        if start_index == stop_index:
            if not silent:
                print >> sys.stderr, ("Insufficient data\n"),
            return None
        data_indices = hic.data['cis_indices'][startfend:(stopfend + 1)]
        data_indices -= data_indices[0]
        data = hic.data['cis_data'][start_index:stop_index, :]
        data[:, :2] -= startfend
    else:
        data_indices = None
        data = None
    # Determine mapping of valid fends to bins
    mapping = numpy.zeros(stopfend - startfend, dtype=numpy.int32) - 1
    valid = numpy.where(hic.filter[startfend:stopfend] > 0)[0]
    mids = hic.fends['fends']['mid'][startfend:stopfend]
    if binsize == 0 and binbounds is None:
        if skipfiltered:
            mapping[valid] = numpy.arange(valid.shape[0])
            num_bins = valid.shape[0]
        else:
            mapping[valid] = valid
            num_bins = mapping.shape[0]
    elif not binbounds is None:
        start_indices = numpy.searchsorted(binbounds[:, 0], mids[valid], side='right') - 1
        stop_indices = numpy.searchsorted(binbounds[:, 1], mids[valid], side='right')
        where = numpy.where(start_indices == stop_indices)[0]
        valid = valid[where]
        mapping[valid] = start_indices[where]
        num_bins = binbounds.shape[0]
    else:
        mapping[valid] = (mids[valid] - start) / binsize
        num_bins = (stop - start) / binsize
    # Find maximum interaction partner for each fend
    if num_bins < 2:
        if not silent:
            print >> sys.stderr, ("Insufficient data\n"),
        return None
    max_fend = numpy.zeros(mapping.shape[0], dtype=numpy.int32)
    find_max_fend(max_fend, mids, hic.fends['fends']['chr'][startfend:stopfend],
                  hic.fends['chr_indices'][...], startfend, maxdistance)
    max_fend = numpy.minimum(max_fend, mapping.shape[0])
    if binsize == 0:
        max_bin = numpy.amax(max_fend - numpy.arange(mapping.shape[0]))
        if max_bin <= 0:
            if not silent:
                print >> sys.stderr, ("Insufficient data.\n"),
            return None
    else:
        if maxdistance == 0:
            max_bin = num_bins - 1
        else:
            max_bin = maxdistance / binsize
    # If correction is required, determine what type and get appropriate data
    if 'binning' not in hic.normalization and datatype != 'raw':
        corrections = hic.corrections[startfend:stopfend]
    elif datatype == 'raw':
        corrections = numpy.ones(stopfend - startfend, dtype=numpy.float32)
    else:
        corrections = None
    if ((hic.normalization in ['express', 'probability'] and
            datatype == 'fend') or datatype == 'raw') and maxdistance == 0:
        if datatype == 'fend':
            correction_sums = numpy.bincount(mapping[valid], weights=corrections[valid],
                                             minlength=num_bins).astype(numpy.float64)
        else:
            correction_sums = numpy.bincount(mapping[valid], minlength=num_bins).astype(numpy.float64)
    else:
        correction_sums = None
    if 'binning' in hic.normalization and datatype not in ['raw', 'distance']:
        binning_corrections = hic.binning_corrections
        binning_num_bins = hic.binning_num_bins
        fend_indices = hic.binning_fend_indices
    else:
        binning_corrections = None
        binning_num_bins = None
        fend_indices = None
    if datatype in ['distance', 'enrichment', 'expected']:
        distance_parameters = hic.distance_parameters
        chrom_mean = hic.chromosome_means[chrint]
    else:
        distance_parameters = None
        chrom_mean = 0.0
    # If proportional is requested, find bin ranges
    if proportional and binsize > 0:
        fends = hic.fends['fends'][startfend:stopfend]
        ranges = numpy.zeros((mapping.shape[0], 2), dtype=numpy.int32)
        overlap = numpy.zeros(ranges.shape, dtype=numpy.float32)
        positions = numpy.arange(1, 1 + num_bins) * binsize + start
        ranges[:, 0] = numpy.searchsorted(positions, fends['start'])
        ranges[:, 1] = numpy.searchsorted(positions[:-1], fends['stop'])
        where = numpy.where(ranges[:, 0] < ranges[:, 1])[0]
        overlap[where, 0] = numpy.minimum(positions[ranges[where, 0]] - fends['start'][where],
                                          binsize) / float(binsize)
        overlap[where, 1] = numpy.minimum(fends['stop'][where] - positions[ranges[where, 1]] + binsize,
                                          binsize) / float(binsize)
        where = numpy.where(ranges[:, 0] == ranges[:, 1])[0]
        overlap[where, 0] = (fends['stop'][where] - fends['start'][where]) / float(binsize)
    else:
        ranges = None
        overlap = None
    # Create requested array
    if arraytype == 'compact':
        data_array = numpy.zeros((num_bins, max_bin, 2), dtype=numpy.float32)
    else:
        data_array = numpy.zeros((num_bins * (num_bins - 1) / 2, 2), dtype=numpy.float32)
    # Fill in data values
    if arraytype == 'compact':
        _hic_binning.find_cis_compact_expected(mapping, corrections, binning_corrections,
                                               binning_num_bins, fend_indices, mids, distance_parameters,
                                               max_fend, data_array, correction_sums, ranges, overlap,
                                               chrom_mean, startfend)
        if datatype != 'expected':
            _hic_binning.find_cis_compact_observed(data, data_indices, mapping, max_fend, data_array, ranges, overlap)
        else:
            data_array[:, :, 0] = data_array[:, :, 1]
            data_array[:, :, 1].fill(0)
            correction_sums = numpy.bincount(mapping[valid], minlength=num_bins).astype(numpy.float64)
            corrections.fill(1)
            _hic_binning.find_cis_compact_expected(mapping, corrections, None, None, None, mids, None,
                                                   max_fend, data_array, correction_sums, ranges, overlap,
                                                   chrom_mean, startfend)
            data_array = data_array[:, :, ::-1]
    else:
        _hic_binning.find_cis_upper_expected(mapping, corrections, binning_corrections,
                                             binning_num_bins, fend_indices, mids, distance_parameters,
                                             max_fend, data_array, correction_sums, ranges, overlap,
                                             chrom_mean, startfend)
        if datatype != 'expected':
            _hic_binning.find_cis_upper_observed(data, data_indices, mapping, max_fend, data_array, ranges, overlap)
        else:
            data_array[:, 0] = data_array[:, 1]
            data_array[:, 1].fill(0)
            correction_sums = numpy.bincount(mapping[valid], minlength=num_bins).astype(numpy.float64)
            corrections.fill(1)
            _hic_binning.find_cis_upper_expected(mapping, corrections, None, None, None, mids, None,
                                                 max_fend, data_array, correction_sums, ranges, overlap,
                                                 chrom_mean, startfend)
            data_array = data_array[:, ::-1]
    # If requesting 'full' array, convert 'upper' array type to 'full'
    if arraytype == 'full':
        indices = numpy.triu_indices(num_bins, 1)
        full_data_array = numpy.zeros((num_bins, num_bins, 2), dtype=numpy.float32)
        full_data_array[indices[1], indices[0], :] = data_array
        full_data_array[indices[0], indices[1], :] = data_array
        del data_array
        data_array = full_data_array
    if returnmapping:
        bin_mapping = numpy.zeros((num_bins, 4), dtype=numpy.int32)
        if binsize == 0 and binbounds is None:
            if skipfiltered:
                bin_mapping[:, 2] = valid + startfend
            else:
                bin_mapping[:, 2] = numpy.arange(startfend, stopfend)
            bin_mapping[:, 3] = bin_mapping[:, 2] + 1
            bin_mapping[:, 0] = hic.fends['fends']['start'][bin_mapping[:, 2]]
            bin_mapping[:, 1] = hic.fends['fends']['stop'][bin_mapping[:, 2]]
        else:
            if binbounds is None:
                bin_mapping[:, 0] = start + binsize * numpy.arange(num_bins)
                bin_mapping[:, 1] = bin_mapping[:, 0] + binsize
            else:
                bin_mapping[:, :2] = binbounds
            bin_mapping[:, 2] = numpy.searchsorted(mids, bin_mapping[:, 0]) + startfend
            bin_mapping[:, 3] = numpy.searchsorted(mids, bin_mapping[:, 1]) + startfend
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return [data_array, bin_mapping]
    else:
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return data_array

def _find_fend_from_coord(hic, chrint, coord):
    """Find the next fend after the coordinate on chromosome 'chrint'."""
    first_fend = hic.fends['chr_indices'][chrint]
    last_fend = hic.fends['chr_indices'][chrint + 1]
    return numpy.searchsorted(hic.fends['fends']['mid'][first_fend:last_fend], coord) + first_fend

def bin_cis_array(data_array, data_mapping, binsize=10000, binbounds=None, start=None, stop=None, arraytype='full',
                  returnmapping=False, **kwargs):
    """
    Create an array of format 'arraytype' and fill 'binsize' bins or bins defined by 'binbounds' with data provided in the array passed by 'data_array'.

    :param data_array: A 2d (upper) or 3d (compact) array containing data to be binned. Array format will be determined from the number of dimensions.
    :type data_array: numpy array
    :param data_mapping: An N x 4 2d integer array containing the start and stop coordinates, and start and stop fends for each of the N bin ranges in 'data_array'.
    :type data_mapping: numpy array
    :param binsize: This is the coordinate width of each bin. If binbounds is not None, this value is ignored.
    :type binsize: int.
    :param binbounds: An array containing start and stop coordinates for a set of user-defined bins. Any bin from 'data_array' not falling in a bin is ignored.
    :type binbounds: numpy array
    :param start: The coordinate at the beginning of the first bin of the binned data. If unspecified, 'start' will be the first multiple of 'binsize' below the first coordinate from 'data_mapping'. If 'binbounds' is given, 'start' is ignored. Optional.
    :type start: int.
    :param stop: The coordinate at the end of the last bin of the binned data. If unspecified, 'stop' will be the first multiple of 'binsize' after the last coordinate from 'data_mapping'. If needed, 'stop' is adjusted upward to create a complete last bin. If 'binbounds' is given, 'stop' is ignored. Optional.
    :type stop: int.
    :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'compact', 'full', and 'upper'. 'compact' means data are arranged in a N x M x 2 array where N is the number of bins, M is the maximum number of steps between included bin pairs, and data are stored such that bin n,m contains the interaction values between n and n + m + 1. 'full' returns a square, symmetric array of size N x N x 2. 'upper' returns only the flattened upper triangle of a full array, excluding the diagonal of size (N * (N - 1) / 2) x 2.
    :type arraytype: str.
    :param returnmapping: If 'True', a list containing the data array and a 2d array containing first coordinate included and excluded from each bin, and the first fend included and excluded from each bin is returned. Otherwise only the data array is returned.
    :type returnmapping: bool.
    :returns: Array in format requested with 'arraytype' containing binned data requested with 'datatype' pulled from 'data_array' or list of binned data array and mapping array.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    # check that arraytype value is acceptable
    if arraytype not in ['full', 'compact', 'upper']:
        if not silent:
            print >> sys.stderr, ("Unrecognized array type. No data returned.\n"),
        return None
    # Determine input array type
    if len(data_array.shape) == 2 and data_mapping.shape[0] * (data_mapping.shape[0] - 1) / 2 == data_array.shape[0]:
        input_type = 'upper'
    elif len(data_array.shape) == 3 and data_array.shape[0] == data_mapping.shape[0]:
        input_type = 'compact'
    else:
        if not silent:
            print >> sys.stderr, ("Unrecognized input array type. No data returned.\n"),
        return None
    # Determine start and stop, if necessary
    if binbounds is None:
        if start is None:
            start = (data_mapping[0, 0] / binsize) * binsize
        if stop is None:
            stop = ((data_mapping[-1, 1] - 1) / binsize + 1) * binsize
        else:
            stop = ((stop - 1 - start) / binsize + 1) * binsize + start
        num_bins = (stop - start) / binsize
        binbounds = numpy.zeros((num_bins, 2), dtype=numpy.int32)
        binbounds[:, 0] = numpy.arange(num_bins) * binsize + start
        binbounds[:, 1] = binbounds[:, 0] + binsize
    else:
        num_bins = binbounds.shape[0]
        start = binbounds[0, 0]
        stop = binbounds[0, 1]
    mids = (data_mapping[:, 0] + data_mapping[:, 1]) / 2
    if not silent:
        print >> sys.stderr, ("Finding binned %s array...") % (arraytype),
    # Find bin mapping for each fend
    mapping = numpy.zeros(mids.shape[0], dtype=numpy.int32) - 1
    fend_ranges = numpy.zeros((binbounds.shape[0], 2), dtype=numpy.int32)
    for i in range(binbounds.shape[0]):
        firstbin = numpy.searchsorted(mids, binbounds[i, 0])
        lastbin = numpy.searchsorted(mids, binbounds[i, 1])
        mapping[firstbin:lastbin] = i
        fend_ranges[i, 0] = data_mapping[firstbin, 2]
        fend_ranges[i, 1] = data_mapping[lastbin, 3]
    # Create requested array
    if arraytype == 'compact':
        max_bin = (stop - start) / binsize + 1
        binned_array = numpy.zeros((num_bins, max_bin, 2), dtype=numpy.float32)
    else:
        binned_array = numpy.zeros((num_bins * (num_bins - 1) / 2, 2), dtype=numpy.float32)
    # Fill in binned data values
    if arraytype == 'compact':
        if input_type == 'compact':
            _hic_binning.bin_compact_to_compact(binned_array, data_array, mapping)
        else:
            _hic_binning.bin_upper_to_compact(binned_array, data_array, mapping)
        # Trim unused bins
        valid = numpy.where(numpy.sum(binned_array[:, :, 1] > 0, axis=0) > 0)[0][-1]
        binned_array = binned_array[:, :(valid + 1), :]
    else:
        if input_type == 'compact':
            _hic_binning.bin_compact_to_upper(binned_array, data_array, mapping, num_bins)
        else:
            _hic_binning.bin_upper_to_upper(binned_array, data_array, mapping, num_bins)
    # If requesting 'full' array, convert 'upper' array type to 'full'
    if arraytype == 'full':
        indices = numpy.triu_indices(num_bins, 1)
        full_binned_array = numpy.zeros((num_bins, num_bins, 2), dtype=numpy.float32)
        full_binned_array[indices[1], indices[0], :] = binned_array
        full_binned_array[indices[0], indices[1], :] = binned_array
        del binned_array
        binned_array = full_binned_array
    # If mapping requested, calculate bin bounds
    if returnmapping:
        mapping = numpy.zeros((num_bins, 4), dtype=numpy.int32)
        mapping[:, 0] = binbounds[:, 0]
        mapping[:, 1] = binbounds[:, 1]
        mapping[:, 2:4] = fend_ranges
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return [binned_array, mapping]
    else:
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return binned_array

def dynamically_bin_cis_array(unbinned, unbinnedpositions, binned, binbounds, minobservations=10,
                              searchdistance=0, removefailed=True, **kwargs):
    """
    Expand bins in 'binned' to include additional data provided in 'unbinned' as necessary to meet 'minobservations', or 'searchdistance' criteria.

    :param unbinned: A 2d or 3d array containing data in either compact or upper format to be used for filling expanding bins. Array format will be determined from the number of dimensions.
    :type unbinned: numpy array
    :param unbinnedpositions: A 2d integer array indicating the first and last coordinate of each bin in 'unbinned' array.
    :type unbinnedpositions: numpy array
    :param binned: A 2d or 3d array containing binned data in either compact or upper format to be dynamically binned. Array format will be determined from the number of dimensions. Data in this array will be altered by this function.
    :type binned: numpy array
    :param binbounds: An integer array indicating the start and end position of each bin in 'binned' array. This array should be N x 2, where N is the number of intervals in 'binned'.
    :type binbounds: numpy array
    :param minobservations: The fewest number of observed reads needed for a bin to counted as valid and stop expanding.
    :type minobservations: int.
    :param searchdistance: The furthest distance from the bin minpoint to expand bounds. If this is set to zero, there is no limit on expansion distance.
    :type searchdistance: int.
    :param removefailed: If a non-zero 'searchdistance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'removefailed' is True, the observed and expected values for that bin are zero.
    :type removefailed: bool.
    :returns: None
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    # Determine unbinned array type
    if len(unbinned.shape) == 2 and (unbinnedpositions.shape[0] * (unbinnedpositions.shape[0] - 1) / 2 ==
                                     unbinned.shape[0]):
        unbinned_type = 'upper'
    elif len(unbinned.shape) == 3 and unbinned.shape[0] == unbinnedpositions.shape[0]:
        unbinned_type = 'compact'
    else:
        if not silent:
            print >> sys.stderr, ("Unrecognized unbinned array type. No data returned.\n"),
        return None
    # Determine binned array type
    if len(binned.shape) == 2 and binbounds.shape[0] * (binbounds.shape[0] - 1) / 2 == binned.shape[0]:
        binned_type = 'upper'
    elif len(binned.shape) == 3 and binned.shape[0] == binbounds.shape[0]:
        binned_type = 'compact'
    else:
        if not silent:
            print >> sys.stderr, ("Unrecognized binned array type. No data returned.\n"),
        return None
    if not silent:
        print >> sys.stderr, ("Dynamically binning data..."),
    # Determine bin edges relative to unbinned positions
    unbinnedmids = (unbinnedpositions[:, 0] + unbinnedpositions[:, 1]) / 2
    binedges = numpy.zeros(binbounds.shape, dtype=numpy.int32)
    binedges[:, 0] = numpy.searchsorted(unbinnedmids, binbounds[:, 0])
    binedges[:, 1] = numpy.searchsorted(unbinnedmids, binbounds[:, 1])
    # Determine bin midpoints
    mids = (binbounds[:, 0] + binbounds[:, 1]) / 2
    # Dynamically bin using appropriate array type combination
    if unbinned_type == 'upper':
        if binned_type == 'upper':
            _hic_binning.dynamically_bin_upper_from_upper(unbinned, unbinnedmids, binned, binedges,
                                                      mids, minobservations, searchdistance, int(removefailed))
        else:
            _hic_binning.dynamically_bin_compact_from_upper(unbinned, unbinnedmids, binned, binedges,
                                                        mids, minobservations, searchdistance, int(removefailed))
    else:
        if binned_type == 'upper':
            _hic_binning.dynamically_bin_upper_from_compact(unbinned, unbinnedmids, binned, binedges,
                                                        mids, minobservations, searchdistance, int(removefailed))
        else:
            _hic_binning.dynamically_bin_compact_from_compact(unbinned, unbinnedmids, binned, binedges,
                                                          mids, minobservations, searchdistance, int(removefailed))
    if not silent:
        print >> sys.stderr, ("Done\n"),
    return None

def find_trans_signal(hic, chrom1, chrom2, binsize=10000, binbounds1=None, binbounds2=None, start1=None, stop1=None,
                      startfend1=None, stopfend1=None, start2=None, stop2=None, startfend2=None, stopfend2=None,
                      datatype='enrichment', skipfiltered=False, returnmapping=False, **kwargs):
    """
    Create an array of format 'arraytype' and fill with data requested in 'datatype'.

    :param hic: A :class:`HiC <hifive.hic.HiC>` class object containing fend and count data.
    :type hic: :class:`HiC <hifive.hic.HiC>`
    :param chrom: The name of a chromosome contained in 'hic'.
    :type chrom: str.
    :param binsize: This is the coordinate width of each bin. A value of zero indicates unbinned. If binbounds is not None, this value is ignored.
    :type binsize: int.
    :param binbounds: An array containing start and stop coordinates for a set of user-defined bins. Any fend not falling in a bin is ignored.
    :type binbounds: numpy array
    :param start: The smallest coordinate to include in the array, measured from fend midpoints or the start of the first bin. If 'binbounds' is given, this value is ignored. If both 'start' and 'startfend' are given, 'start' will override 'startfend'. If unspecified, this will be set to the midpoint of the first fend for 'chrom', adjusted to the first multiple of 'binsize' if not zero. Optional.
    :type start: int.
    :param stop: The largest coordinate to include in the array, measured from fend midpoints or the end of the last bin. If 'binbounds' is given, this value is ignored. If both 'stop' and 'stopfend' are given, 'stop' will override 'stopfend'. If unspecified, this will be set to the midpoint of the last fend plus one for 'chrom', adjusted to the last multiple of 'start' + 'binsize' if not zero. Optional.
    :type stop: int.
    :param startfend: The first fend to include in the array. If 'binbounds' is given, this value is ignored. If unspecified and 'start' is not given, this is set to the first valid fend in 'chrom'. In cases where 'start' is specified and conflicts with 'startfend', 'start' is given preference. Optional
    :type startfend: int.
    :param stopfend: The first fend not to include in the array. If 'binbounds' is given, this value is ignored. If unspecified and 'stop' is not given, this is set to the last valid fend in 'chrom' plus one. In cases where 'stop' is specified and conflicts with 'stopfend', 'stop' is given preference. Optional.
    :type stopfend: str.
    :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend', 'enrichment', and 'expected'. Observed values are always in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, unfiltered fends return value of one. Expected values are returned for 'distance', 'fend', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fend' uses only fend correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
    :type datatype: str.
    :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'compact', 'full', and 'upper'. 'compact' means data are arranged in a N x M x 2 array where N is the number of bins, M is the maximum number of steps between included bin pairs, and data are stored such that bin n,m contains the interaction values between n and n + m + 1. 'full' returns a square, symmetric array of size N x N x 2. 'upper' returns only the flattened upper triangle of a full array, excluding the diagonal of size (N * (N - 1) / 2) x 2.
    :type arraytype: str.
    :param maxdistance: This specifies the maximum coordinate distance between bins that will be included in the array. If set to zero, all distances are included.
    :type maxdistance: str.
    :param skipfiltered: If 'True', all interaction bins for filtered out fends are removed and a reduced-size array is returned.
    :type skipfiltered: bool.
    :param returnmapping: If 'True', a list containing the data array and two 2d array containing first coordinate included and excluded from each bin, and the first fend included and excluded from each bin for the first and second axis is returned. Otherwise only the data array is returned.
    :type returnmapping: bool.
    :returns: Array in format requested with 'arraytype' containing data requested with 'datatype'.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    # check that all values are acceptable
    if datatype not in ['raw', 'fend', 'distance', 'enrichment', 'expected']:
        if not silent:
            print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
        return None
    elif datatype in ['fend', 'enrichment'] and hic.normalization == 'none':
        if not silent:
            print >> sys.stderr, ("Normalization has not been performed yet on this project. Select either 'raw' or 'distance' for datatype. No data returned\n"),
        return None
    # Determine start, stop, startfend, and stopfend
    chrint1 = hic.chr2int[chrom1.strip('chr')]
    chrint2 = hic.chr2int[chrom2.strip('chr')]
    if not binbounds1 is None:
        start1 = binbounds1[0, 0]
        stop1 = binbounds1[-1, 1]
        startfend1 = _find_fend_from_coord(hic, chrint1, start1)
        stopfend1 = _find_fend_from_coord(hic, chrint1, stop1) + 1
    else:
        if start1 is None and startfend1 is None:
            startfend1 = hic.fends['chr_indices'][chrint1]
            while startfend1 < hic.fends['chr_indices'][chrint1 + 1] and hic.filter[startfend1] == 0:
                startfend1 += 1
            if startfend1 == hic.fends['chr_indices'][chrint1 + 1]:
                if not silent:
                    print >> sys.stderr, ("Insufficient data.\n"),
                return None
            start1 = hic.fends['fends']['mid'][startfend1]
            if binsize > 0:
                start1 = (start1 / binsize) * binsize
        elif start1 is None:
            start1 = hic.fends['fends']['mid'][startfend1]
            if binsize > 0:
                start1 = (start1 / binsize) * binsize
        else:
            startfend1 = _find_fend_from_coord(hic, chrint1, start1)
        if (stop1 is None or stop1 == 0) and stopfend1 is None:
            stopfend1 = hic.fends['chr_indices'][chrint1 + 1]
            while stopfend1 > hic.fends['chr_indices'][chrint1] and hic.filter[stopfend1 - 1] == 0:
                stopfend1 -= 1
            stop1 = hic.fends['fends']['mid'][stopfend1 - 1]
            if binsize > 0:
                stop1 = ((stop1 - 1 - start1) / binsize + 1) * binsize + start1
        elif stop1 is None or stop1 == 0:
            stop1 = hic.fends['fends']['mid'][stopfend1 - 1]
            if binsize > 0:
                stop1 = ((stop1 - 1 - start1) / binsize + 1) * binsize + start1
        else:
            if binsize > 0:
                stop1 = ((stop1 - 1 - start1) / binsize + 1) * binsize + start1
            stopfend1 = _find_fend_from_coord(hic, chrint1, stop1) + 1
    if not binbounds1 is None:
        start2 = binbounds1[0, 0]
        stop2 = binbounds1[-1, 1]
        startfend2 = _find_fend_from_coord(hic, chrint2, start2)
        stopfend2 = _find_fend_from_coord(hic, chrint2, stop2) + 1
    else:
        if start2 is None and startfend2 is None:
            startfend2 = hic.fends['chr_indices'][chrint2]
            while startfend2 < hic.fends['chr_indices'][chrint2 + 1] and hic.filter[startfend2] == 0:
                startfend2 += 1
            if startfend2 == hic.fends['chr_indices'][chrint2 + 1]:
                if not silent:
                    print >> sys.stderr, ("Insufficient data.\n"),
                return None
            start2 = hic.fends['fends']['mid'][startfend2]
            if binsize > 0:
                start2 = (start2 / binsize) * binsize
        elif start2 is None:
            start2 = hic.fends['fends']['mid'][startfend2]
            if binsize > 0:
                start2 = (start2 / binsize) * binsize
        else:
            startfend2 = _find_fend_from_coord(hic, chrint2, start2)
        if (stop2 is None or stop2 == 0) and stopfend2 is None:
            stopfend2 = hic.fends['chr_indices'][chrint2 + 1]
            while stopfend2 > hic.fends['chr_indices'][chrint2] and hic.filter[stopfend2 - 1] == 0:
                stopfend2 -= 1
            stop2 = hic.fends['fends']['mid'][stopfend2 - 1]
            if binsize > 0:
                stop2 = ((stop2 - 1 - start2) / binsize + 1) * binsize + start2
        elif stop2 is None or stop2 == 0:
            stop2 = hic.fends['fends']['mid'][stopfend2 - 1]
            if binsize > 0:
                stop2 = ((stop2 - 1 - start2) / binsize + 1) * binsize + start2
        else:
            if binsize > 0:
                stop2 = ((stop2 - 1 - start2) / binsize + 1) * binsize + start2
            stopfend2 = _find_fend_from_coord(hic, chrint2, stop2) + 1
    if not silent:
        print >> sys.stderr, ("Finding %s array for %s:%i-%i by %s:%i-%i...") % (datatype,  chrom1,
                                                                                 start1, stop1, chrom2, start2,
                                                                                 stop2),
    # If datatype is not 'expected', pull the needed slice of data
    if datatype != 'expected':
        if chrint1 < chrint2:
            start_index = hic.data['trans_indices'][startfend1]
            stop_index = hic.data['trans_indices'][stopfend1]
        else:
            start_index = hic.data['trans_indices'][startfend2]
            stop_index = hic.data['trans_indices'][stopfend2]
        if start_index == stop_index:
            if not silent:
                print >> sys.stderr, ("Insufficient data\n"),
            return None
        if chrint1 < chrint2:
            data_indices = hic.data['trans_indices'][startfend1:(stopfend1 + 1)]
        else:
            data_indices = hic.data['trans_indices'][startfend2:(stopfend2 + 1)]
        data_indices -= data_indices[0]
        data = hic.data['trans_data'][start_index:stop_index, :]
        if chrint1 < chrint2:
            data[:, 0] -= startfend1
            data[:, 1] -= startfend2
        else:
            data[:, 0] -= startfend2
            data[:, 1] -= startfend1
    else:
        data_indices = None
        data = None
    # Determine mapping of valid fends to bins
    mapping1 = numpy.zeros(stopfend1 - startfend1, dtype=numpy.int32) - 1
    mapping2 = numpy.zeros(stopfend2 - startfend2, dtype=numpy.int32) - 1
    valid1 = numpy.where(hic.filter[startfend1:stopfend1] > 0)[0].astype(numpy.int32)
    valid2 = numpy.where(hic.filter[startfend2:stopfend2] > 0)[0].astype(numpy.int32)
    mids1 = hic.fends['fends']['mid'][startfend1:stopfend1]
    mids2 = hic.fends['fends']['mid'][startfend2:stopfend2]
    if binsize == 0 and binbounds1 is None:
        if skipfiltered:
            mapping1[valid1] = numpy.arange(valid1.shape[0])
            num_bins1 = valid1.shape[0]
        else:
            mapping1[valid1] = valid1
            num_bins1 = mapping1.shape[0]
    elif not binbounds1 is None:
        start_indices = numpy.searchsorted(binbounds1[:, 0], mids1[valid1], side='right') - 1
        stop_indices = numpy.searchsorted(binbounds1[:, 1], mids1[valid1], side='right')
        where = numpy.where(start_indices == stop_indices)[0]
        valid1 = valid1[where]
        mapping1[valid1] = start_indices[where]
        num_bins1 = binbounds1.shape[0]
    else:
        mapping1[valid1] = (mids1[valid1] - start1) / binsize
        num_bins1 = (stop1 - start1) / binsize
    if binsize == 0 and binbounds2 is None:
        if skipfiltered:
            mapping2[valid2] = numpy.arange(valid2.shape[0])
            num_bins2 = valid2.shape[0]
        else:
            mapping2[valid2] = valid2
            num_bins2 = mapping2.shape[0]
    elif not binbounds2 is None:
        start_indices = numpy.searchsorted(binbounds2[:, 0], mids2[valid2], side='right') - 1
        stop_indices = numpy.searchsorted(binbounds2[:, 1], mids2[valid2], side='right')
        where = numpy.where(start_indices == stop_indices)[0]
        valid2 = valid2[where]
        mapping2[valid2] = start_indices[where]
        num_bins2 = binbounds2.shape[0]
    else:
        mapping2[valid2] = (mids2[valid2] - start2) / binsize
        num_bins2 = (stop2 - start2) / binsize
    # Find maximum interaction partner for each fend
    if num_bins1 < 1 or num_bins2 < 1:
        if not silent:
            print >> sys.stderr, ("Insufficient data\n"),
        return None
    # If correction is required, determine what type and get appropriate data
    if hic.normalization != 'binning' and datatype != 'raw':
        corrections1 = hic.corrections[startfend1:stopfend1]
        corrections2 = hic.corrections[startfend2:stopfend2]
    elif datatype == 'raw':
        corrections1 = numpy.ones(stopfend1 - startfend1, dtype=numpy.float32)
        corrections2 = numpy.ones(stopfend2 - startfend2, dtype=numpy.float32)
    else:
        corrections1 = None
        corrections2 = None
    if ((hic.normalization in ['express', 'probability'] and
            datatype == 'fend') or datatype == 'raw'):
        correction_sums1 = numpy.zeros(num_bins1, dtype=numpy.float64)
        correction_sums2 = numpy.zeros(num_bins2, dtype=numpy.float64)
        if datatype == 'fend':
            correction_sums1[:] = numpy.bincount(mapping1[valid1], weights=corrections1[valid1], minlength=num_bins1)
            correction_sums2[:] = numpy.bincount(mapping2[valid2], weights=corrections2[valid2], minlength=num_bins2)
        else:
            correction_sums1[:] = numpy.bincount(mapping1[valid1], minlength=num_bins1)
            correction_sums2[:] = numpy.bincount(mapping2[valid2], minlength=num_bins2)
    else:
        correction_sums1 = None
        correction_sums2 = None
    if (hic.normalization in ['binning', 'binning-express', 'binning-probability'] and
            datatype not in ['raw', 'distance']):
        binning_corrections = hic.binning_corrections
        binning_num_bins = hic.binning_num_bins
        fend_indices = hic.binning_fend_indices
    else:
        binning_corrections = None
        binning_num_bins = None
        fend_indices = None
    if datatype in ['distance', 'enrichment', 'expected']:
        if 'trans_means' not in hic.__dict__.keys():
            hic.find_trans_means()
        if chrint1 < chrint2:
            index = chrint1 * (hic.fends['chromosomes'].shape[0] - 1) - chrint1 * (chrint1 + 1) / 2 - 1 + chrint2
        else:
            index = chrint2 * (hic.fends['chromosomes'].shape[0] - 1) - chrint2 * (chrint2 + 1) / 2 - 1 + chrint1
        trans_mean = hic.trans_means[index]
    else:
        trans_mean = 1.0
    # Create data array
    if chrint1 < chrint2:
        data_array = numpy.zeros((num_bins1, num_bins2, 2), dtype=numpy.float32)
    else:
        data_array = numpy.zeros((num_bins2, num_bins1, 2), dtype=numpy.float32)
    # Fill in data values
    if chrint1 < chrint2:
        _hic_binning.find_trans_expected(mapping1, mapping2, corrections1, corrections2, binning_corrections,
                                         binning_num_bins, fend_indices, data_array,
                                         correction_sums1, correction_sums2, trans_mean, startfend1, startfend2)
        if datatype != 'expected':
            _hic_binning.find_trans_observed(data, data_indices, mapping1, mapping2, data_array)
        else:
            data_array[:, :, 0] = data_array[:, :, 1]
            data_array[:, :, 1].fill(0)
            corrections1.fill(1.0)
            corrections2.fill(1.0)
            correction_sums1 = numpy.bincount(mapping1[valid1], minlength=num_bins1).astype(numpy.float64)
            correction_sums2 = numpy.bincount(mapping2[valid2], minlength=num_bins2).astype(numpy.float64)
            _hic_binning.find_trans_expected(mapping1, mapping2, corrections1, corrections2, None, None, None,
                                             data_array, correction_sums1, correction_sums2, 1.0, startfend1,
                                             startfend2)
            temp = data_array[:, :, 0]
            data_array[:, :, 0] = data_array[:, :, 1]
            data_array[:, :, 1] = temp
    else:
        _hic_binning.find_trans_expected(mapping2, mapping1, corrections2, corrections1, binning_corrections,
                                         binning_num_bins, fend_indices, data_array,
                                         correction_sums2, correction_sums1, trans_mean, startfend2, startfend1)
        if datatype != 'expected':
            _hic_binning.find_trans_observed(data, data_indices, mapping2, mapping1, data_array)
        else:
            data_array[:, :, 0] = data_array[:, :, 1]
            data_array[:, :, 1].fill(0)
            corrections1.fill(1.0)
            corrections2.fill(1.0)
            correction_sums1 = numpy.bincount(mapping1[valid1], minlength=num_bins1).astype(numpy.float64)
            correction_sums2 = numpy.bincount(mapping2[valid2], minlength=num_bins2).astype(numpy.float64)
            _hic_binning.find_trans_expected(mapping2, mapping1, corrections2, corrections1, None, None, None,
                                             data_array, correction_sums2, correction_sums1, 1.0, startfend2,
                                             startfend1)
            temp = data_array[:, :, 0]
            data_array[:, :, 0] = data_array[:, :, 1]
            data_array[:, :, 1] = temp
    if chrint2 < chrint1:
        data_array = numpy.transpose(data_array, (1, 0, 2))
    if returnmapping:
        bin_mapping1 = numpy.zeros((num_bins1, 4), dtype=numpy.int32)
        if binsize == 0 and binbounds1 is None:
            if skipfiltered:
                bin_mapping1[:, 2] = valid1 + startfend1
            else:
                bin_mapping1[:, 2] = numpy.arange(startfend1, stopfend1)
            bin_mapping1[:, 3] = bin_mapping1[:, 2] + 1
            bin_mapping1[:, 0] = hic.fends['fends']['start'][bin_mapping1[:, 2]]
            bin_mapping1[:, 1] = hic.fends['fends']['stop'][bin_mapping1[:, 2]]
        else:
            if binbounds1 is None:
                bin_mapping1[:, 0] = start1 + binsize * numpy.arange(num_bins1)
                bin_mapping1[:, 1] = bin_mapping1[:, 0] + binsize
            else:
                bin_mapping1[:, :2] = binbounds1
            bin_mapping1[:, 2] = numpy.searchsorted(mids1, bin_mapping1[:, 0]) + startfend1
            bin_mapping1[:, 3] = numpy.searchsorted(mids1, bin_mapping1[:, 1]) + startfend1
        bin_mapping2 = numpy.zeros((num_bins2, 4), dtype=numpy.int32)
        if binsize == 0 and binbounds2 is None:
            if skipfiltered:
                bin_mapping2[:, 2] = valid2 + startfend2
            else:
                bin_mapping2[:, 2] = numpy.arange(startfend2, stopfend2)
            bin_mapping2[:, 3] = bin_mapping2[:, 2] + 1
            bin_mapping2[:, 0] = hic.fends['fends']['start'][bin_mapping2[:, 2]]
            bin_mapping2[:, 1] = hic.fends['fends']['stop'][bin_mapping2[:, 2]]
        else:
            if binbounds2 is None:
                bin_mapping2[:, 0] = start2 + binsize * numpy.arange(num_bins2)
                bin_mapping2[:, 1] = bin_mapping2[:, 0] + binsize
            else:
                bin_mapping2[:, :2] = binbounds2
            bin_mapping2[:, 2] = numpy.searchsorted(mids2, bin_mapping2[:, 0]) + startfend2
            bin_mapping2[:, 3] = numpy.searchsorted(mids2, bin_mapping2[:, 1]) + startfend2
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return [data_array, bin_mapping1, bin_mapping2]
    else:
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return data_array

def bin_trans_array(data_array, data_mapping1, data_mapping2, binsize=10000, binbounds1=None, start1=None, stop1=None,
                    binbounds2=None, start2=None, stop2=None, returnmapping=False, **kwargs):
    """
    Create an array of format 'arraytype' and fill 'binsize' bins or bins defined by 'binbounds' with data provided in the array passed by 'unbinned'.

    :param hic: A :class:`HiC <hifive.hic.HiC>` class object containing fend and count data.
    :type hic: :class:`HiC <hifive.hic.HiC>`
    :param data_array: A 3d array containing data to be binned.
    :type data_array: numpy array
    :param data_mapping1: An N x 4 2d integer array containing the start and stop coordinates, and start and stop fends for each of the N bin ranges along the first axis in 'data_array'.
    :type data_mapping1: numpy array
    :param data_mapping2: An N x 4 2d integer array containing the start and stop coordinates, and start and stop fends for each of the N bin ranges along the second axis in 'data_array'.
    :type data_mapping2: numpy array
    :param binsize: This is the coordinate width of each bin. If binbounds is not None, this value is ignored.
    :type binsize: int.
    :param binbounds1: An array containing start and stop coordinates for a set of user-defined bins along the first axis. Any bin from 'data_array' not falling in a bin is ignored.
    :type binbounds1: numpy array
    :param start1: The coordinate at the beginning of the first bin for the first axis of the binned data. If unspecified, 'start1' will be the first multiple of 'binsize' below the first coordinate from 'data_mapping1'. If 'binbounds1' is given, 'start1' is ignored. Optional.
    :type start1: int.
    :param stop1: The coordinate at the end of the last bin for the first axis of the binned data. If unspecified, 'stop1' will be the first multiple of 'binsize' after the last coordinate from 'data_mapping1'. If needed, 'stop1' is adjusted upward to create a complete last bin. If 'binbounds1' is given, 'stop1' is ignored. Optional.
    :type stop1: int.
    :param binbounds2: An array containing start and stop coordinates for a set of user-defined bins along the second axis. Any bin from 'data_array' not falling in a bin is ignored.
    :type binbounds2: numpy array
    :param start2: The coordinate at the beginning of the first bin for the second axis of the binned data. If unspecified, 'start2' will be the first multiple of 'binsize' below the first coordinate from 'data_mapping2'. If 'binbounds2' is given, 'start2' is ignored. Optional.
    :type start2: int.
    :param stop2: The coordinate at the end of the last bin for the second axis of the binned data. If unspecified, 'stop2' will be the first multiple of 'binsize' after the last coordinate from 'data_mapping2'. If needed, 'stop2' is adjusted upward to create a complete last bin. If 'binbounds2' is given, 'stop2' is ignored. Optional.
    :type stop2: int.
    :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend', 'enrichment', and 'expected'. Observed values are always in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, unfiltered fends return value of one. Expected values are returned for 'distance', 'fend', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fend' uses only fend correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
    :type datatype: str.
    :param returnmapping: If 'True', a list containing the data array and a 2d array containing first coordinate included and excluded from each bin, and the first fend included and excluded from each bin is returned. Otherwise only the data array is returned.
    :type returnmapping: bool.
    :returns: Array in format requested with 'arraytype' containing binned data requested with 'datatype' pulled from 'unbinned'.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    # Determine start and stop, if necessary
    if binbounds1 is None:
        if start1 is None:
            start1 = (data_mapping1[0, 0] / binsize) * binsize
        if stop1 is None:
            stop1 = ((data_mapping1[-1, 1] - 1) / binsize + 1) * binsize
        else:
            stop1 = ((stop1 - 1 - start1) / binsize + 1) * binsize + start1
        num_bins1 = (stop1 - start1) / binsize
        binbounds1 = numpy.zeros((num_bins1, 2), dtype=numpy.int32)
        binbounds1[:, 0] = numpy.arange(num_bins1) * binsize + start1
        binbounds1[:, 1] = binbounds1[:, 0] + binsize
    else:
        num_bins1 = binbounds1.shape[0]
        start1 = binbounds1[0, 0]
        stop1 = binbounds1[0, 1]
    if binbounds2 is None:
        if start2 is None:
            start2 = (data_mapping2[0, 0] / binsize) * binsize
        if stop2 is None:
            stop2 = ((data_mapping2[-1, 1] - 1) / binsize + 1) * binsize
        else:
            stop2 = ((stop2 - 1 - start2) / binsize + 1) * binsize + start2
        num_bins2 = (stop2 - start2) / binsize
        binbounds2 = numpy.zeros((num_bins2, 2), dtype=numpy.int32)
        binbounds2[:, 0] = numpy.arange(num_bins2) * binsize + start2
        binbounds2[:, 1] = binbounds2[:, 0] + binsize
    else:
        num_bins2 = binbounds2.shape[0]
        start2 = binbounds2[0, 0]
        stop2 = binbounds2[0, 1]
    mids1 = (data_mapping1[:, 0] + data_mapping1[:, 1]) / 2
    mids2 = (data_mapping2[:, 0] + data_mapping2[:, 1]) / 2
    if not silent:
        print >> sys.stderr, ("Finding binned trans array..."),
    # Find bin mapping for each fend
    mapping1 = numpy.zeros(mids1.shape[0], dtype=numpy.int32) - 1
    fend_ranges1 = numpy.zeros((binbounds1.shape[0], 2), dtype=numpy.int32)
    for i in range(binbounds1.shape[0]):
        firstbin = numpy.searchsorted(mids1, binbounds1[i, 0])
        lastbin = numpy.searchsorted(mids1, binbounds1[i, 1])
        mapping1[firstbin:lastbin] = i
        fend_ranges1[i, 0] = data_mapping1[firstbin, 2]
        fend_ranges1[i, 1] = data_mapping1[lastbin, 3]
    valid1 = numpy.where(mapping1 >= 0)[0]
    mapping2 = numpy.zeros(mids2.shape[0], dtype=numpy.int32) - 1
    fend_ranges2 = numpy.zeros((binbounds2.shape[0], 2), dtype=numpy.int32)
    for i in range(binbounds2.shape[0]):
        firstbin = numpy.searchsorted(mids2, binbounds2[i, 0])
        lastbin = numpy.searchsorted(mids2, binbounds2[i, 1])
        mapping2[firstbin:lastbin] = i
        fend_ranges2[i, 0] = data_mapping2[firstbin, 2]
        fend_ranges2[i, 1] = data_mapping2[lastbin, 3]
    valid2 = numpy.where(mapping2 >= 0)[0]
    # Create requested array
    binned_array = numpy.zeros((num_bins1, num_bins2, 2), dtype=numpy.float32)
    # Fill in binned data values
    for i in range(valid1.shape[0]):
        binned_array[i, :, 0] = numpy.bincount(mapping2[valid2], weights=data_array[valid1[i], valid2, 0],
                                            minlength=num_bins2)
        binned_array[i, :, 1] = numpy.bincount(mapping2[valid2], weights=data_array[valid1[i], valid2, 1],
                                            minlength=num_bins2)
    # If mapping requested, calculate bin bounds
    if returnmapping:
        mapping1 = numpy.zeros((num_bins1, 4), dtype=numpy.int32)
        mapping1[:, 0] = binbounds1[:, 0]
        mapping1[:, 1] = binbounds1[:, 1]
        mapping1[:, 2:4] = fend_ranges1
        mapping2 = numpy.zeros((num_bins2, 4), dtype=numpy.int32)
        mapping2[:, 0] = binbounds2[:, 0]
        mapping2[:, 1] = binbounds2[:, 1]
        mapping2[:, 2:4] = fend_ranges2
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return [binned_array, mapping1, mapping2]
    else:
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return binned_array

def dynamically_bin_trans_array(unbinned, unbinnedpositions1, unbinnedpositions2, binned, binbounds1, binbounds2,
                                minobservations=10, searchdistance=0, removefailed=False, **kwargs):
    """
    Expand bins in 'binned' to include additional data provided in 'unbinned' as necessary to meet 'minobservations', or 'searchdistance' criteria.

    :param unbinned: A 3d array containing data to be used for filling expanding bins. This array should be  N x M x 2, where N is the number of bins or fends from the first chromosome and M is the number of bins or fends from the second chromosome.
    :type unbinned: numpy array
    :param unbinnedpositions1: A 2d integer array indicating the first and last coordinate of each bin along the first axis in 'unbinned' array.
    :type unbinnedpositions1: numpy array
    :param unbinnedpositions2: A 2d integer array indicating the first and last coordinate of each bin along the first axis in 'unbinned' array.
    :type unbinnedpositions2: numpy array
    :param binned: A 3d array containing binned data to be dynamically binned. This array should be  N x M x 2, where N is the number of bins from the first chromosome and M is the number of bins from the second chromosome. Data in this array will be altered by this function.
    :type binned: numpy array
    :param binbounds1: An integer array indicating the start and end position of each bin from the first chromosome in the 'binned' array. This array should be N x 2, where N is the size of the first dimension of 'binned'.
    :type binbounds1: numpy array
    :param binbounds2: An integer array indicating the start and end position of each bin from the second chromosome in the 'binned' array. This array should be N x 2, where N is the size of the second dimension of 'binned'.
    :type binbounds2: numpy array
    :param minobservations: The fewest number of observed reads needed for a bin to counted as valid and stop expanding.
    :type minobservations: int.
    :param searchdistance: The furthest distance from the bin minpoint to expand bounds. If this is set to zero, there is no limit on expansion distance.
    :type searchdistance: int.
    :param removefailed: If a non-zero 'searchdistance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'removefailed' is True, the observed and expected values for that bin are zero.
    :type removefailed: bool.
    :returns: None
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    if not silent:
        print >> sys.stderr, ("Dynamically binning data..."),
    # Determine bin edges relative to unbinned positions
    unbinnedmids1 = (unbinnedpositions1[:, 0] + unbinnedpositions1[:, 1]) / 2
    unbinnedmids2 = (unbinnedpositions2[:, 0] + unbinnedpositions2[:, 1]) / 2
    binedges1 = numpy.zeros(binbounds1.shape, dtype=numpy.int32)
    binedges1[:, 0] = numpy.searchsorted(unbinnedmids1, binbounds1[:, 0])
    binedges1[:, 1] = numpy.searchsorted(unbinnedmids1, binbounds1[:, 1])
    binedges2 = numpy.zeros(binbounds2.shape, dtype=numpy.int32)
    binedges2[:, 0] = numpy.searchsorted(unbinnedmids2, binbounds2[:, 0])
    binedges2[:, 1] = numpy.searchsorted(unbinnedmids2, binbounds2[:, 1])
    # Determine bin midpoints
    mids1 = (binbounds1[:, 0] + binbounds1[:, 1]) / 2
    mids2 = (binbounds2[:, 0] + binbounds2[:, 1]) / 2
    # Dynamically bin using appropriate array type combination
    _hic_binning.dynamically_bin_trans(unbinned, unbinnedmids1, unbinnedmids2, binned, binedges1,
                                   binedges2, mids1, mids2, minobservations, searchdistance, int(removefailed))
    if not silent:
        print >> sys.stderr, ("Done\n"),
    return None

def write_heatmap_dict(hic, filename, binsize, includetrans=True, datatype='enrichment', chroms=[], 
                       dynamically_binned=False, minobservations=0, searchdistance=0, expansion_binsize=0,
                       removefailed=False, **kwargs):
    """
    Create an h5dict file containing binned interaction arrays, bin positions, and an index of included chromosomes. This function is MPI compatible.

    :param hic: A :class:`HiC <hifive.hic.HiC>` class object containing fend and count data.
    :type hic: :class:`HiC <hifive.hic.HiC>`
    :param filename: Location to write h5dict object to.
    :type filename: str.
    :param binsize: Size of bins for interaction arrays.
    :type binsize: int.
    :param includetrans: Indicates whether trans interaction arrays should be calculated and saved.
    :type includetrans: bool.
    :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend', 'enrichment', and 'expected'. Observed values are always in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, unfiltered fends return value of one. Expected values are returned for 'distance', 'fend', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fend' uses only fend correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
    :type datatype: str.
    :param chroms: A list of chromosome names indicating which chromosomes should be included. If left empty, all chromosomes are included. Optional.
    :type chroms: list
    :param dynamically_binned: If 'True', return dynamically binned data.
    :type dynamically_binned: bool.
    :param minobservations: The fewest number of observed reads needed for a bin to counted as valid and stop expanding.
    :type minobservations: int.
    :param searchdistance: The furthest distance from the bin minpoint to expand bounds. If this is set to zero, there is no limit on expansion distance.
    :type searchdistance: int.
    :param expansion_binsize: The size of bins to use for data to pull from when expanding dynamic bins. If set to zero, unbinned data is used.
    :type expansion_binsize: int.
    :param removefailed: If a non-zero 'searchdistance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'removefailed' is True, the observed and expected values for that bin are zero.
    :type removefailed: bool.
    :returns: None
    """
    # check if MPI is available
    if 'mpi4py' in sys.modules.keys():
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        num_procs = comm.Get_size()
    else:
        comm = None
        rank = 0
        num_procs = 1
    if ('silent' in kwargs and kwargs['silent']) or rank > 0:
        silent = True
    else:
        silent = False
    # Check if trans mean is needed and calculate if not already done
    if includetrans and datatype in ['distance', 'enrichment'] and 'trans_mean' not in hic.__dict__.keys():
        hic.find_trans_means()
    # Check if filename already exists, and remove if it does
    if rank == 0:
        if os.path.exists(filename):
            if not silent:
                print >> sys.stderr, ("%s already exists, overwriting.") % filename
            subprocess.call('rm %s' % filename, shell=True)
        if not silent:
            print >> sys.stderr, ("Creating binned heatmap...\n"),
        output = h5py.File(filename, 'w')
        output.attrs['resolution'] = binsize
        # If chromosomes not specified, fill list
        if len(chroms) == 0:
            chroms = list(hic.fends['chromosomes'][...])
        # Assemble list of requested arrays
        needed = []
        chr_indices = hic.fends['chr_indices'][...]
        for i in range(len(chroms))[::-1]:
            chrom = chroms[i]
            chrint = hic.chr2int[chrom]
            if numpy.sum(hic.filter[chr_indices[chrint]:chr_indices[chrint + 1]]) > 0:
                needed.append((chrom,))
            else:
                del chroms[i]
        if includetrans:
            for i in range(len(chroms)-1):
                for j in range(i + 1, len(chroms)):
                    needed.append((chroms[i],chroms[j]))
        if num_procs == 1:
            node_needed = needed
        else:
            node_ranges = numpy.round(numpy.linspace(0, len(needed), num_procs + 1)).astype(numpy.int32)
            for i in range(1, num_procs):
                comm.send(needed[node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
            node_needed = needed[node_ranges[0]:node_ranges[1]]
    else:
        node_needed = comm.recv(source=0, tag=11)
    heatmaps = {}
    # Find heatmaps
    for chrom in node_needed:
        if len(chrom) == 1:
            # Find cis heatmap
            # determine if data is to be dynamically binned
            if not dynamically_binned:
                heatmaps[chrom] = find_cis_signal(hic, chrom[0], binsize=binsize, datatype=datatype,
                                                  arraytype='upper', returnmapping=True, silent=silent,
                                                  skipfiltered=True)
            else:
                temp = find_cis_signal(hic, chrom[0], binsize=expansion_binsize, datatype=datatype, arraytype='upper',
                                       returnmapping=True, silent=silent)
                if temp is None:
                    continue
                expansion, exp_mapping = temp
                binned, mapping = find_cis_signal(hic, chrom[0], binsize=binsize, datatype=datatype,
                                                  arraytype='upper', returnmapping=True, silent=silent)
                dynamically_bin_cis_array(expansion, exp_mapping, binned, mapping, minobservations=minobservations,
                                          searchdistance=searchdistance, removefailed=removefailed, silent=silent)

                heatmaps[chrom] = [binned, mapping]
        else:
            # Find trans heatmap
            # determine if data is to be dynamically binned
            if not dynamically_binned:
                heatmaps[chrom] = find_trans_signal(hic, chrom[0], chrom[1],  binsize=binsize, datatype=datatype,
                                                    returnmapping=False, silent=silent, skipfiltered=True)
            else:
                temp = find_trans_signal(hic, chrom[0], chrom[1], binsize=expansion_binsize, datatype=datatype, 
                                         returnmapping=True, silent=silent)
                if temp is None:
                    continue
                expansion, exp_mapping1, exp_mapping2 = temp
                binned, mapping1, mapping2 = find_trans_signal(hic, chrom[0], chrom[1], binsize=binsize,
                                                               datatype=datatype, returnmapping=True, silent=silent)
                dynamically_bin_trans_array(expansion, exp_mapping1, exp_mapping2, binned, mapping1, mapping2,
                                            minobservations=minobservations, searchdistance=searchdistance,
                                            removefailed=removefailed, silent=silent)

                heatmaps[chrom] = binned
        # Check if array contains data
        if heatmaps[chrom] is None or heatmaps[chrom][0].shape[0] == 0:
            del heatmaps[chrom]
    # Collect heatmaps at node 0 and write to h5dict
    if rank == 0:
        if num_procs > 1:
            for i in range(1, num_procs):
                if node_ranges[i + 1] - node_ranges[i] > 0:
                    temp = comm.recv(source=i, tag=11)
                    heatmaps.update(temp)
            del temp
        for chrom in heatmaps.keys():
            if len(chrom) == 1:
                output.create_dataset('%s.counts' % chrom[0], data=heatmaps[chrom][0][:, 0])
                output.create_dataset('%s.expected' % chrom[0], data=heatmaps[chrom][0][:, 1])
                output.create_dataset('%s.positions' % chrom[0], data=heatmaps[chrom][1][:, :2])
            else:
                output.create_dataset('%s_by_%s.counts' % (chrom[0], chrom[1]), data=heatmaps[chrom][:, :, 0])
                output.create_dataset('%s_by_%s.expected' % (chrom[0], chrom[1]), data=heatmaps[chrom][:, :, 1])
        output.create_dataset('chromosomes', data=numpy.array(chroms))
        if 'history' in kwargs:
            output.attrs['history'] = kwargs['history']
        output.close()
        if not silent:
            print >> sys.stderr, ("Creating binned heatmap...Done\n"),
    else:
        if len(heatmaps) > 0:
            comm.send(heatmaps, dest=0, tag=11)
        del heatmaps
    return None


def find_multiresolution_heatmap(hic, chrom, start, stop, chrom2=None, start2=None, stop2=None, minbinsize=5000,
                                 maxbinsize=12800000, minobservations=5, datatype='fend', midbinsize=40000,
                                 silent=True):
    """
    Create a multi-resolution data and index heatmap array for a chromosome or chromosome pair.

    :param hic: A :class:`HiC <hifive.hic.HiC>` class object containing fend and count data.
    :type hic: :class:`HiC <hifive.hic.HiC>`
    :param chrom: The first (or only) chromosome to find the multi-resolution heatmap for.
    :type chrom: str.
    :param start: The first bin start coordinate.
    :type start: int.
    :param stop: The last bin stop coordinate. The difference between start and stop must be a multiple of maxbinsize.
    :type stop: int.
    :param chrom2: The second chromosome to find the multi-resolution heatmap for. If None, an intra-chromosomal multi-resolution heatmap is returned for chrom.
    :type chrom2: str.
    :param start2: The first bin start coordinate for the second chromosome.
    :type start2: int.
    :param stop2: The last bin stop coordinate for the second chromosome. The difference between start and stop must be a multiple of maxbinsize.
    :type stop2: int.
    :param maxbinsize: The maximum sized bin (lowest resolution) heatmap to be produced for each chromosome.
    :type maxbinsize: int.
    :param minbinsize: The minimum sized bin (highest resolution) heatmap to be produced for each chromosome.
    :type minbinsize: int.
    :param minobservations: The minimum number of reads needed for a bin to be considered valid and be included in the heatmap.
    :type minobservations: int.
    :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend', and 'enrichment'. Observed values are always in the first index along the last axis. If 'raw' is specified, unfiltered fends return value of one. Expected values are returned for 'distance', 'fend', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fend' uses only fend correction values, and 'enrichment' uses both correction and distance mean values.
    :type datatype: str.
    :param midbinsize: This is used to determine the smallest bin size (highest resolution) complete heatmap to generate in producing the multi-resolution heatmap. It does not affect the resulting output but can be used to limit the total memory usage, with higher values using less memory but more time.
    :type midbinsize: int.
    :param silent: Indicates whether to display messages or not.
    :type silent: bool.
    """
    # check that all values are acceptable
    if datatype not in ['raw', 'fend', 'distance', 'enrichment']:
        if not silent:
            print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
        return None
    if not chrom2 is None and (start2 is None or stop2 is None):
        if not silent:
            print >> sys.stderr, ("Need values for start2 and stop2. No data returned\n"),
        return None
    if (stop - start) % maxbinsize != 0 or (not chrom2 is None and (stop2 - start2) % maxbinsize != 0):
        if not silent:
            print >> sys.stderr, ("Genomic intervals must be multiples of maxbinsize. No data returned\n"),
        return None
    res_levels = numpy.round(numpy.log(maxbinsize / minbinsize) / numpy.log(2.0)).astype(numpy.int32)
    if maxbinsize != minbinsize * 2 ** res_levels:
        if not silent:
            print >> sys.stderr, ("Maxbinsize must be a multiple of 2^N and minbinsize for an integer N. No data returned\n"),
        return None
    if not silent:
        if chrom2 is None:
            target = chrom
        else:
            target = '%s by %s' % (chrom, chrom2)
        print >> sys.stderr, ("\r%s\rFinding multi-resolution heatmap for %s...") % (' ' * 80, target),
    # determine if finding cis or trans multi-resolution heatmap
    chrint = hic.chr2int[chrom]
    chrint2 = None
    span = stop - start
    startfend = _find_fend_from_coord(hic, chrint, start)
    stopfend = _find_fend_from_coord(hic, chrint, stop)
    if chrom2 is None:
        trans = False
    else:
        span2 = stop2 - start2
        chrint2 = hic.chr2int[chrom2]
        trans = True
        startfend2 = _find_fend_from_coord(hic, chrint2, start2)
        stopfend2 = _find_fend_from_coord(hic, chrint2, stop2)
    # determine actual midresolution limit
    temp = maxbinsize
    while temp / 2 >= max(midbinsize, minbinsize):
        temp /= 2
    midbinsize = temp
    # pull relevant data
    n = span / midbinsize
    valid = numpy.where(hic.filter[startfend:stopfend])[0].astype(numpy.int32)
    fend_nums = valid + startfend
    mids = hic.fends['fends']['mid'][fend_nums] - start
    binbounds = numpy.round(numpy.linspace(0, span, n + 1)).astype(numpy.int32)
    bin_mids = (binbounds[:-1] + binbounds[1:]) / 2
    mapping = numpy.empty(stopfend - startfend, dtype=numpy.int32)
    mapping.fill(-1)
    mapping[valid] = numpy.arange(valid.shape[0])
    binmapping = mids / midbinsize
    obs_indices = numpy.searchsorted(mids, binbounds).astype(numpy.int32)
    if hic.normalization in ['express', 'probability', 'binning-express', 'binning-probability']:
        corrections = hic.corrections[fend_nums]
        correction_sums = numpy.bincount(binmapping, weights=corrections, minlength=n).astype(numpy.float32)
    else:
        corrections = None
        correction_sums = None
    if hic.normalization in ['binning', 'binning-express', 'binning-probability']:
        binning_corrections = hic.binning_corrections
        fend_indices = hic.binning_fend_indices[fend_nums, :, :]
    else:
        binning_corrections = None
        fend_indices = None
    if datatype in ['distance', 'enrichment']:
        distance_parameters = hic.distance_parameters
        chrom_mean = hic.chromosome_means[chrint]
    else:
        distance_parameters = None
        chrom_mean = 0.0
    if trans:
        m = span2 / midbinsize
        valid2 = numpy.where(hic.filter[startfend2:stopfend2])[0]
        fend_nums2 = valid2 + startfend2
        mids2 = hic.fends['fends']['mid'][fend_nums2] - start2
        binbounds2 = numpy.round(numpy.linspace(0, span2, m + 1)).astype(numpy.int32)
        bin_mids2 = (binbounds2[:-1] + binbounds2[1:]) / 2
        obs_indices2 = numpy.searchsorted(mids2, binbounds2).astype(numpy.int32)
        mapping2 = numpy.empty(stopfend2 - startfend2, dtype=numpy.int32)
        mapping2.fill(-1)
        mapping2[valid2] = numpy.arange(valid2.shape[0])
        binmapping2 = mids2 / midbinsize
        if hic.normalization in ['express', 'probability', 'binning-express', 'binning-probability']:
            corrections2 = hic.corrections[fend_nums2]
            correction_sums2 = numpy.bincount(binmapping2, weights=corrections2, minlength=m).astype(numpy.float32)
        else:
            corrections2 = None
            correction_sums2 = None
        if hic.normalization in ['binning', 'binning-express', 'binning-probability']:
            fend_indices2 = hic.binning_fend_indices[fend_nums2, :, :]
        else:
            fend_indices2 = None
        if datatype in ['distance', 'enrichment']:
            if 'trans_means' not in hic.__dict__.keys():
                hic.find_trans_means()
            if chrint < chrint2:
                index = chrint * (hic.fends['chromosomes'].shape[0] - 1) - chrint * (chrint + 1) / 2 - 1 + chrint2
            else:
                index = chrint2 * (hic.fends['chromosomes'].shape[0] - 1) - chrint2 * (chrint2 + 1) / 2 - 1 + chrint
            chrom_mean = hic.trans_means[index]
        # pull relevant trans observations and remap
        if chrint2 < chrint:
            start_index = hic.data['trans_indices'][startfend2]
            stop_index = hic.data['trans_indices'][stopfend2]
            data = hic.data['trans_data'][start_index:stop_index, :]
            data_indices = hic.data['trans_indices'][startfend2:(stopfend2 + 1)]
            data_indices -= data_indices[0]
            num_data = _hic_binning.remap_mrh_data(
                        data,
                        data_indices,
                        mapping2,
                        mapping,
                        startfend,
                        stopfend,
                        startfend2,
                        stopfend2 - startfend2,
                        1)
        else:
            start_index = hic.data['trans_indices'][startfend]
            stop_index = hic.data['trans_indices'][stopfend]
            data = hic.data['trans_data'][start_index:stop_index, :]
            data_indices = hic.data['trans_indices'][startfend:(stopfend + 1)]
            data_indices -= data_indices[0]
            num_data = _hic_binning.remap_mrh_data(
                        data,
                        data_indices,
                        mapping,
                        mapping2,
                        startfend2,
                        stopfend2,
                        startfend,
                        stopfend - startfend,
                        0)
    else:
        # pull relevant cis observations
        start_index = hic.data['cis_indices'][startfend]
        stop_index = hic.data['cis_indices'][stopfend]
        data = hic.data['cis_data'][start_index:stop_index, :]
        data_indices = hic.data['cis_indices'][startfend:(stopfend + 1)]
        data_indices -= data_indices[0]
        num_data = _hic_binning.remap_mrh_data(
                    data,
                    data_indices,
                    mapping,
                    None,
                    startfend,
                    stopfend,
                    startfend,
                    stopfend - startfend,
                    0)
    if trans and chrint2 < chrint:
        data = data[numpy.lexsort((data[:, 1], data[:, 0])), :]
    data_indices = numpy.r_[0, numpy.bincount(data[:num_data, 0], minlength=valid.shape[0])].astype(numpy.int64)
    for i in range(1, data_indices.shape[0]):
        data_indices[i] += data_indices[i - 1]
    data = data[:data_indices[-1], 1:]
    # convert observations into binned matrix
    if trans:
        observed = numpy.zeros((n, m), dtype=numpy.int32)
    else:
        observed = numpy.zeros((n, n), dtype=numpy.int32)
        binmapping2 = None
    _hic_binning.find_mrh_observed(
        data,
        data_indices,
        observed,
        binmapping,
        binmapping2)
    expected = numpy.zeros(observed.shape, dtype=numpy.float32)
    datatype_int = {'raw':0, 'fend':1, 'distance':2, 'enrichment':3}
    dt_int = datatype_int[datatype]
    if trans:
        _hic_binning.find_mrh_trans_expected(
            expected,
            binmapping,
            binmapping2,
            obs_indices,
            obs_indices2,
            corrections,
            corrections2,
            correction_sums,
            correction_sums2,
            binning_corrections,
            fend_indices,
            fend_indices2,
            chrom_mean,
            dt_int)
    else:
        _hic_binning.find_mrh_cis_expected(
            expected,
            fend_nums,
            binmapping,
            mapping,
            mids,
            obs_indices,
            corrections,
            correction_sums,
            binning_corrections,
            fend_indices,
            distance_parameters,
            chrom_mean,
            dt_int)
    # find features for largest binned data array
    n_bins = span / maxbinsize
    m_bins = 0
    binbounds = numpy.linspace(0, span, n_bins + 1)
    if trans:
        m_bins = span2 / maxbinsize
        binbounds2 = numpy.linspace(0, span2, m_bins + 1)
    # find fend assignments for largest bin sizes
    obs_indices = numpy.searchsorted(bin_mids, binbounds).astype(numpy.int32)
    if trans:
        obs_indices2 = numpy.searchsorted(bin_mids2, binbounds2).astype(numpy.int32)
    else:
        obs_indices2 = None
    # make data arrays to hold output
    if trans:
        current_level_data = numpy.zeros(n_bins * m_bins, dtype=numpy.float32)
    else:
        current_level_data = numpy.zeros((n_bins * (n_bins + 1)) / 2, dtype=numpy.float32)
    current_level_indices = numpy.empty(current_level_data.shape, dtype=numpy.int32)
    current_level_indices.fill(-1)
    current_level_shapes = numpy.zeros(current_level_data.shape, dtype=numpy.int32)
    bin_position = numpy.empty(current_level_data.shape, dtype=numpy.int32)
    bin_position.fill(-1)
    # find largest binned data array
    if trans:
        _hic_binning.make_trans_mrh_toplevel(observed,
                                             expected,
                                             current_level_data,
                                             obs_indices,
                                             obs_indices2,
                                             bin_position,
                                             minobservations)
    else:
        _hic_binning.make_cis_mrh_toplevel(observed,
                                           expected,
                                           current_level_data,
                                           obs_indices,
                                           bin_position,
                                           minobservations)
    all_data = [current_level_data]
    all_indices = [current_level_indices]
    all_shapes = [current_level_shapes]
    # find subpartitioning for all valid bins for each resolution level
    resolution = maxbinsize / 2
    if trans:
        pos = n_bins * m_bins
    else:
        pos = (n_bins * (n_bins + 1)) / 2
    # find levels below the first but above or equal to midbinsize
    while resolution >= midbinsize:
        prev_bin_position = bin_position
        bin_position = numpy.empty(prev_bin_position.shape[0] * 4, dtype=numpy.int32)
        bin_position.fill(-1)
        prev_level_data = all_data[-1]
        current_level_data = numpy.empty(prev_level_data.shape[0] * 4, dtype=numpy.float32)
        current_level_data.fill(numpy.nan)
        prev_level_indices = all_indices[-1]
        prev_level_shapes = all_shapes[-1]
        prev_n_bins = n_bins
        prev_m_bins = 0
        n_bins = span / resolution
        binbounds = numpy.linspace(0, span, n_bins + 1)
        obs_indices = numpy.searchsorted(bin_mids, binbounds).astype(numpy.int32)
        if trans:
            prev_m_bins = m_bins
            m_bins = span2 / resolution
            binbounds2 = numpy.linspace(0, span2, m_bins + 1)
            obs_indices2 = numpy.searchsorted(bin_mids2, binbounds2).astype(numpy.int32)
        if trans:
            _hic_binning.make_trans_mrh_midlevel(observed,
                                                 expected,
                                                 current_level_data,
                                                 prev_level_data,
                                                 prev_level_indices,
                                                 prev_level_shapes,
                                                 obs_indices,
                                                 obs_indices2,
                                                 prev_bin_position,
                                                 bin_position,
                                                 prev_m_bins,
                                                 m_bins,
                                                 minobservations,
                                                 pos)
        else:
            _hic_binning.make_cis_mrh_midlevel(observed,
                                               expected,
                                               current_level_data,
                                               prev_level_data,
                                               prev_level_indices,
                                               prev_level_shapes,
                                               obs_indices,
                                               prev_bin_position,
                                               bin_position,
                                               prev_n_bins,
                                               n_bins,
                                               minobservations,
                                               pos)
        where = numpy.where(bin_position >= 0)[0]
        pos += where.shape[0]
        bin_position = bin_position[where]
        all_data.append(current_level_data[where])
        if resolution > minbinsize:
            all_indices.append(numpy.empty(all_data[-1].shape[0], dtype=numpy.int32))
            all_indices[-1].fill(-1)
            all_shapes.append(numpy.zeros(all_data[-1].shape[0], dtype=numpy.int32))
        resolution /= 2
    # find levels below midbinsize
    if midbinsize > minbinsize:
        while resolution >= minbinsize:
            prev_bin_position = bin_position
            bin_position = numpy.empty(prev_bin_position.shape[0] * 4, dtype=numpy.int32)
            bin_position.fill(-1)
            prev_level_data = all_data[-1]
            current_level_data = numpy.empty(prev_level_data.shape[0] * 4, dtype=numpy.float32)
            current_level_data.fill(numpy.nan)
            prev_level_indices = all_indices[-1]
            prev_level_shapes = all_shapes[-1]
            prev_n_bins = n_bins
            prev_m_bins = 0
            n_bins = span / resolution
            binbounds = numpy.linspace(0, span, n_bins + 1)
            obs_indices = numpy.searchsorted(mids, binbounds).astype(numpy.int32)
            correction_sums = numpy.zeros(n_bins, dtype=numpy.float32)
            for i in range(n_bins):
                correction_sums[i] = numpy.sum(corrections[obs_indices[i]:obs_indices[i + 1]])
            if trans:
                prev_m_bins = m_bins
                m_bins = span2 / resolution
                binbounds2 = numpy.linspace(0, span2, m_bins + 1)
                obs_indices2 = numpy.searchsorted(mids2, binbounds2).astype(numpy.int32)
                correction_sums2 = numpy.zeros(m_bins, dtype=numpy.float32)
                for i in range(m_bins):
                    correction_sums2[i] = numpy.sum(corrections2[obs_indices2[i]:obs_indices2[i + 1]])
            if trans:
                _hic_binning.make_trans_mrh_lowerlevel(data,
                                                       data_indices,
                                                       correction_sums,
                                                       correction_sums2,
                                                       current_level_data,
                                                       prev_level_data,
                                                       prev_level_indices,
                                                       prev_level_shapes,
                                                       obs_indices,
                                                       obs_indices2,
                                                       prev_bin_position,
                                                       bin_position,
                                                       prev_m_bins,
                                                       m_bins,
                                                       minobservations,
                                                       pos)
            else:
                _hic_binning.make_cis_mrh_lowerlevel(data,
                                                     data_indices,
                                                     corrections,
                                                     correction_sums,
                                                     fend_nums,
                                                     current_level_data,
                                                     prev_level_data,
                                                     prev_level_indices,
                                                     prev_level_shapes,
                                                     obs_indices,
                                                     prev_bin_position,
                                                     bin_position,
                                                     prev_n_bins,
                                                     n_bins,
                                                     minobservations,
                                                     pos)
            where = numpy.where(bin_position >= 0)[0]
            pos += where.shape[0]
            bin_position = bin_position[where]
            all_data.append(current_level_data[where])
            if resolution > minbinsize:
                all_indices.append(numpy.empty(all_data[-1].shape[0], dtype=numpy.int32))
                all_indices[-1].fill(-1)
                all_shapes.append(numpy.zeros(all_data[-1].shape[0], dtype=numpy.int32))
            resolution /= 2
    data = all_data[0]
    for i in range(1, len(all_data)):
        where = numpy.where(numpy.logical_not(numpy.isnan(all_data[i])))
        data = numpy.hstack((data, all_data[i]))
        all_data[i] = None
    indices = all_indices[0]
    for i in range(1, len(all_indices)):
        indices = numpy.hstack((indices, all_indices[i]))
        all_indices[i] = None
    shapes = all_shapes[0]
    for i in range(1, len(all_shapes)):
        shapes = numpy.hstack((shapes, all_shapes[i]))
        all_shapes[i] = None
    if not silent:
        print >> sys.stderr, ("Done\n"),
    return [data, indices, shapes]
