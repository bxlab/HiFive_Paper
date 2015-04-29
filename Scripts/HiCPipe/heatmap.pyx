#(c) 2013 Emory University. All Rights Reserved

"""These functions provide increased speed in handling the signal-binning
functions necessary for creating HiCPipe heatmaps.
"""

import cython
cimport numpy as np
import numpy

ctypedef np.float32_t DTYPE_t
ctypedef np.float64_t DTYPE_64_t
ctypedef np.int32_t DTYPE_int_t
ctypedef np.int64_t DTYPE_int64_t
ctypedef np.uint32_t DTYPE_uint_t
ctypedef np.int8_t DTYPE_int8_t
cdef double Inf = numpy.inf

cdef extern from "math.h":
    double exp(double x) nogil
    double log(double x) nogil
    double log10(double x) nogil
    double sqrt(double x) nogil
    double pow(double x, double x) nogil
    double abs(double x) nogil
    double round(double x) nogil
    double floor(double x) nogil
    double ceil(double x) nogil


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bin_cis(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int64_t, ndim=1] indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=2] fends not None,
        np.ndarray[DTYPE_t, ndim=2] map_matrix not None,
        np.ndarray[DTYPE_t, ndim=2] gc_matrix not None,
        np.ndarray[DTYPE_t, ndim=2] len_matrix not None,
        np.ndarray[DTYPE_int_t, ndim=1] observed not None,
        np.ndarray[DTYPE_t, ndim=1] expected not None,
        int start_fend,
        int num_bins,
        int binsize):
    cdef long long int i, j, fend1, fend2, bin1, bin2, index
    cdef double value
    cdef long long int num_fends = indices.shape[0] - 1
    with nogil:
        # find observed bin values
        for i in range(num_fends - 2):
            fend1 = start_fend + i
            if fends[fend1, 0] == 0:
                continue
            bin1 = mids[fend1] / binsize
            index = bin1 * num_bins - bin1 * (bin1 + 1) / 2 - bin1 - 1
            for j in range(indices[i], indices[i + 1]):
                fend2 = data[j, 1]
                if fends[fend2, 0] == 0:
                    continue
                bin2 = mids[fend2] / binsize
                if bin2 <= bin1:
                    continue
                observed[index + bin2] += data[j, 2]
        # find expected bin values
        for i in range(num_fends - 2):
            fend1 = start_fend + i
            if fends[fend1, 0] == 0:
                continue
            bin1 = mids[fend1] / binsize
            index = bin1 * num_bins - bin1 * (bin1 + 1) / 2 - bin1 - 1
            fend2 = fend1 + 2
            bin2 = mids[fend2] / binsize
            if fends[fend2, 0] != 0 and bin2 != bin1:
                value = 1.0
                value *= map_matrix[fends[fend1, 1], fends[fend2, 1]]
                value *= gc_matrix[fends[fend1, 2], fends[fend2, 2]]
                value *= len_matrix[fends[fend1, 3], fends[fend2, 3]]
                expected[index + bin2] += value
            for j in range((i / 2 + 2) * 2, num_fends):
                fend2 = j + start_fend
                if fends[fend2, 0] == 0:
                    continue
                bin2 = mids[fend2] / binsize
                if bin2 <= bin1:
                    continue
                value = 1.0
                value *= map_matrix[fends[fend1, 1], fends[fend2, 1]]
                value *= gc_matrix[fends[fend1, 2], fends[fend2, 2]]
                value *= len_matrix[fends[fend1, 3], fends[fend2, 3]]
                expected[index + bin2] += value
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bin_trans(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int64_t, ndim=1] indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=2] fends not None,
        np.ndarray[DTYPE_t, ndim=2] map_matrix not None,
        np.ndarray[DTYPE_t, ndim=2] gc_matrix not None,
        np.ndarray[DTYPE_t, ndim=2] len_matrix not None,
        np.ndarray[DTYPE_int_t, ndim=2] observed not None,
        np.ndarray[DTYPE_t, ndim=2] expected not None,
        int start_fend1,
        int start_fend2,
        int stop_fend2,
        int binsize):
    cdef long long int i, j, fend1, fend2, bin1, bin2
    cdef double value
    cdef long long int num_fends = indices.shape[0] - 1
    with nogil:
        # find observed bin values
        for i in range(num_fends):
            fend1 = start_fend1 + i
            if fends[fend1, 0] == 0:
                continue
            bin1 = mids[fend1] / binsize
            for j in range(indices[i], indices[i + 1]):
                fend2 = data[j, 1]
                if fend2 < start_fend2 or fend2 >= stop_fend2 or fends[fend2, 0] == 0:
                    continue
                bin2 = mids[fend2] / binsize
                observed[bin1, bin2] += data[j, 2]
        # find expected bin values
        for i in range(num_fends):
            fend1 = start_fend1 + i
            if fends[fend1, 0] == 0:
                continue
            bin1 = mids[fend1] / binsize
            for fend2 in range(start_fend2, stop_fend2):
                if fends[fend2, 0] == 0:
                    continue
                bin2 = mids[fend2] / binsize
                value = 1.0
                value *= map_matrix[fends[fend1, 1], fends[fend2, 1]]
                value *= gc_matrix[fends[fend1, 2], fends[fend2, 2]]
                value *= len_matrix[fends[fend1, 3], fends[fend2, 3]]
                expected[bin1, bin2] += value
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_compact_array(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int64_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=2] fends not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=2] gc_mat not None,
        np.ndarray[DTYPE_t, ndim=2] len_mat not None,
        np.ndarray[DTYPE_t, ndim=2] map_mat not None,
        np.ndarray[DTYPE_int_t, ndim=1] max_fend not None,
        np.ndarray[DTYPE_t, ndim=3] signal not None):
    cdef long long int fend1, fend2, i, j, k, gc_bin, len_bin, map_bin
    cdef double value
    cdef long long int num_fends = mapping.shape[0]
    with nogil:
        # if finding anything but expected, fill in actual signal
        for i in range(num_fends - 1):
            fend1 = mapping[i]
            if fends[fend1, 0] == 0:
                continue
            j = i + 1
            k = indices[fend1]
            while j < max_fend[i]:
                fend2 = mapping[j]
                if fends[fend2, 0] > 0:
                    while k < indices[fend1 + 1] and data[k, 1] < fend2:
                        k += 1
                    if k < indices[fend1 + 1] and data[k, 1] == fend2:
                        signal[i, j - i - 1, 0] = data[k, 2]
                j += 1
        # fill in expected signal
        for i in range(num_fends - 1):
            fend1 = mapping[i]
            if fends[fend1, 0] == 0:
                continue
            map_bin = fends[fend1, 1]
            gc_bin = fends[fend1, 2]
            len_bin = fends[fend1, 3]
            j = i + 1
            fend2 = mapping[j]
            k = 0
            # find opposite strand adjacents, skipping same fragment and same strand adjacents
            while j < num_fends and fend2 < (fend1 / 2 + 2) * 2:
                if fend2 == fend1 + 2 and fends[fend2, 0] == 1:
                    value = 1.0
                    value *= map_mat[map_bin, fends[fend2, 1]]
                    value *= gc_mat[gc_bin, fends[fend2, 2]]
                    value *= len_mat[len_bin, fends[fend2, 3]]
                    signal[i, j - i - 1, 1] = value
                j += 1
                if j < num_fends:
                    fend2 = mapping[j]
            while j < max_fend[i]:
                fend2 = mapping[j]
                if fends[fend2, 0] > 0:
                    value = 1.0
                    value *= map_mat[map_bin, fends[fend2, 1]]
                    value *= gc_mat[gc_bin, fends[fend2, 2]]
                    value *= len_mat[len_bin, fends[fend2, 3]]
                    signal[i, j - i - 1, 1] = value
                j += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bin_bounds_cis(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int64_t, ndim=1] indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] bins not None,
        np.ndarray[DTYPE_int_t, ndim=2] fends not None,
        np.ndarray[DTYPE_t, ndim=2] map_matrix not None,
        np.ndarray[DTYPE_t, ndim=2] gc_matrix not None,
        np.ndarray[DTYPE_t, ndim=2] len_matrix not None,
        np.ndarray[DTYPE_int_t, ndim=1] observed not None,
        np.ndarray[DTYPE_t, ndim=1] expected not None,
        int start_fend,
        int stop_fend,
        int num_bins):
    cdef long long int i, j, fend1, fend2, bin1, bin2, index
    cdef double value
    cdef long long int num_fends = bins.shape[0]
    with nogil:
        # find observed bin values
        for i in range(num_fends - 2):
            fend1 = start_fend + i
            bin1 = bins[i]
            if fends[fend1, 0] == 0 or bin1 < 0:
                continue
            index = bin1 * num_bins - bin1 * (bin1 + 1) / 2 - bin1 - 1
            j = indices[i]
            while j < indices[i + 1] and data[j, 1] < stop_fend:
                fend2 = data[j, 1]
                j += 1
                bin2 = bins[fend2 - start_fend]
                if fends[fend2, 0] == 0 or bin2 <= bin1:
                    continue
                observed[index + bin2] += data[j, 2]
        # find expected bin values
        for i in range(num_fends - 2):
            fend1 = start_fend + i
            bin1 = bins[i]
            if fends[fend1, 0] == 0 or bin1 < 0:
                continue
            index = bin1 * num_bins - bin1 * (bin1 + 1) / 2 - bin1 - 1
            fend2 = fend1 + 2
            bin2 = bins[fend2 - start_fend]
            if fends[fend2, 0] != 0 and bin2 > bin1 and bin2 > 0:
                value = 1.0
                value *= map_matrix[fends[fend1, 1], fends[fend2, 1]]
                value *= gc_matrix[fends[fend1, 2], fends[fend2, 2]]
                value *= len_matrix[fends[fend1, 3], fends[fend2, 3]]
                expected[index + bin2] += value
            for j in range((i / 2 + 2) * 2, num_fends):
                fend2 = j + start_fend
                bin2 = bins[fend2 - start_fend]
                if fends[fend2, 0] == 0 or bin2 < 0:
                    continue
                if bin2 <= bin1:
                    continue
                value = 1.0
                value *= map_matrix[fends[fend1, 1], fends[fend2, 1]]
                value *= gc_matrix[fends[fend1, 2], fends[fend2, 2]]
                value *= len_matrix[fends[fend1, 3], fends[fend2, 3]]
                expected[index + bin2] += value
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_HiCPipe_corrected_counts(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int_t, ndim=2] fends not None,
        np.ndarray[DTYPE_t, ndim=2] map_matrix not None,
        np.ndarray[DTYPE_t, ndim=2] gc_matrix not None,
        np.ndarray[DTYPE_t, ndim=2] len_matrix not None,
        np.ndarray[DTYPE_t, ndim=1] corrected not None):
    cdef long long int i, j, fend1, fend2
    cdef double value
    cdef long long int num_data = data.shape[0]
    with nogil:
        # find observed bin values
        for i in range(num_data):
            fend1 = data[i, 0]
            fend2 = data[i, 1]
            if fends[fend1, 0] == 0 or fends[fend2, 0] == 0:
                continue
            value = 1.0
            value *= map_matrix[fends[fend1, 1], fends[fend2, 1]]
            value *= gc_matrix[fends[fend1, 2], fends[fend2, 2]]
            value *= len_matrix[fends[fend1, 3], fends[fend2, 3]]
            corrected[i] = data[i, 2] / value
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_fend_distance_bin_values(
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int64_t, ndim=1] data_indices not None,
        np.ndarray[DTYPE_t, ndim=1] corrected not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=1] distance_bins not None,
        np.ndarray[DTYPE_64_t, ndim=1] bin_sums not None,
        np.ndarray[DTYPE_int64_t, ndim=1] bin_counts not None,
        int fend,
        int stopfend):
    cdef long long int i, j, distance, fend2
    cdef long long int num_bins = distance_bins.shape[0]
    with nogil:
        j = 0
        for i in range(data_indices[fend], data_indices[fend + 1]):
            fend2 = data[i, 1]
            if filter[fend2] == 0 or fend2 >= stopfend:
                continue
            distance = mids[fend2] - mids[fend]
            while j < num_bins - 1 and distance > distance_bins[j]:
                j += 1
            bin_sums[j] += corrected[i]
        fend2 = fend + 2
        j = 0
        if fend2 < stopfend and filter[fend2] == 1:
            distance = mids[fend2] - mids[fend]
            while j < num_bins - 1 and distance > distance_bins[j]:
                j += 1
            bin_counts[j] += 1
        for fend2 in range((fend / 2 + 2) * 2, stopfend):
            if filter[fend2] == 0:
                continue
            distance = mids[fend2] - mids[fend]
            while j < num_bins - 1 and distance > distance_bins[j]:
                j += 1
            bin_counts[j] += 1
    return None