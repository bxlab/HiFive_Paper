# normalization_functions.pyx
# cython: profile=False

import cython
cimport numpy as np
import numpy

ctypedef np.float32_t DTYPE_t
ctypedef np.float64_t DTYPE_64_t
ctypedef np.int32_t DTYPE_int_t
ctypedef np.int64_t DTYPE_int64_t
ctypedef np.uint32_t DTYPE_uint_t
ctypedef np.int8_t DTYPE_int8_t


cdef extern from "math.h":
	double exp( double x ) nogil
	double log( double x ) nogil
	double log10( double x ) nogil
	double sqrt( double x ) nogil
	double pow( double x, double x ) nogil
	double abs( double x) nogil
	double round( double x) nogil
	double floor( double x) nogil
	double ceil( double x) nogil


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def BinCounts(
	np.ndarray[DTYPE_int_t, ndim=1] binned_counts not None,  
	np.ndarray[DTYPE_int_t, ndim=2] data not None,  
	np.ndarray[DTYPE_int64_t, ndim=1] indices not None,  
	np.ndarray[DTYPE_t, ndim=2] features not None,  
	np.ndarray[DTYPE_int_t, ndim=1] mapping not None,  
	int start,
	int bin_size,
	int num_bins ):
	cdef int i, fend1, fend2, bin1, mapped_bin1, index, bin2, mapped_bin2
	cdef int num_fends = features.shape[0]
	with nogil:
		for fend1 in range(num_fends - 1):
			bin1 = int(floor((features[fend1, 0] - start) / bin_size))
			if mapping[bin1] == -1:
				continue
			mapped_bin1 = mapping[bin1]
			index = num_bins * mapped_bin1 - mapped_bin1 * (mapped_bin1 + 1) / 2 - mapped_bin1 - 1
			for i in range(indices[fend1], indices[fend1 + 1]):
				fend2 = data[i, 1]
				bin2 = int(floor((features[fend2, 0] - start) / bin_size))
				if mapping[bin2] == -1:
					continue
				mapped_bin2 = mapping[bin2]
				if mapped_bin2 != mapped_bin1:
					binned_counts[index + mapped_bin2] += data[i, 2]
	return None
			
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def BinTransCounts(
	np.ndarray[DTYPE_int_t, ndim=2] binned_counts not None,  
	np.ndarray[DTYPE_int_t, ndim=2] data not None,  
	np.ndarray[DTYPE_int64_t, ndim=1] indices not None,  
	np.ndarray[DTYPE_t, ndim=2] features1 not None,  
	np.ndarray[DTYPE_t, ndim=2] features2 not None,  
	np.ndarray[DTYPE_int_t, ndim=1] mapping1 not None,  
	np.ndarray[DTYPE_int_t, ndim=1] mapping2 not None,  
	int start1,
	int start2,
	int bin_size):
	cdef int i, fend1, fend2, bin1, mapped_bin1, bin2, mapped_bin2
	cdef int num_fends1 = features1.shape[0]
	cdef int num_fends2 = features2.shape[0]
	with nogil:
		for fend1 in range(num_fends1):
			bin1 = int(floor((features1[fend1, 0] - start1) / bin_size))
			if mapping1[bin1] == -1:
				continue
			mapped_bin1 = mapping1[bin1]
			for i in range(indices[fend1], indices[fend1 + 1]):
				fend2 = data[i, 1]
				bin2 = int(floor((features2[fend2, 0] - start2) / bin_size))
				if mapping2[bin2] == -1:
					continue
				mapped_bin2 = mapping2[bin2]
				binned_counts[mapped_bin1, mapped_bin2] += data[i, 2]
	return None
			

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def FindUpperTriFeatures(
	np.ndarray[DTYPE_t, ndim=1] feature_vector not None,  
	np.ndarray[DTYPE_t, ndim=1] features not None):
	cdef int i, j, k
	cdef int num_bins = features.shape[0]
	with nogil:
		k = 0
		for i in range(num_bins - 1):
			for j in range(i + 1, num_bins):
				feature_vector[k] = log(features[i] * features[j])
				k += 1
	return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def FindExpectedCounts(
	np.ndarray[DTYPE_t, ndim=1] len_vec not None,  
	np.ndarray[DTYPE_t, ndim=1] gc_vec not None,  
	np.ndarray[DTYPE_t, ndim=1] map_vec not None,  
	np.ndarray[DTYPE_t, ndim=1] expected not None,  
	np.ndarray[DTYPE_t, ndim=1] coeff not None):
	cdef int i
	cdef int num_bins = expected.shape[0]
	with nogil:
		for i in range(num_bins):
			expected[i] = exp(coeff[0] + coeff[1] * len_vec[i] + coeff[2] * gc_vec[i] + map_vec[i])
	return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def FindTransExpectedCounts(
	np.ndarray[DTYPE_t, ndim=2] len_mat not None,  
	np.ndarray[DTYPE_t, ndim=2] gc_mat not None,  
	np.ndarray[DTYPE_t, ndim=2] map_mat not None,  
	np.ndarray[DTYPE_t, ndim=2] expected not None,  
	np.ndarray[DTYPE_t, ndim=1] coeff not None):
	cdef int i, j
	cdef int xdim = expected.shape[0]
	cdef int ydim = expected.shape[1]
	with nogil:
		for i in range(xdim):
			for j in range(ydim):
				expected[i, j] = exp(coeff[0] + coeff[1] * len_mat[i, j] + coeff[2] * gc_mat[i, j] + map_mat[i, j])
	return None
