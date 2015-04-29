#!/usr/bin/env python

"""A class for handling HiC analysis."""

import os
import sys

import numpy
import h5py
from scipy.optimize import fmin_l_bfgs_b as bfgs
try:
    from mpi4py import MPI
except:
    pass
try:
    import mlpy
except:
    pass

import hic_binning
import libraries._hic_binning as _binning
import libraries._hic_distance as _distance
import libraries._hic_interactions as _interactions
import libraries._hic_optimize as _optimize
import plotting


class HiC(object):

    """
    This is the class for handling HiC analysis.

    This class relies on :class:`Fend <hifive.fend.Fend>` and :class:`HiCData <hifive.hic_data.HiCData>` for genomic position and interaction count data. Use this class to perform filtering of fends based on coverage, model fend bias and distance dependence, and downstream analysis and manipulation. This includes binning of data, plotting of data, modeling of data, and statistical analysis.

        .. note::
          This class is also available as hifive.HiC

    When initialized, this class creates an h5dict in which to store all data associated with this object.

    :param filename: The file name of the h5dict. This should end with the suffix '.hdf5'
    :type filename: str.
    :param mode: The mode to open the h5dict with. This should be 'w' for creating or overwriting an h5dict with name given in filename.
    :type mode: str.
    :param silent: Indicates whether to print information about function execution for this object.
    :type silent: bool.
    :returns: :class:`HiC <hifive.hic.HiC>` class object.

    :attributes: * **file** (*str.*) - A string containing the name of the file passed during object creation for saving the object to.
                 * **silent** (*bool.*) - A boolean indicating whether to suppress all of the output messages.
                 * **history** (*str.*) - A string containing all of the commands executed on this object and their outcome.
                 * **normalization** (*str.*) - A string stating which type of normalization has been performed on this object. This starts with the value 'none'.
                 * **comm** (*class*) - A link to the MPI.COMM_WORLD class from the mpi4py package. If this package isn't present, this is set to 'None'.
                 * **rank** (*int.*) - The rank integer of this process, if running with mpi, otherwise set to zero.
                 * **num_procs** (*int.*) - The number of processes being executed in parallel. If mpi4py package is not present, this is set to one.

    In addition, many other attributes are initialized to the 'None' state.
    """

    def __init__(self, filename, mode='r', silent=False):
        """Create a HiC object."""
        self.file = os.path.abspath(filename)
        if 'mpi4py' in sys.modules.keys():
            self.comm = MPI.COMM_WORLD
            self.rank = self.comm.Get_rank()
            self.num_procs = self.comm.Get_size()
        else:
            self.comm = None
            self.rank = 0
            self.num_procs = 1
        if self.rank == 0:
            self.silent = silent
        else:
            self.silent = True
        self.binning_corrections = None
        self.binning_correction_indices = None
        self.binning_fend_indices = None
        self.binning_num_bins = None
        self.model_parameters = None
        self.corrections = None
        self.distance_parameters = None
        self.bin_distance_parameters = None
        self.chromosome_means = None
        self.normalization = 'none'
        self.history = ''
        if mode != 'w':
            self.load()
        return None

    def __getitem__(self, key):
        """Dictionary-like lookup."""
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            return None

    def __setitem__(self, key, value):
        """Dictionary-like value setting."""
        self.__dict__[key] = value
        return None

    def load_data(self, filename):
        """
        Load fend-pair counts and fend object from :class:`HiCData <hifive.hic_data.HiCData>` object.

        :param filename: Specifies the file name of the :class:`HiCData <hifive.hic_data.HiCData>` object to associate with this analysis.
        :type filename: str.
        :returns: None

        :Attributes: * **datafilename** (*str.*) - A string containing the relative path of the HiCData file.
                     * **fendfilename** (*str.*) - A string containing the relative path of the Fend file associated with the HiCData file.
                     * **fends** (*filestream*) - A filestream to the hdf5 Fragment file such that all saved Fend attributes can be accessed through this class attribute.
                     * **data** (*filestream*) - A filestream to the hdf5 FiveCData file such that all saved HiCData attributes can be accessed through this class attribute.
                     * **chr2int** (*dict.*) - A dictionary that converts chromosome names to chromosome indices.
                     * **filter** (*ndarray*) - A numpy array of type int32 and size N where N is the number of fends. This contains the inclusion status of each fend with a one indicating included and zero indicating excluded and is initialized with all fends included.

        When a HiCData object is associated with the project file, the 'history' attribute is updated with the history of the HiCData object.
        """
        self.history += "HiC.load_data(filename='%s') - " % (filename)
        if self.rank > 0:
            return None
        filename = os.path.abspath(filename)
        # ensure data h5dict exists
        if not os.path.exists(filename):
            if not self.silent:
                print >> sys.stderr, ("Could not find %s. No data loaded.\n") % (filename),
            self.history += "Error: '%s' not found\n" % filename
            return None
        self.datafilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(filename)),
                                       os.path.dirname(self.file)), os.path.basename(filename))
        self.data = h5py.File(filename, 'r')
        self.history = self.data['/'].attrs['history'] + self.history
        fendfilename = self.data['/'].attrs['fendfilename']
        if fendfilename[:2] == './':
            fendfilename = fendfilename[2:]
        parent_count = fendfilename.count('../')
        fendfilename = '/'.join(os.path.abspath(filename).split('/')[:-(1 + parent_count)] +
                                fendfilename.lstrip('/').split('/')[parent_count:])
        self.fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(fendfilename),
                                       os.path.dirname(self.file)), os.path.basename(fendfilename))
        # ensure fend h5dict exists
        if not os.path.exists(fendfilename):
            if not self.silent:
                print >> sys.stderr, ("Could not find %s.\n") % (fendfilename),
            self.history += "Error: '%s' not found\n" % fendfilename
            return None
        self.fends = h5py.File(fendfilename, 'r')
        # create dictionary for converting chromosome names to indices
        self.chr2int = {}
        for i, chrom in enumerate(self.fends['chromosomes']):
            self.chr2int[chrom] = i
        # create arrays
        self.filter = numpy.ones(self.fends['fends'].shape[0], dtype=numpy.int32)
        self.history += "Succcess\n"
        return None

    def save(self, out_fname=None):
        """
        Save analysis parameters to h5dict.

        :param filename: Specifies the file name of the :class:`HiC <hifive.hic.HiC>` object to save this analysis to.
        :type filename: str.
        :returns: None
        """
        self.history.replace("'None'", "None")
        if self.rank > 0:
            return None
        if not out_fname is None:
            original_file = os.path.abspath(self.file)
            out_fname = os.path.abspath(out_fname)
            if 'datafilename' in self.__dict__:
                datafilename = self.datafilename
                if datafilename[:2] == './':
                    datafilename = datafilename[2:]
                parent_count = datafilename.count('../')
                datafilename = '/'.join(original_file.split('/')[:-(1 + parent_count)] +
                                        datafilename.lstrip('/').split('/')[parent_count:])
                datafilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(datafilename)),
                                               os.path.dirname(out_fname)), os.path.basename(datafilename))
            if 'fendfilename' in self.__dict__:
                fendfilename = self.fendfilename
                if fendfilename[:2] == './':
                    fendfilename = fendfilename[2:]
                parent_count = fendfilename.count('../')
                fendfilename = '/'.join(original_file.split('/')[:-(1 + parent_count)] +
                                        fendfilename.lstrip('/').split('/')[parent_count:])
                fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                               os.path.dirname(out_fname)), os.path.basename(fendfilename))
        else:
            out_fname = self.file
            if 'datafilename' in self.__dict__:
                datafilename = self.datafilename
            if 'fendfilename' in self.__dict__:
                fendfilename = self.fendfilename
        datafile = h5py.File(out_fname, 'w')
        for key in self.__dict__.keys():
            if key in ['data', 'fends', 'file', 'chr2int', 'comm', 'rank', 'num_procs', 'silent']:
                continue
            elif self[key] is None:
                continue
            elif key == 'fendfilename':
                datafile.attrs[key] = fendfilename
            elif key == 'datafilename':
                datafile.attrs[key] = datafilename
            elif isinstance(self[key], numpy.ndarray):
                datafile.create_dataset(key, data=self[key])
            elif not isinstance(self[key], dict):
                datafile.attrs[key] = self[key]
        datafile.close()
        return None

    def load(self):
        """
        Load analysis parameters from h5dict specified at object creation and open h5dicts for associated :class:`HiCData <hifive.hic_data.HiCData>` and :class:`Fend <hifive.fend.Fend>` objects.

        Any call of this function will overwrite current object data with values from the last :func:`save` call.

        :returns: None
        """
        # set parameters to init state
        self.binning_corrections = None
        self.binning_correction_indices = None
        self.binning_fend_indices = None
        self.binning_num_bins = None
        self.model_parameters = None
        self.corrections = None
        self.distance_parameters = None
        self.bin_distance_parameters = None
        self.chromosome_means = None
        self.normalization = 'none'
        self.history = ''
        # load data hdf5 dict 
        datafile = h5py.File(self.file, 'r')
        for key in datafile.keys():
            self[key] = numpy.copy(datafile[key])
        for key in datafile['/'].attrs.keys():
            self[key] = datafile['/'].attrs[key]
        # ensure data h5dict exists
        if 'datafilename' in self.__dict__:
            datafilename = self.datafilename
            if datafilename[:2] == './':
                datafilename = datafilename[2:]
            parent_count = datafilename.count('../')
            datafilename = '/'.join(self.file.split('/')[:-(1 + parent_count)] +
                                datafilename.lstrip('/').split('/')[parent_count:])
            if not os.path.exists(datafilename):
                if not self.silent:
                    print >> sys.stderr, ("Could not find %s. No data loaded.\n") % (datafilename),
            else:
                self.data = h5py.File(datafilename, 'r')
        # ensure fend h5dict exists
        if 'fendfilename' in self.__dict__:
            fendfilename = self.fendfilename
            if fendfilename[:2] == './':
                fendfilename = fendfilename[2:]
            parent_count = fendfilename.count('../')
            fendfilename = '/'.join(self.file.split('/')[:-(1 + parent_count)] +
                                fendfilename.lstrip('/').split('/')[parent_count:])
            if not os.path.exists(fendfilename):
                if not self.silent:
                    print >> sys.stderr, ("Could not find %s. No fends loaded.\n") % (fendfilename),
            else:
                self.fends = h5py.File(fendfilename, 'r')
        # create dictionary for converting chromosome names to indices
        self.chr2int = {}
        for i, chrom in enumerate(self.fends['chromosomes']):
            self.chr2int[chrom] = i
        datafile.close()
        return None

    def reset_filter(self):
        """
        Return all fends to a valid filter state.

        :returns: None
        """
        self.history += "HiC.reset_filter() - "
        self.filter.fill(1)
        self.history += "Success\n"
        return None

    def filter_fends(self, mininteractions=10, mindistance=0, maxdistance=0):
        """
        Iterate over the dataset and remove fends that do not have 'minobservations' within 'maxdistance' of themselves using only unfiltered fends.

        In order to create a set of fends that all have the necessary number of interactions, after each round of filtering, fend interactions are retallied using only interactions that have unfiltered fends at both ends.

        :param mininteractions: The required number of interactions for keeping a fend in analysis.
        :type mininteractions: int.
        :param mindistance: The minimum inter-fend distance used to count fend interactions.
        :type mindistance: int.
        :param maxdistance: The maximum inter-fend distance used to count fend interactions. A value of 0 indicates no maximum should be used.
        :type maxdistance: int.
        :returns: None
        """
        self.history += "HiC.filter_fends(mininteractions=%i, mindistance=%s, maxdistance=%s) - " % (mininteractions, str(mindistance), str(maxdistance))
        if not self.silent:
            print >> sys.stderr, ("Filtering fends..."),
        self.mininteractions = mininteractions
        self.mindistance = mindistance
        if maxdistance is None:
            maxdistance = 0
        self.maxdistance = maxdistance
        original_count = numpy.sum(self.filter)
        previous_valid = original_count + 1
        current_valid = original_count
        coverage = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
        # determine maximum ranges of valid interactions for each fend
        max_fend = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
        min_fend = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
        _interactions.find_max_fend(max_fend,
                                    self.fends['fends']['mid'][...],
                                    self.fends['fends']['chr'][...],
                                    self.fends['chr_indices'][...],
                                    0,
                                    maxdistance)
        _interactions.find_min_fend(min_fend,
                                    self.fends['fends']['mid'][...],
                                    self.fends['fends']['chr'][...],
                                    self.fends['chr_indices'][...],
                                    mindistance)
         # copy needed arrays
        data = self.data['cis_data'][...]
        indices = self.data['cis_indices'][...]
        # repeat until all remaining fends have mininteraction valid interactions
        while current_valid < previous_valid:
            previous_valid = current_valid
            coverage.fill(0)
            _interactions.find_fend_coverage(data,
                                             indices,
                                             self.filter,
                                             min_fend,
                                             max_fend,
                                             coverage,
                                             mininteractions)
            current_valid = numpy.sum(self.filter)
        if not self.silent:
            print >> sys.stderr, ("Removed %i of %i fends\n") % (original_count - current_valid, original_count),
        self.history += "Success\n"
        return None

    def find_distance_parameters(self, numbins=90, minsize=200, maxsize=0, corrected=False):
        """
        Count reads and possible interactions from valid fend pairs in each distance bin to find mean bin signals. This function is MPI compatible.

        This partitions the range of interaction distances (measured from mipoints of the involved fends) from the 'minsize' to 'maxsize' into a number of partitions equal to 'numbins'. The first bin contains all distances less than or equal to 'minsize'. The remaining bins are defined such that their log ranges are equal to one another. The curve defined by the mean interaction value of each bin can be smoothed using a  triangular smoothing operation.

        :param numbins: The number of bins to divide the distance range into. The first bin extends from zero to 'minsize', while the remaining bins are divided into evenly-spaced log-sized bins from 'minsize' to 'maxsize' or the maximum inter-fend distance, whichever is greater.
        :type numbins: int.
        :param minsize: The upper size limit of the smallest distance bin.
        :type minsize: int.
        :param maxsize: If this value is larger than the largest included chromosome, it will extend bins out to maxsize. If this value is smaller, it is ignored.
        :type maxsize: int.
        :param corrected: If True, correction values are applied to counts prior to summing.
        :type corrected: bool.
        :returns: None

        :Attributes: * **distance_parameters** (*ndarray*) - A numpy array of type float32 and size of N x 3 where N is one less than the number of distance bins containing at least one valid observation out of the 'numbins' number of bins that the distance range was divided into. The First column contains upper distance cutoff for each bin, the second column contains the slope associated with each bin line segment, and the third column contains the line segment intercepts. Line segments describe the relationship of  observation counts versus distance.
                     * **bin_distance_parameters** (*ndarray*) - A numpy array of type float32 and size of N x 3 where N is one less than the number of distance bins containing at least one valid observation out of the 'numbins' number of bins that the distance range was divided into. The First column contains upper distance cutoff for each bin, the second column contains the slope associated with each bin line segment, and the third column contains the line segment intercepts. Line segments describe the relationship of binary observations versus distance.
                     * **chromosome_means** (*ndarray*) - A numpy array of type float32 and length equal to the number of chromosomes. This is initialized to zeros until fend correction values are found.
        """
        if numbins == 0:
            return None
        self.history += "HiC.find_distance_parameters(numbins=%i, minsize=%i, maxsize=%s, corrected=%s) - " % (numbins, minsize, str(maxsize), corrected)
        if not self.silent:
            print >> sys.stderr, ('Finding distance arrays...'),
        # determine max distance range of data
        chr_indices = self.fends['chr_indices'][...]
        if maxsize is None:
            maxsize = 0
        max_dist = 0
        valid_chroms = []
        for i in range(chr_indices.shape[0] - 1):
            start_fend = chr_indices[i]
            stop_fend = chr_indices[i + 1]
            while start_fend < stop_fend and self.filter[start_fend] == 0:
                start_fend += 1
            while stop_fend > start_fend and self.filter[stop_fend - 1] == 0:
                stop_fend -= 1
            if stop_fend - 1 > start_fend:
                valid_chroms.append(i)
                max_dist = max(max_dist,
                               self.fends['fends']['mid'][stop_fend - 1] - self.fends['fends']['mid'][start_fend])
        valid_chroms = numpy.asarray(valid_chroms, dtype=numpy.int32)
        # create cutoff values evenly spaced in log space
        cutoffs = numpy.linspace(numpy.log(max(minsize, 1.0)), numpy.log(max(maxsize, max_dist)),
                                 numbins).astype(numpy.float32)
        cutoffs[-1] += 1.0
        bin_size = numpy.zeros((numbins, 2), dtype=numpy.int64)
        count_sum = numpy.zeros(numbins, dtype=numpy.float64)
        logdistance_sum = numpy.zeros(numbins, dtype=numpy.float64)
        # for each chromosome, find counts, possible interactions, and distance sums for each bin
        for h, i in enumerate(valid_chroms):
            start_fend = chr_indices[i]
            stop_fend = chr_indices[i + 1]
            rev_mapping = numpy.where(self.filter[start_fend:stop_fend] == 1)[0].astype(numpy.int32)
            num_valid = rev_mapping.shape[0]
            if not self.silent:
                print >> sys.stderr, ('\r%s\rFinding distances for chromosome %s...') % \
                                     (' ' * 80, self.fends['chromosomes'][i]),
            # partition total possible interactions into roughly even-sized groups to spread across nodes
            total_pairs = num_valid * (num_valid - 1) / 2
            node_cutoffs = numpy.linspace(0, total_pairs, self.num_procs + 1).astype(numpy.float64)
            pair_sums = numpy.r_[0, num_valid - numpy.arange(1, num_valid + 1)].astype(numpy.int64)
            for j in range(1, pair_sums.shape[0]):
                pair_sums[j] += pair_sums[j - 1]
            node_ranges = numpy.searchsorted(pair_sums, node_cutoffs, side='left')
            node_start = node_ranges[self.rank]
            if (node_start > 0 and node_cutoffs[self.rank] - pair_sums[node_start - 1] <
                                   pair_sums[node_start] - node_cutoffs[self.rank]):
                node_start -= 1
            node_stop = node_ranges[self.rank + 1]
            if (node_stop > 0 and node_cutoffs[self.rank + 1] - pair_sums[node_stop - 1] <
                                   pair_sums[node_stop] - node_cutoffs[self.rank + 1]):
                node_stop -= 1
            mapping = numpy.zeros(stop_fend - start_fend, dtype=numpy.int32) - 1
            mapping[rev_mapping] = numpy.arange(num_valid, dtype=numpy.int32)
            mids = self.fends['fends']['mid'][rev_mapping + start_fend]
            # pull relevant data
            start_index = self.data['cis_indices'][rev_mapping[node_start] + start_fend]
            stop_index = self.data['cis_indices'][max(rev_mapping[node_stop - 1] + 1,
                                                      rev_mapping[node_start]) + start_fend]
            indices = self.data['cis_data'][start_index:stop_index, :]
            counts = indices[:, 2].astype(numpy.float64)
            if corrected:
                counts /= (self.corrections[indices[:, 0]] * self.corrections[indices[:, 1]])
            indices = indices[:, :2]
            indices -= start_fend
            # find bin sums
            _distance.find_distance_bin_sums(mapping,
                                             rev_mapping,
                                             cutoffs,
                                             mids,
                                             counts,
                                             indices,
                                             bin_size,
                                             count_sum,
                                             logdistance_sum,
                                             node_start,
                                             node_stop)
        if self.rank == 0:
            # exchange arrays
            if not self.silent:
                print >> sys.stderr, ('\r%s\rCalculating distance function...') % (' ' * 80),
            for i in range(1, self.num_procs):
                bin_size += self.comm.recv(source=i, tag=11)
                count_sum += self.comm.recv(source=i, tag=11)
                logdistance_sum += self.comm.recv(source=i, tag=11)
            valid = numpy.where(count_sum > 0)[0]
            count_means = numpy.log(count_sum[valid].astype(numpy.float64) / bin_size[valid, 1])
            binary_means = numpy.log(bin_size[valid, 0].astype(numpy.float64) / bin_size[valid, 1])
            distance_means = logdistance_sum[valid] / bin_size[valid, 1]
            # find distance line parameters, cutoffs, slopes and intercepts
            distance_parameters = numpy.zeros((valid.shape[0] - 1, 3), dtype=numpy.float32)
            distance_parameters[:-1, 0] = distance_means[1:-1]
            distance_parameters[-1, 0] = numpy.inf
            distance_parameters[:, 1] = ((count_means[1:] - count_means[:-1]) /
                                         (distance_means[1:] - distance_means[:-1]))
            distance_parameters[:, 2] = (count_means[1:] - distance_parameters[:, 1] * distance_means[1:])
            bin_distance_parameters = numpy.zeros((valid.shape[0] - 1, 3), dtype=numpy.float32)
            bin_distance_parameters[:-1, 0] = distance_means[1:-1]
            bin_distance_parameters[-1, 0] = numpy.inf
            bin_distance_parameters[:, 1] = ((binary_means[1:] - binary_means[:-1]) /
                                         (distance_means[1:] - distance_means[:-1]))
            bin_distance_parameters[:, 2] = (binary_means[1:] - bin_distance_parameters[:, 1] * distance_means[1:])
            # distribute distance parameters to all nodes
            for i in range(1, self.num_procs):
                self.comm.send(distance_parameters, dest=i, tag=11)
                self.comm.send(bin_distance_parameters, dest=i, tag=11)
        else:
            self.comm.send(bin_size, dest=0, tag=11)
            self.comm.send(count_sum, dest=0, tag=11)
            self.comm.send(logdistance_sum, dest=0, tag=11)
            distance_parameters = self.comm.recv(source=0, tag=11)
            bin_distance_parameters = self.comm.recv(source=0, tag=11)
        self.distance_parameters = distance_parameters
        self.bin_distance_parameters = bin_distance_parameters
        if self.chromosome_means is None:
            self.chromosome_means = numpy.zeros(self.fends['chr_indices'].shape[0] - 1, dtype=numpy.float32)
        if not self.silent:
            print >> sys.stderr, ('\r%s\rFinding distance curve... Done\n') % (' ' * 80),
        self.history += "Success\n"
        return None

    def find_probability_fend_corrections(self, mindistance=0, maxdistance=0, minchange=0.0001,
                                          max_iterations=1000, learningstep=0.5, chroms=[], precalculate=True,
                                          precorrect=False):
        """
        Using gradient descent, learn correction values for each valid fend based on a Poisson distribution of observations. This function is MPI compatible.

        :param mindistance: The minimum inter-fend distance to be included in modeling.
        :type mindistance: int.
        :param maxdistance: The maximum inter-fend distance to be included in modeling.
        :type maxdistance: int.
        :param minchange: The cutoff threshold for early learning termination for the maximum absolute gradient value.
        :type minchange: float
        :param max_iterations: The maximum number of iterations to carry on gradient descent for.
        :type max_iterations: int.
        :param learningstep: The scaling factor for decreasing learning rate by if step doesn't meet armijo criterion.
        :type learningstep: float
        :param chroms: A list of chromosomes to calculate corrections for. If set as None, all chromosome corrections are found.
        :type chroms: list
        :param precalculate: Specifies whether the correction values should be initialized at the fend means.
        :type precalculate: bool.
        :param precorrect: Use binning-based corrections in expected value calculations, resulting in a chained normalization approach.
        :type precorrect: bool.
        :returns: None

        :Attributes: * **corrections** (*ndarray*) - A numpy array of type float32 and length equal to the number of fends. All invalid fends have an associated correction value of zero.

        The 'normalization' attribute is updated to 'probability' or 'binning-probability', depending on if the 'precorrect' option is selected. In addition, the 'chromosome_means' attribute is updated such that the mean correction (sum of all valid chromosomal correction value pairs) is adjusted to zero and the corresponding chromosome mean is adjusted the same amount but the opposite sign. 
        """
        self.history += "HiC.find_probability_fend_corrections(mindistance=%i, maxdistance=%s, minchange=%f, max_iterations=%i, learningstep=%f, chroms=%s, precalculate=%s, precorrect=%s) - " % (mindistance, str(maxdistance), minchange, max_iterations, learningstep, str(chroms), precalculate, precorrect)
        if precorrect and self.binning_corrections is None:
            if not self.silent:
                print >> sys.stderr, ("Precorrection can only be used in project has previously run 'find_binning_fend_corrections'.\n"),
            self.history += "Error: 'find_binning_fend_corrections()' not run yet\n"
            return None
        # make sure distance parameters have been calculated
        if self.distance_parameters is None:
            if not self.silent:
                print >> sys.stderr, ("This normalization requires a project that has already had 'find_distance_parameters' run.\n"),
            self.history += "Error: 'find_distance_parameters()' not run yet\n"
            return None
        if self.corrections is None:
            self.corrections = numpy.ones(self.fends['fends'].shape[0], dtype=numpy.float32)
        if (chroms is None or
                (isinstance(chroms, list) and
                (len(chroms) == 0 or
                (len(chroms) == 1 and chroms[0] == ''))) or
                chroms == ''):
            chroms = self.chr2int.keys()
        chrints = numpy.zeros(len(chroms), dtype=numpy.int32)
        chr_indices = self.fends['chr_indices'][...]
        for i in range(len(chroms)):
            chrints[i] = self.chr2int[chroms[i]]
        chroms = list(numpy.array(chroms)[numpy.argsort(chrints)])
        if self.chromosome_means is None:
            self.chromosome_means = numpy.zeros(self.fends['chr_indices'].shape[0] - 1, dtype=numpy.float32)
        if maxdistance == 0 or maxdistance is None:
            maxdistance = 1999999999
        for chrom in chroms:
            chrint = self.chr2int[chrom]
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding fend correction arrays for chromosome %s...") %\
                    (' ' * 80, chrom),
            start_fend = chr_indices[chrint]
            stop_fend = chr_indices[chrint + 1]
            num_fends = stop_fend - start_fend
            rev_mapping = numpy.where(self.filter[start_fend:stop_fend] == 1)[0].astype(numpy.int32)
            num_valid = rev_mapping.shape[0]
            if num_valid == 0:
                continue
            node_ranges = numpy.round(numpy.linspace(0, num_valid, self.num_procs + 1)).astype(numpy.int32)
            mapping = numpy.zeros(num_fends, dtype=numpy.int32) - 1
            mapping[rev_mapping] = numpy.arange(num_valid).astype(numpy.int32)
            fend_ranges = numpy.zeros((num_valid, 3), dtype=numpy.int64)
            mids = self.fends['fends']['mid'][rev_mapping + start_fend]
            # find number of downstream interactions for each fend using distance limits
            _interactions.find_fend_ranges(rev_mapping,
                                           mids,
                                           fend_ranges,
                                           mindistance,
                                           maxdistance,
                                           node_ranges[self.rank],
                                           node_ranges[self.rank + 1],
                                           start_fend)
            if self.rank == 0:
                for i in range(1, self.num_procs):
                    fend_ranges[node_ranges[i]:node_ranges[i + 1], :] = self.comm.recv(source=i, tag=11)
                for i in range(1, self.num_procs):
                    self.comm.Send(fend_ranges, dest=i, tag=13)
            else:
                self.comm.send(fend_ranges[node_ranges[self.rank]:node_ranges[self.rank + 1], :], dest=0, tag=11)
                self.comm.Recv(fend_ranges, source=0, tag=13)
            # find range of fends for each node, creating as even spacing as possible
            total_pairs = numpy.sum(fend_ranges[:, 0])
            node_ranges = numpy.round(numpy.linspace(0, total_pairs, self.num_procs + 1)).astype(numpy.int64)
            temp = numpy.r_[0, fend_ranges[:, 0]]
            for i in range(1, temp.shape[0]):
                temp[i] += temp[i - 1]
            start = numpy.searchsorted(temp, node_ranges[self.rank])
            if start > 0:
                if node_ranges[self.rank] - temp[start - 1] < temp[start] - node_ranges[self.rank]:
                    start -= 1
            stop = numpy.searchsorted(temp, node_ranges[self.rank + 1])
            if stop > 0:
                if node_ranges[self.rank + 1] - temp[stop - 1] < temp[stop] - node_ranges[self.rank + 1]:
                    stop -= 1
            num_pairs = temp[stop] - temp[start]
            del temp
            # pull needed data for each node and determine number of nonzero pairs
            start_index = self.data['cis_indices'][start_fend + rev_mapping[start]]
            stop_index = max(start_index, self.data['cis_indices'][start_fend + rev_mapping[stop - 1] + 1])
            temp_data = self.data['cis_data'][start_index:stop_index, :]
            temp_data = temp_data[numpy.where(temp_data[:, 1] < stop_fend)[0], :]
            temp_data[:, 0] = mapping[temp_data[:, 0] - start_fend]
            temp_data[:, 1] = mapping[temp_data[:, 1] - start_fend]
            nonzero_pairs = _interactions.find_nonzeros_in_range(fend_ranges, temp_data)
            # allocate nonzero index and count arrays and fill and find number of interactions
            nonzero_indices0 = numpy.zeros(nonzero_pairs, dtype=numpy.int32)
            nonzero_indices1 = numpy.zeros(nonzero_pairs, dtype=numpy.int32)
            counts = numpy.zeros(nonzero_pairs, dtype=numpy.int32)
            interactions = numpy.zeros(rev_mapping.shape[0], dtype=numpy.int32)
            _interactions.find_nonzero_node_indices(fend_ranges,
                                                    nonzero_indices0,
                                                    nonzero_indices1,
                                                    counts,
                                                    temp_data,
                                                    interactions)
            del temp_data
            # allocate zero index arrays and fill and find number of interactions
            zero_indices0 = numpy.zeros(num_pairs - nonzero_pairs, dtype=numpy.int32)
            zero_indices1 = numpy.zeros(num_pairs - nonzero_pairs, dtype=numpy.int32)
            _interactions.find_zero_node_indices(rev_mapping,
                                                 fend_ranges,
                                                 nonzero_indices0,
                                                 nonzero_indices1,
                                                 zero_indices0,
                                                 zero_indices1,
                                                 interactions,
                                                 start,
                                                 stop,
                                                 start_fend)
            # find priors based on binary distance depedence function
            nonzero_means = numpy.zeros(nonzero_indices0.shape[0], dtype=numpy.float32)
            _distance.find_remapped_distance_means(nonzero_indices0,
                                                   nonzero_indices1,
                                                   mids,
                                                   nonzero_means,
                                                   self.bin_distance_parameters,
                                                   0.0)
            zero_means = numpy.zeros(zero_indices0.shape[0], dtype=numpy.float32)
            _distance.find_remapped_distance_means(zero_indices0,
                                                   zero_indices1,
                                                   mids,
                                                   zero_means,
                                                   self.bin_distance_parameters,
                                                   0.0)
            if self.rank == 0:
                for i in range(1, self.num_procs):
                    interactions += self.comm.recv(source=i, tag=11)
                temp = numpy.zeros(rev_mapping.shape[0], dtype=numpy.float64)
            else:
                self.comm.send(interactions, dest=0, tag=11)
                temp = None
            # if precorrecting using binning correction values, find correction matrices and adjust distance means
            if precorrect:
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rFinding binning corrections for chromosome %s...") % (' ' * 80,
                                                                                                          chrom),
                _optimize.find_binning_correction_adjustment(nonzero_means,
                                                             nonzero_indices0,
                                                             nonzero_indices1,
                                                             self.binning_corrections,
                                                             self.binning_num_bins,
                                                             self.binning_fend_indices[rev_mapping + start_fend, :])
                _optimize.find_binning_correction_adjustment(zero_means,
                                                             zero_indices0,
                                                             zero_indices1,
                                                             self.binning_corrections,
                                                             self.binning_num_bins,
                                                             self.binning_fend_indices[rev_mapping + start_fend, :])

            # if precalculating, find approximate corrections
            if precalculate:
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rPrecalculating corrections for chromosome %s...") %\
                        (' ' * 80, chrom),
                expected = numpy.zeros(rev_mapping.shape[0], dtype=numpy.float64)
                observed = numpy.zeros(rev_mapping.shape[0], dtype=numpy.float64)
                _interactions.sum_weighted_indices(nonzero_indices0, nonzero_indices1, nonzero_means, expected)
                _interactions.sum_weighted_indices(zero_indices0, zero_indices1, zero_means, expected)
                _interactions.sum_weighted_indices(nonzero_indices0, nonzero_indices1, None, observed)
                corrections = numpy.zeros(rev_mapping.shape[0], dtype=numpy.float32)
                if self.rank == 0:
                    temp1 = numpy.empty(expected.shape, dtype=expected.dtype)
                    for i in range(1, self.num_procs):
                        self.comm.Recv(temp1, source=i, tag=13)
                        expected += temp1
                        self.comm.Recv(temp1, source=i, tag=13)
                        observed += temp1
                    corrections[:] = numpy.minimum(100.0, numpy.maximum(0.01, observed / expected))
                    for i in range(1, self.num_procs):
                        self.comm.Send(corrections, dest=i, tag=13)
                else:
                    self.comm.Send(expected, dest=0, tag=13)
                    self.comm.Send(observed, dest=0, tag=13)
                    self.comm.Recv(corrections, source=0, tag=13)
            else:
                corrections = self.corrections[numpy.where(mapping >= 0)[0] + start_fend]
            new_corrections = numpy.copy(corrections)
            # calculate correction gradients
            gradients = numpy.zeros(corrections.shape[0], dtype=numpy.float64)
            cont = True
            if not self.silent:
                print >> sys.stderr, ("\r%s\rLearning corrections...") % (' ' * 80),            
            nonzero_means = numpy.log(nonzero_means)
            log_corrections = numpy.log(corrections).astype(numpy.float32)
            start_cost = _optimize.calculate_binom_cost(zero_indices0,
                                                        zero_indices1,
                                                        nonzero_indices0,
                                                        nonzero_indices1,
                                                        nonzero_means,
                                                        zero_means,
                                                        corrections,
                                                        log_corrections)
            if self.rank == 0:
                for i in range(1, self.num_procs):
                    start_cost += self.comm.recv(source=i, tag=11)
            else:
                self.comm.send(start_cost, dest=0, tag=11)
            previous_cost = start_cost
            change = 0.0
            cont = True
            iteration = 0
            while cont:
                gradients.fill(0.0)
                inv_corrections = (1.0 / corrections).astype(numpy.float32)
                _optimize.calculate_binom_gradients(zero_indices0,
                                                   zero_indices1,
                                                   nonzero_indices0,
                                                   nonzero_indices1,
                                                   zero_means,
                                                   corrections,
                                                   inv_corrections,
                                                   gradients)
                self._exchange_gradients(gradients, temp)
                if self.rank == 0:
                    gradients /= interactions
                    gradient_norm = numpy.sum(gradients ** 2.0)
                # find best step size
                armijo = numpy.inf
                t = 1.0
                while armijo > 0.0:
                    # if using multiple cores, pass gradients to root
                    if self.rank == 0:
                        new_corrections = numpy.minimum(100.0, numpy.maximum(0.01,
                                                        corrections - t * gradients)).astype(numpy.float32)
                        for i in range(1, self.num_procs):
                            self.comm.Send(new_corrections, dest=i, tag=13)
                    else:
                        self.comm.Recv(new_corrections, source=0, tag=13)
                    log_corrections = numpy.log(new_corrections)
                    cost = _optimize.calculate_binom_cost(zero_indices0,
                                                          zero_indices1,
                                                          nonzero_indices0,
                                                          nonzero_indices1,
                                                          nonzero_means,
                                                          zero_means,
                                                          new_corrections,
                                                          log_corrections)
                    if self.rank == 0:
                        for i in range(1, self.num_procs):
                            cost += self.comm.recv(source=i, tag=11)
                        if numpy.isnan(cost):
                            cost = numpy.inf
                            armijo = numpy.inf
                        else:
                            armijo = cost - previous_cost + t * gradient_norm
                        for i in range(1, self.num_procs):
                            self.comm.send(armijo, dest=i, tag=11)
                        if not self.silent:
                            print >> sys.stderr, ("\r%s iteration:%i cost:%f change:%f armijo: %f %s") %\
                                                 ('Learning corrections...', iteration, previous_cost,
                                                  change, armijo, ' ' * 20),
                    else:
                        self.comm.send(cost, dest=0, tag=11)
                        armijo = self.comm.recv(source=0, tag=11)
                    t *= learningstep
                previous_cost = cost
                corrections = new_corrections
                # find change
                if self.rank == 0:
                    change = numpy.amax(numpy.abs(gradients / corrections))
                    for i in range(1, self.num_procs):
                        self.comm.send(change, dest=i, tag=11)
                else:
                    change = self.comm.recv(source=0, tag=11)
                if not self.silent:
                    print >> sys.stderr, ("\r%s iteration:%i cost:%f change:%f %s") %\
                                         ('Learning corrections...', iteration, cost, change, ' ' * 40),
                iteration += 1
                if iteration >= max_iterations or change <= minchange:
                    cont = False
            # calculate chromosome mean
            chrom_mean = numpy.sum(corrections)
            chrom_mean = chrom_mean ** 2.0 - numpy.sum(corrections ** 2.0)
            chrom_mean /= corrections.shape[0] * (corrections.shape[0] - 1)
            self.chromosome_means[self.chr2int[chrom]] += numpy.log(chrom_mean)
            log_corrections = numpy.log(corrections).astype(numpy.float32)
            cost = _optimize.calculate_binom_cost(zero_indices0,
                                                 zero_indices1,
                                                 nonzero_indices0,
                                                 nonzero_indices1,
                                                 nonzero_means,
                                                 zero_means,
                                                 corrections,
                                                 log_corrections)
            self.corrections[rev_mapping + start_fend] = corrections / (chrom_mean ** 0.5)
            if self.rank == 0:
                for i in range(1, self.num_procs):
                    cost += self.comm.recv(source=i, tag=11)
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rLearning corrections... chromosome %s  Initial Cost:%f  Final Cost:%f  Done\n") % \
                        (' ' * 80, chrom, start_cost, cost),
            else:
                self.comm.send(cost, dest=0, tag=11)
            del counts, nonzero_indices0, nonzero_indices1, nonzero_means, zero_indices0, zero_indices1, zero_means
        if not self.silent:
            print >> sys.stderr, ("\rLearning corrections... Done%s\n") % (' ' * 80),
        if precorrect:
            self.normalization = 'binning-probability'
        else:
            self.normalization = 'probability'
        self.history += "Success\n"
        return None

    def _exchange_gradients(self, gradients, temp):
        if self.rank == 0:
            for i in range(1, self.num_procs):
                self.comm.Recv(temp, source=i, tag=13)
                gradients += temp
        else:
            self.comm.Send(gradients, dest=0, tag=13)
        return None

    def find_express_fend_corrections(self, iterations=100, mindistance=0, maxdistance=0, remove_distance=True, 
                                      usereads='cis', mininteractions=0, minchange=0.0001, chroms=[], precorrect=False,
                                      binary=False, kr=True):
        """
        Using iterative matrix-balancing approximation, learn correction values for each valid fend. This function is MPI compatible.

        :param iterations: The minimum number of iterations to use for learning fend corrections.
        :type iterations: int.
        :param mindistance: This is the minimum distance between fend midpoints needed to be included in the analysis. All possible and observed interactions with a distance shorter than this are ignored. If 'usereads' is set to 'trans', this value is ignored.
        :param maxdistance: The maximum inter-fend distance to be included in modeling. If 'usereads' is set to 'trans', this value is ignored.
        :type maxdistance: int.
        :type mindistance: int.
        :param remove_distance: Specifies whether the estimated distance-dependent portion of the signal is removed prior to learning fend corrections.
        :type remove_distance: bool.
        :param usereads: Specifies which set of interactions to use, 'cis', 'trans', or 'all'.
        :type usereads: str.
        :param mininteractions: If a non-zero 'mindistance' is specified or only 'trans' interactions are used, fend filtering will be performed again to ensure that the data being used is sufficient for analyzed fends. This parameter may specify how many interactions are needed for valid fends. If not given, the value used for the last call to :func:`filter_fends` is used or, barring that, one.
        :type mininteractions: int.
        :param minchange: The minimum mean change in fend correction parameter values needed to keep running past 'iterations' number of iterations. If using the Knight-Ruiz algorithm this is the residual cutoff.
        :type minchange: float
        :param chroms: A list of chromosomes to calculate corrections for. If set as None, all chromosome corrections are found.
        :type chroms: list
        :param precorrect: Use binning-based corrections in expected value calculations, resulting in a chained normalization approach.
        :type precorrect: bool.
        :param binary: Use binary indicator instead of counts.
        :type binary: bool.
        :param kr: Use the Knight Ruiz matrix balancing algorithm instead of weighted matrix balancing. This option ignores 'iterations'.
        :type kr: bool.
        :returns: None

        :Attributes: * **corrections** (*ndarray*) - A numpy array of type float32 and length equal to the number of fends. All invalid fends have an associated correction value of zero.

        The 'normalization' attribute is updated to 'express' or 'binning-express', depending on if the 'precorrect' option is selected. In addition, the 'chromosome_means' attribute is updated such that the mean correction (sum of all valid chromosomal correction value pairs) is adjusted to zero and the corresponding chromosome mean is adjusted the same amount but the opposite sign. 
        """
        self.history += "HiC.find_express_fend_corrections(iterations=%i, mindistance=%i, maxdistance=%s, remove_distance=%s, usereads='%s', mininteractions=%i, chroms=%s, precorrect=%s, binary=%s, kr=%s) - " % (iterations, mindistance, str(maxdistance), remove_distance, usereads, mininteractions, str(chroms), precorrect, binary, kr)
        if mininteractions is None:
            if 'mininteractions' in self.__dict__.keys():
                mininteractions = self.mininteractions
            else:
                mininteractions = 1
        # make sure usereads has a valid value
        if usereads not in ['cis', 'trans', 'all']:
            if not self.silent:
                print >> sys.stderr, ("'usereads' does not have a valid value.\n"),
            self.history += "Error: '%s' not a valid value for 'usereads'\n" % usereads
            return None
        if usereads == 'cis':
            useread_int = 0
        elif usereads == 'all':
            useread_int = 1
        else:
            useread_int = 2
        # make sure distance parameters have been calculated, if removing distance
        if remove_distance and self.distance_parameters is None:
            if not self.silent:
                print >> sys.stderr, ("'remove_distance' requires a project that has already had 'find_distance_parameters' run.\n"),
            self.history += "Error: 'find_distance_parameters()' not run yet\n"
            return None
        if precorrect and self.binning_corrections is None:
            if not self.silent:
                print >> sys.stderr, ("Precorrection can only be used in project has previously run 'find_binning_fend_corrections'.\n"),
            self.history += "Error: 'find_binning_fend_corrections()' not run yet\n"
            return None
        if self.corrections is None:
            self.corrections = numpy.ones(self.fends['fends'].shape[0], dtype=numpy.float32)
        if kr:
            self._find_kr_corrections(mindistance, maxdistance, remove_distance, 
                                      usereads, mininteractions, minchange, chroms, precorrect,
                                      binary)
            return None
        # create needed arrays
        mids = self.fends['fends']['mid'][...]
        filt = numpy.copy(self.filter)
        chr_indices = self.fends['chr_indices'][...]
        if (chroms is None or
                (isinstance(chroms, list) and
                (len(chroms) == 0 or
                (len(chroms) == 1 and chroms[0] == ''))) or
                chroms == ''):
            chroms = self.chr2int.keys()
            chroms.sort()
        elif not chroms is None and isinstance(chroms, str):
            chroms = [chroms]
        chrints = numpy.zeros(len(chroms), dtype=numpy.int32)
        for i in range(len(chroms)):
            chrints[i] = self.chr2int[chroms[i]]
        chroms = list(numpy.array(chroms)[numpy.argsort(chrints)])
        chrints = chrints[numpy.argsort(chrints)]
        for chrm, i in self.chr2int.iteritems():
            if chrm not in chroms:
                filt[chr_indices[i]:chr_indices[i + 1]] = 0
        # copy needed arrays from h5dict
        if not self.silent:
            print >> sys.stderr, ("\r%s\rLoading needed data...") % (' ' * 80),
        if usereads in ['cis', 'all']:
            cis_ranges = numpy.round(numpy.linspace(0, self.data['cis_data'].shape[0],
                                                    self.num_procs + 1)).astype(numpy.int64)
            data = self.data['cis_data'][cis_ranges[self.rank]:cis_ranges[self.rank + 1], :]
            distances = mids[data[:, 1]] - mids[data[:, 0]]
            if maxdistance == 0 or maxdistance is None:
                maxdistance = numpy.amax(distances) + 1
                if self.rank == 0:
                    for i in range(1, self.num_procs):
                        maxdistance = max(maxdistance, self.comm.recv(source=i, tag=11))
                    for i in range(1, self.num_procs):
                        self.comm.send(maxdistance, dest=i, tag=11)
                else:
                    self.comm.send(maxdistance, dest=0, tag=11)
                    maxdistance = self.comm.recv(source=0, tag=11)
            valid = numpy.where(filt[data[:, 0]] * filt[data[:, 1]] *
                                (distances >= mindistance) * (distances < maxdistance))[0]
            data = data[valid, :]
        else:
            data = None
        if usereads in ['trans', 'all']:
            trans_ranges = numpy.round(numpy.linspace(0, self.data['trans_data'].shape[0],
                                                      self.num_procs + 1)).astype(numpy.int64)
            trans_data = self.data['trans_data'][trans_ranges[self.rank]:trans_ranges[self.rank + 1], :]
            valid = numpy.where(filt[trans_data[:, 0]] * filt[trans_data[:, 1]])[0]
            trans_data = trans_data[valid, :]
        else:
            trans_data = None
        trans_means = None
        if not self.silent:
            print >> sys.stderr, ("\r%s\rChecking for fend interaction count...") % (' ' * 80),
        # double check that, given the type of reads being used for learning, there are enough for each fend
        # to meet the mininteraction criteria
        observed_interactions = numpy.zeros(filt.shape[0], dtype=numpy.int32)
        if not trans_data is None:
            observed_interactions += numpy.bincount(trans_data[:, 0], minlength=filt.shape[0])
            observed_interactions += numpy.bincount(trans_data[:, 1], minlength=filt.shape[0])
        if not data is None:
            observed_interactions += numpy.bincount(data[:, 0], minlength=filt.shape[0])
            observed_interactions += numpy.bincount(data[:, 1], minlength=filt.shape[0])
        if self.rank == 0:
            for i in range(1, self.num_procs):
                observed_interactions += self.comm.recv(source=i, tag=11)
            minobs = numpy.amin(observed_interactions[numpy.where(filt)])
            for i in range(1, self.num_procs):
                self.comm.send(minobs, dest=i, tag=11)
        else:
            self.comm.send(observed_interactions, dest=0, tag=11)
            minobs = self.comm.recv(source=0, tag=11)
        if minobs < mininteractions:
            if not self.silent:
                print >> sys.stderr, ("\nInsufficient interactions for one or more fends.\n"),
                print >> sys.stderr, ("Try resetting and refiltering fends or expanding distance range.\n"),
            self.history += "Error: Too few interactions for give settings\n"
            return None
        if self.rank == 0:
            chrints = [self.chr2int[chrom] for chrom in chroms]
            chrints = numpy.array(chrints, dtype=numpy.int32)
            numpy.random.shuffle(chrints)
            chrint_ranges = numpy.round(numpy.linspace(0, chrints.shape[0], self.num_procs + 1)).astype(numpy.int32)
            for i in range(1, self.num_procs):
                self.comm.send(chrints[chrint_ranges[i]:chrint_ranges[i + 1]], dest=i, tag=11)
            chrints = chrints[:chrint_ranges[1]]
        else:
            chrints = self.comm.recv(source=0, tag=11)
        interactions = numpy.zeros(filt.shape[0], dtype=numpy.int64)
        _interactions.find_distancebound_possible_interactions(interactions,
                                                               chr_indices,
                                                               mids,
                                                               filt,
                                                               chrints,
                                                               useread_int,
                                                               mindistance,
                                                               maxdistance)
        if self.rank == 0:
            for i in range(1, self.num_procs):
                interactions += self.comm.recv(source=i, tag=11)
            for i in range(1, self.num_procs):
                self.comm.Send(interactions, dest=i, tag=13)
        else:
            self.comm.send(interactions, dest=0, tag=11)
            self.comm.Recv(interactions, source=0, tag=13)
        # precalculate interaction distance means for all included interactions
        if not remove_distance or data is None:
            distance_means = None
            if trans_data is None:
                if binary:
                    datasum = data.shape[0]
                else:
                    datasum = numpy.sum(data[:, 2])
            elif data is None:
                if binary:
                    datasum = trans_data.shape[0]
                else:
                    datasum = numpy.sum(trans_data[:, 2])
            else:
                if binary:
                    datasum = data.shape[0] + trans_data.shape[0]
                else:
                    datasum = numpy.sum(data[:, 2]) + numpy.sum(trans_data[:, 2])
            if self.rank == 0:
                for i in range(1, self.num_procs):
                    datasum += self.comm.recv(source=i, tag=11)
                temp_mu = (2.0 * datasum) / numpy.sum(interactions)
                for i in range(1, self.num_procs):
                    self.comm.send(temp_mu, dest=i, tag=11)
            else:
                self.comm.send(datasum, dest=0, tag=11)
                temp_mu = self.comm.recv(source=0, tag=11)
            if trans_data is None:
                trans_mu = 1.0
                mu = temp_mu
            elif data is None:
                mu = 1.0
                trans_mu = temp_mu
            else:
                mu = temp_mu
                trans_mu = temp_mu
        else:
            if not self.silent:
                print >> sys.stderr, ("\r%s\rPrecalculating distances...") % (' ' * 80),
            mu = 1.0
            all_chrints = self.fends['fends']['chr'][data[:, 0]]
            chrint_indices = numpy.r_[0, numpy.bincount(all_chrints)]
            for i in range(1, chrint_indices.shape[0]):
                chrint_indices[i] += chrint_indices[i - 1]
            distance_means = numpy.zeros(data.shape[0], dtype=numpy.float32)
            for i in range(chrint_indices.shape[0] - 1):
                if chrint_indices[i] < chrint_indices[i + 1]:
                    if binary:
                        _distance.find_remapped_distance_means(data[chrint_indices[i]:chrint_indices[i + 1], 0],
                                                               data[chrint_indices[i]:chrint_indices[i + 1], 1],
                                                               mids,
                                                               distance_means[chrint_indices[i]:chrint_indices[i + 1]],
                                                               self.bin_distance_parameters,
                                                               self.chromosome_means[i])
                    else:
                        _distance.find_remapped_distance_means(data[chrint_indices[i]:chrint_indices[i + 1], 0],
                                                               data[chrint_indices[i]:chrint_indices[i + 1], 1],
                                                               mids,
                                                               distance_means[chrint_indices[i]:chrint_indices[i + 1]],
                                                               self.distance_parameters,
                                                               self.chromosome_means[i])
            if not trans_data is None:
                total_possible = numpy.sum(filt).astype(numpy.int64) ** 2
                for i in range(self.fends['chr_indices'].shape[0] - 1):
                    start = self.fends['chr_indices'][i]
                    stop = self.fends['chr_indices'][i + 1]
                    total_possible -= numpy.sum(filt[start:stop].astype(numpy.int64)) ** 2
                if self.rank == 0:
                    if binary:
                        trans_sum = trans_data.shape[0]
                    else:
                        trans_sum = numpy.sum(trans_data[:, 2])
                    for i in range(1, self.num_procs):
                        trans_sum += self.comm.recv(source=i, tag=11)
                    for i in range(1, self.num_procs):
                        self.comm.send(trans_sum, dest=i, tag=11)
                else:
                    self.comm.send(numpy.sum(trans_data[:, 2]), dest=0, tag=11)
                    trans_sum = self.comm.recv(source=0, tag=11)
                trans_mu = (2.0 * trans_sum) / total_possible
            else:
                trans_mu = 1.0
        if precorrect:
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding binning corrections...") % (' ' * 80),
            if not data is None:
                if distance_means is None:
                    distance_means = numpy.ones(data.shape[0], dtype=numpy.float32)
                _optimize.find_binning_correction_adjustment(distance_means,
                                                             data[:, 0],
                                                             data[:, 1],
                                                             self.binning_corrections,
                                                             self.binning_num_bins,
                                                             self.binning_fend_indices)
            if not trans_data is None:
                trans_means = numpy.ones(trans_data.shape[0], dtype=numpy.float32)
                _optimize.find_binning_correction_adjustment(trans_means,
                                                             trans_data[:, 0],
                                                             trans_data[:, 1],
                                                             self.binning_corrections,
                                                             self.binning_num_bins,
                                                             self.binning_fend_indices)
        if not self.silent and self.rank == 0:
            print >> sys.stderr, ("\r%s\rFinding fend corrections...") % (' ' * 80),
        # calculate corrections
        fend_means = numpy.zeros(filt.shape[0], dtype=numpy.float64)
        if self.rank == 0:
            temp = numpy.zeros(filt.shape[0], dtype=numpy.float64)
        corrections = numpy.copy(self.corrections)
        cont = True
        iteration = 0
        change = numpy.zeros(1, dtype=numpy.float64)
        while cont:
            iteration += 1
            _optimize.find_fend_means(distance_means,
                                      trans_means,
                                      fend_means,
                                      data,
                                      trans_data,
                                      corrections,
                                      mu,
                                      trans_mu,
                                      int(binary))
            if self.rank == 0:
                for i in range(1, self.num_procs):
                    self.comm.Recv(temp, source=i, tag=13)
                    fend_means += temp
                cost = _optimize.update_express_corrections(filt,
                                                            interactions,
                                                            fend_means,
                                                            corrections,
                                                            change)
                if iteration >= iterations or change[0] < minchange:
                    cont = False
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rFinding fend corrections  Iteration: %i  Cost: %f  Change: %f") % (' ' * 80,
                                          iteration, cost, change),
                for i in range(1, self.num_procs):
                    self.comm.Send(corrections, dest=i, tag=13)
                    self.comm.send(cont, dest=i, tag=11)
            else:
                self.comm.Send(fend_means, dest=0, tag=13)
                self.comm.Recv(corrections, source=0, tag=13)
                cont = self.comm.recv(source=0, tag=11)
        # calculate chromosome mean
        if self.chromosome_means is None:
            self.chromosome_means = numpy.zeros(self.fends['chr_indices'].shape[0] - 1, dtype=numpy.float32)
        for chrom in chroms:
            chrint = self.chr2int[chrom]
            valid = numpy.where(self.filter[chr_indices[chrint]:chr_indices[chrint + 1]])[0] + chr_indices[chrint]
            if valid.shape[0] == 0:
                continue
            chrom_mean = numpy.sum(corrections[valid])
            chrom_mean = chrom_mean ** 2.0 - numpy.sum(corrections[valid] ** 2.0)
            chrom_mean /= valid.shape[0] * (valid.shape[0] - 1)
            if remove_distance:
                self.chromosome_means[chrint] += numpy.log(chrom_mean)
            corrections[valid] /= chrom_mean ** 0.5
        self.corrections = corrections
        if not self.silent and self.rank == 0:
            print >> sys.stderr, ("\r%s\rCompleted learning express corrections. Final cost: %f\n") % (' ' * 80, cost),
        if precorrect:
            self.normalization = 'binning-express'
        else:
            self.normalization = 'express'
        self.history += "Succcess\n"
        return None

    def _find_kr_corrections(self, mindistance=0, maxdistance=0, remove_distance=True, 
                             usereads='cis', mininteractions=0, minchange=0.0001, chroms=[], precorrect=False,
                             binary=False):
        if (chroms is None or
                (isinstance(chroms, list) and
                (len(chroms) == 0 or
                (len(chroms) == 1 and chroms[0] == ''))) or
                chroms == ''):
            chroms = self.chr2int.keys()
            chroms.sort()
        all_chroms = list(chroms)
        filt = numpy.copy(self.filter)
        if maxdistance == 0 or maxdistance is None:
            maxdistance = 99999999999
        if usereads != 'cis':
            for chrom, i in self.chr2int.iteritems():
                if chrom not in chroms:
                    filt[self.fends['chr_indices'][i]:self.fends['chr_indices'][i + 1]] = 0
            chroms = ['all']
        for chrom in chroms:
            if chrom == 'all':
                startfend = 0
                stopfend = self.fends['chr_indices'][-1]
                chrfilt = filt
            else:
                chrint = self.chr2int[chrom]
                startfend = self.fends['chr_indices'][chrint]
                stopfend = self.fends['chr_indices'][chrint + 1]
                chrfilt = filt[startfend:stopfend]
            # create needed arrays
            if not self.silent:
                print >> sys.stderr, ("\r%s\rLoading needed data...") % (' ' * 80),
            mids = self.fends['fends']['mid'][startfend:stopfend]
            if usereads in ['cis', 'all']:
                start_index = self.data['cis_indices'][startfend]
                stop_index = self.data['cis_indices'][stopfend]
                cis_ranges = numpy.round(numpy.linspace(start_index, stop_index,
                                                        self.num_procs + 1)).astype(numpy.int64)
                data = self.data['cis_data'][cis_ranges[self.rank]:cis_ranges[self.rank + 1], :]
                distances = mids[data[:, 1] - startfend] - mids[data[:, 0] - startfend]
                valid = numpy.where(chrfilt[data[:, 0] - startfend] * chrfilt[data[:, 1] - startfend] *
                                    (distances >= mindistance) * (distances < maxdistance))[0]
                data = data[valid, :]
            else:
                data = None
            if usereads in ['trans', 'all']:
                trans_ranges = numpy.round(numpy.linspace(0, self.data['trans_data'].shape[0],
                                                          self.num_procs + 1)).astype(numpy.int64)
                trans_data = self.data['trans_data'][trans_ranges[self.rank]:trans_ranges[self.rank + 1], :]
                valid = numpy.where(filt[trans_data[:, 0]] * filt[trans_data[:, 1]])[0]
                trans_data = trans_data[valid, :]
            else:
                trans_data = None
            trans_means = None
            # remapped data
            rev_mapping = numpy.where(chrfilt)[0]
            mapping = numpy.zeros(chrfilt.shape[0], dtype=numpy.int32) - 1
            mapping[rev_mapping] = numpy.arange(rev_mapping.shape[0])
            if not data is None:
                data[:, 0] = mapping[data[:, 0] - startfend]
                data[:, 1] = mapping[data[:, 1] - startfend]
            if not trans_data is None:
                trans_data[:, 0] = mapping[trans_data[:, 0]]
                trans_data[:, 1] = mapping[trans_data[:, 1]]
            mids = mids[rev_mapping]
            if not self.silent:
                print >> sys.stderr, ("\r%s\rChecking for fend interaction count...") % (' ' * 80),
            # double check that, given the type of reads being used for learning, there are enough for each fend
            # to meet the mininteraction criteria
            observed_interactions = numpy.zeros(rev_mapping.shape[0], dtype=numpy.int32)
            if not trans_data is None:
                observed_interactions += numpy.bincount(trans_data[:, 0], minlength=rev_mapping.shape[0])
                observed_interactions += numpy.bincount(trans_data[:, 1], minlength=rev_mapping.shape[0])
            if not data is None:
                observed_interactions += numpy.bincount(data[:, 0], minlength=rev_mapping.shape[0])
                observed_interactions += numpy.bincount(data[:, 1], minlength=rev_mapping.shape[0])
            if self.rank == 0:
                for i in range(1, self.num_procs):
                    observed_interactions += self.comm.recv(source=i, tag=11)
                minobs = numpy.amin(observed_interactions)
                for i in range(1, self.num_procs):
                    self.comm.send(minobs, dest=i, tag=11)
            else:
                self.comm.send(observed_interactions, dest=0, tag=11)
                minobs = self.comm.recv(source=0, tag=11)
            if minobs < mininteractions:
                if not self.silent:
                    print >> sys.stderr, ("\nInsufficient interactions for one or more fends.\n"),
                    print >> sys.stderr, ("Try resetting and refiltering fends or expanding distance range.\n"),
                self.history += "Error: Too few interactions for give settings\n"
                return None
            # precalculate interaction distance means for all included interactions
            if not data is None:
                if binary:
                    counts = numpy.ones(data.shape[0], dtype=numpy.float64)
                else:
                    counts = data[:, 2].astype(numpy.float64)
            else:
                counts = None
            if not trans_data is None:
                if binary:
                    trans_counts = numpy.ones(trans_data.shape[0], dtype=numpy.float64)
                else:
                    trans_counts = trans_data[:, 2].astype(numpy.float64)
            else:
                trans_counts = None
            trans_means = None
            distance_means = None
            if remove_distance:
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rPrecalculating distances...") % (' ' * 80),
                if usereads != 'cis':
                    trans_mean = numpy.sum(trans_counts).astype(numpy.float64)
                    interactions = rev_mapping.shape[0] ** 2
                    all_chrints = self.fends['fends']['chr'][rev_mapping]
                    interactions -= numpy.sum(numpy.bincount(all_chrints,
                                     minlength=self.fends['chr_indices'].shape[0]) ** 2)
                    trans_mean /= interactions
                    trans_means = numpy.empty(trans_data.shape[0], dtype=numpy.float32).fill(trans_mean)
                    if not data is None:
                        distance_means = numpy.zeros(data.shape[0], dtype=numpy.float32)
                        indices = numpy.r_[0, numpy.bincount(all_chrints)]
                        for i in range(1, indices.shape[0]):
                            indices[i] += indices[i - 1]
                        for i in range(indices.shape[0] - 1):
                            if indices[i] < indices[i + 1]:
                                if binary:
                                    _distance.find_remapped_distance_means(data[indices[i]:indices[i + 1], 0],
                                                                           data[indices[i]:indices[i + 1], 1],
                                                                           mids,
                                                                           distance_means[indices[i]:indices[i + 1]],
                                                                           self.bin_distance_parameters,
                                                                           self.chromosome_means[i])
                                else:
                                    _distance.find_remapped_distance_means(data[indices[i]:indices[i + 1], 0],
                                                                           data[indices[i]:indices[i + 1], 1],
                                                                           mids,
                                                                           distance_means[indices[i]:indices[i + 1]],
                                                                           self.distance_parameters,
                                                                           self.chromosome_means[i])
                else:
                    distance_means = numpy.zeros(data.shape[0], dtype=numpy.float32)
                    chrint = self.chr2int[chrom]
                    if binary:
                        _distance.find_remapped_distance_means(data[:, 0],
                                                               data[:, 1],
                                                               mids,
                                                               distance_means,
                                                               self.bin_distance_parameters,
                                                               self.chromosome_means[chrint])
                    else:
                        _distance.find_remapped_distance_means(data[:, 0],
                                                               data[:, 1],
                                                               mids,
                                                               distance_means,
                                                               self.distance_parameters,
                                                               self.chromosome_means[chrint])
            if precorrect:
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rFinding binning corrections...") % (' ' * 80),
                if not data is None:
                    if distance_means is None:
                        distance_means = numpy.ones(data.shape[0], dtype=numpy.float32)
                    _optimize.find_binning_correction_adjustment(distance_means,
                                                                 data[:, 0],
                                                                 data[:, 1],
                                                                 self.binning_corrections,
                                                                 self.binning_num_bins,
                                                                 self.binning_fend_indices[rev_mapping, :, :])
                if not trans_data is None:
                    if trans_means is None:
                        trans_means = numpy.ones(trans_data.shape[0], dtype=numpy.float32)
                    _optimize.find_binning_correction_adjustment(trans_means,
                                                                 trans_data[:, 0],
                                                                 trans_data[:, 1],
                                                                 self.binning_corrections,
                                                                 self.binning_num_bins,
                                                                 self.binning_fend_indices[rev_mapping, :, :])
            if not distance_means is None:
                counts /= distance_means
            if not trans_means is None:
                trans_counts /= trans_means
            if not self.silent and self.rank == 0:
                print >> sys.stderr, ("\r%s\rFinding fend corrections...") % (' ' * 80),
            # calculate corrections
            if self.rank == 0:
                temp = numpy.zeros((rev_mapping.shape[0], 1), dtype=numpy.float64)
            corrections = numpy.ones((rev_mapping.shape[0], 1), dtype=numpy.float64)
            g = 0.9
            eta = etamax = 0.1
            stop_tol = minchange * 0.5
            rt = minchange ** 2.0
            delta = 0.1
            Delta = 3
            v = numpy.zeros((corrections.shape[0], 1), dtype=numpy.float64)
            w = numpy.zeros((corrections.shape[0], 1), dtype=numpy.float64)
            _optimize.calculate_v(data, trans_data, counts, trans_counts, corrections, v)
            if self.rank == 0:
                for i in range(1, self.num_procs):
                    self.Recv(temp, source=i, tag=13)
                    v += temp
                for i in range(1, self.num_procs):
                    self.Send(v, dest=i, tag=13)
            else:
                self.Send(v, dest=0, tag=13)
                self.Recv(v, source=0, tag=13)
            rk = 1.0 - v
            rho_km1 = numpy.dot(rk.T, rk)[0, 0]
            rho_km2 = rho_km1
            rold = rout = rho_km1
            i = MVP = 0
            while rout > rt:
                i += 1
                k = 0
                y = numpy.ones((rev_mapping.shape[0], 1), dtype=numpy.float64)
                innertol = max(eta ** 2.0 * rout, rt)
                while rho_km1 > innertol:
                    k += 1
                    if k == 1:
                        Z = rk / v
                        p = numpy.copy(Z)
                        rho_km1 = numpy.dot(rk.T, Z)
                    else:
                        beta = rho_km1 / rho_km2
                        p = Z + beta * p
                    # Update search direction efficiently
                    w.fill(0.0)
                    _optimize.calculate_w(data, trans_data, counts, trans_counts, corrections, p, w)
                    if self.rank == 0:
                        for j in range(1, self.num_procs):
                            self.Recv(temp, source=j, tag=13)
                            w += temp
                        for j in range(1, self.num_procs):
                            self.Send(w, dest=j, tag=13)
                    else:
                        self.Send(w, dest=0, tag=13)
                        self.Recv(w, source=0, tag=13)
                    w += v * p
                    alpha = rho_km1 / numpy.dot(p.T, w)[0, 0]
                    ap = alpha * p
                    # Test distance to boundary of cone
                    ynew = y + ap
                    if numpy.amin(ynew) <= delta:
                        if delta == 0:
                            break
                        ind = numpy.where(ap < 0.0)[0]
                        gamma = numpy.amin((delta - y[ind]) / ap[ind])
                        y += gamma * ap
                        break
                    if numpy.amax(ynew) >= Delta:
                        ind = numpy.where(ynew > Delta)[0]
                        gamma = numpy.amin((Delta - y[ind]) / ap[ind])
                        y += gamma * ap
                        break
                    y = numpy.copy(ynew)
                    rk -= alpha * w
                    rho_km2 = rho_km1
                    Z = rk / v
                    rho_km1 = numpy.dot(rk.T, Z)[0, 0]
                corrections *= y
                v.fill(0.0)
                _optimize.calculate_v(data, trans_data, counts, trans_counts, corrections, v)
                if self.rank == 0:
                    for j in range(1, self.num_procs):
                        self.Recv(temp, source=j, tag=13)
                        v += temp
                    for j in range(1, self.num_procs):
                        self.Send(v, dest=j, tag=13)
                else:
                    self.Send(v, dest=0, tag=13)
                    self.Recv(v, source=0, tag=13)
                rk = 1.0 - v
                rho_km1 = numpy.dot(rk.T, rk)[0, 0]
                rout = rho_km1
                MVP += k + 1
                # Update inner iteration stopping criterion
                rat = rout / rold
                rold = rout
                res_norm = rout ** 0.5
                eta_o = eta
                eta = g * rat
                if g * eta_o ** 2.0 > 0.1:
                    eta = max(eta, g * eta_o ** 2.0)
                eta = max(min(eta, etamax), stop_tol / res_norm)
                if not self.silent and self.rank == 0:
                    print >> sys.stderr, ("\r%s\rIteration %i Residual: %e") % (" " * 80, i, rout),
            if not self.silent and self.rank == 0:
                print >> sys.stderr, ("\r%s\rFinding fend corrections... Chrom: %s Done\n") % (' ' * 80, chrom),
            self.corrections[rev_mapping + startfend] = 1.0 / corrections
        # calculate chromosome mean
        if self.chromosome_means is None:
            self.chromosome_means = numpy.zeros(self.fends['chr_indices'].shape[0] - 1, dtype=numpy.float32)
        chr_indices = self.fends['chr_indices'][...]
        for chrom in all_chroms:
            chrint = self.chr2int[chrom]
            valid = numpy.where(self.filter[chr_indices[chrint]:chr_indices[chrint + 1]])[0] + chr_indices[chrint]
            if valid.shape[0] == 0:
                continue
            chrom_mean = numpy.sum(self.corrections[valid])
            chrom_mean = chrom_mean ** 2.0 - numpy.sum(self.corrections[valid] ** 2.0)
            chrom_mean /= valid.shape[0] * (valid.shape[0] - 1)
            if remove_distance:
                self.chromosome_means[chrint] += numpy.log(chrom_mean)
            self.corrections[valid] /= chrom_mean ** 0.5
        if not self.silent and self.rank == 0:
            print >> sys.stderr, ("\r%s\rCompleted learning express corrections.\n") % (' ' * 80),
        if precorrect:
            self.normalization = 'binning-express'
        else:
            self.normalization = 'express'
        self.history += "Succcess\n"
        return None

    def find_binning_fend_corrections(self, mindistance=0, maxdistance=0, chroms=[], num_bins=[20, 20, 20],
                                      parameters=['even', 'even', 'even-const'], model=['gc', 'len', 'distance'],
                                      learning_threshold=1.0, max_iterations=10, usereads='cis', pseudocounts=0):
        """
        Using a multivariate binning model, learn correction values for combinations of model parameter bins. This function is MPI compatible.

        :param mindistance: The minimum inter-fend distance to be included in modeling.
        :type mindistance: int.
        :param maxdistance: The maximum inter-fend distance to be included in modeling.
        :type maxdistance: int.
        :param chroms: A list of chromosomes to calculate corrections for. If set as None, all chromosome corrections are found.
        :type chroms: list
        :param remove_distance: Use distance dependence curve in prior probability calculation for each observation.
        :type remove_distance: bool.
        :param model: A list of fend features to be used in model. Valid values are 'len', 'distance', and any features included in the creation of the associated Fend object. The 'distance' parameter is only good with 'cis' or 'all' reads. If used with 'all', distances will be partitioned into n - 1 bins and the final distance bin will contain all trans data.
        :type model: list
        :param num_bins: A list of the number of approximately equal-sized bins two divide model components into.
        :type num_bins: list
        :param parameters: A list of types, one for each model parameter. Types can be either 'even' or 'fixed', indicating whether each parameter bin should contain approximately even numbers of interactions or be of fixed width spanning 1 / Nth of the range of the parameter's values, respectively. Parameter types can also have the suffix '-const' to indicate that the parameter should not be optimized.
        :type parameters: list
        :param learning_threshold: The minimum change in log-likelihood needed to continue iterative learning process.
        :type learning_threshold: float
        :param max_iterations: The maximum number of iterations to use for learning model parameters.
        :type max_iterations: int.
        :param usereads: Specifies which set of interactions to use, 'cis', 'trans', and 'all'.
        :type usereads: str.
        :param pseudocounts: The number of pseudo-counts to add to each bin prior to seeding and learning normalization values.
        :type pseudocounts: int.
        :returns: None
 
        :Attributes: * **model_parameters** (*ndarray*) - A numpy array of strings containing model parameter names. If distance was included in the 'model' option, it is not included in this array since it is only for learning values, not for subsequent corretion.
                     * **binning_num_bins** (*ndarray*) - A numpy array of type int32 containing the number of bins for each non-distance model parameter.
                     * **binning corrections** (*ndarray*) - A numpy array of type float32 and length equal to the sum of binning_num_bins * (binning_num_bins - 1) / 2. This array contains a 1D stack of correction values, ordered according to the parameter order in the 'model_parameters' attribute.
                     * **binning_correction_indices** (*ndarray*) - A numpy array of type int32 and length equal to the number of non-distance model parameters plus one. This array contains the first position in 'binning_corrections' for the first bin of the model parameter in the corresponding position in the 'model_parameters' array. The last position in the array contains the total number of binning correction values.
                     * **binning_fend_indices** (*ndarray*) - A numpy array of type int32 and size N x M x 2 where M is the number of non-distance model parameters and N is the number of fends. This array contains the binning index for each parameter for each fend for the first and the second position in the correction array.
        
        The 'normalization' attribute is updated to 'binning'.
        """
        self.history += "HiC.find_binning_fend_corrections(max_iterations=%i, mindistance=%i, maxdistance=%s, usereads='%s', chroms=%s num_bins=%s, parameters=%s, model=%s, learning_threshold=%f, pseudocounts=%s) - " % (max_iterations, mindistance, str(maxdistance), usereads, str(chroms), num_bins, parameters, model, learning_threshold, str(pseudocounts))
        for parameter in model:
            if not parameter in ['len', 'distance'] and parameter not in self.fends['fends'].dtype.names:
                if not self.silent:
                    print >> sys.stderr, ("Fend feature %s not found in fend object. Try removing it from model or creating a new fend object with feature data.\n") % (parameter),
                self.history += "Error: model parameter '%s' not found in fend data\n" % parameter
                return None
        for parameter in parameters:
            if parameter not in ['even', 'fixed', 'even-const', 'fixed-const']:
                if not self.silent:
                    print >> sys.stderr, ("Fend feature type %s is not valid.") % (parameter),
                self.history += "Error: model feature type '%s' not valid\n" % parameter
                return None
        if 'distance' in model and usereads == 'trans':
            if not self.silent:
                print >> sys.stderr, ("The 'distance' parameter can only be used in conjuction with 'cis' or 'all' reads.\n"),
            self.history += "Error: model parameter 'distance' cannot be used with 'trans'\n"
            return None
        if len(model) != len(num_bins):
            if not self.silent:
                print >> sys.stderr, ("The number of items in the 'model' parameter must be the same as the number in the 'num_bins' parameter.\n"),
            self.history += "Error: mismatch between lengths of 'num_bins' and 'model'\n"
            return None
        # make sure usereads has a valid value
        if usereads not in ['cis', 'trans', 'all']:
            if not self.silent:
                print >> sys.stderr, ("'usereads' does not have a valid value.\n"),
            self.history += "Error: '%s' not a valid value for 'usereads'\n" % usereads
            return None
        if (chroms is None or
                (isinstance(chroms, list) and
                (len(chroms) == 0 or
                (len(chroms) == 1 and chroms[0] == ''))) or
                chroms == ''):
            chroms = self.chr2int.keys()
        chrints = numpy.zeros(len(chroms), dtype=numpy.int32)
        for i in range(len(chroms)):
            chrints[i] = self.chr2int[chroms[i]]
        chroms = list(numpy.array(chroms)[numpy.argsort(chrints)])
        # Create requested bin cutoffs for each feature of the binning model
        if not self.silent:
            print >> sys.stderr, ("\r%s\rPartitioning features into bins...") % (' ' * 80),
        filt = numpy.copy(self.filter)
        chr_indices = self.fends['chr_indices'][...]
        for chrom, i in self.chr2int.iteritems():
            if chrom not in chroms:
                filt[chr_indices[i]:chr_indices[i + 1]] = 0
        if maxdistance == 0 or maxdistance is None:
            for chrom in chroms:
                chrint = self.chr2int[chrom]
                maxdistance = max(maxdistance, self.fends['fends']['mid'][chr_indices[chrint + 1] - 1] -
                                               self.fends['fends']['mid'][chr_indices[chrint]] + 1)
        valid = numpy.where(filt == 1)[0]
        if not self.silent:
            print >> sys.stderr, ("\r%s\rPartitioning features into bins...") % (' ' * 80),
        if 'distance' in model:
            use_distance = True
            index = model.index('distance')
            model = model[:index] + model[(index + 1):]
            distance_bins = num_bins[index]
            num_bins = num_bins[:index] + num_bins[(index + 1):]
            parameters = parameters[:index] + parameters[(index + 1):]
        else:
            use_distance = False
        num_bins = numpy.array(num_bins, dtype=numpy.int32)
        total_bins = 1
        all_bins = numpy.zeros(0, dtype=numpy.float32)
        all_corrections = numpy.ones(0, dtype=numpy.float64)
        all_indices = numpy.zeros((filt.shape[0], len(model), 2), dtype=numpy.int32)
        bin_indices = numpy.zeros(len(model) + 1, dtype=numpy.int32)
        correction_indices = numpy.zeros(len(model) + 1, dtype=numpy.int32)
        bin_divs = numpy.zeros(len(model), dtype=numpy.int32)
        for i in range(len(model)):
            if model[i] == 'len':
                values = (self.fends['fends']['stop'][...] -
                          self.fends['fends']['start'][...]).astype(numpy.float32)
            else:
                values = self.fends['fends'][model[i]][...]
            if parameters[i].count('even') > 0:
                temp = numpy.copy(values)
                temp.sort()
                all_bins = numpy.hstack((all_bins, temp[numpy.round(numpy.linspace(0, values.shape[0],
                                        num_bins[i] + 1)).astype(numpy.int32)[1:] - 1])).astype(numpy.float32)
            else:
                all_bins = numpy.hstack((all_bins, numpy.linspace(numpy.amin(values),
                                        numpy.amax(values), num_bins[i] + 1)[1:])).astype(numpy.float32)
            all_bins[-1] = numpy.inf
            bin_indices[i + 1] = all_bins.shape[0]
            all_corrections = numpy.hstack((all_corrections, numpy.zeros(num_bins[i] * (num_bins[i] + 1) / 2,
                                           dtype=numpy.float64)))
            correction_indices[i + 1] = all_corrections.shape[0]
            bin_divs[i] = total_bins
            total_bins *= num_bins[i] * (num_bins[i] + 1) / 2
            all_indices[:, i, 0] = numpy.searchsorted(all_bins[bin_indices[i]:bin_indices[i + 1]],
                                                      values).astype(numpy.int32)
            all_indices[:, i, 1] = (all_indices[:, i, 0] * (num_bins[i] - 1) - all_indices[:, i, 0] *
                                    (all_indices[:, i, 0] - 1) / 2)
        self.binning_num_bins = num_bins
        mids = self.fends['fends']['mid'][...]
        if use_distance:
            if maxdistance == 0:
                max_dist = 0
                dists = mids[chr_indices[1:] - 1] - mids[chr_indices[:-1]]
                for chrom in chroms:
                    max_dist = max(max_dist, dists[self.chr2int[chrom]])
            else:
                max_dist = maxdistance
            if usereads == 'all':
                distance_cutoffs = numpy.ceil(numpy.exp(numpy.linspace(numpy.log(max(1, mindistance)),
                                              numpy.log(max_dist), distance_bins))[1:]).astype(numpy.int32)
            else:
                distance_cutoffs = numpy.ceil(numpy.exp(numpy.linspace(numpy.log(max(1, mindistance)),
                                              numpy.log(max_dist), distance_bins + 1))[1:]).astype(numpy.int32)
            distance_cutoffs[-1] += 1
            distance_corrections = numpy.zeros(distance_bins, dtype=numpy.float64)
            distance_div = total_bins
            total_bins *= distance_bins
        else:
            distance_cutoffs = None
            distance_corrections = None
            distance_div = 0
            distance_bins = 0
        bin_counts = numpy.zeros((total_bins, 2), dtype=numpy.int64)
        temp_counts = None
        if self.rank == 0:
            temp_counts = numpy.copy(bin_counts)
        for h, chrom in enumerate(chroms):
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding bin counts... chr%s") % (' ' * 80, chrom),
            chrint = self.chr2int[chrom]
            if usereads in ['cis', 'all']:
                # Find number of observations in each bin
                start_index = self.data['cis_indices'][chr_indices[chrint]]
                stop_index = self.data['cis_indices'][chr_indices[chrint + 1]]
                if start_index <= stop_index:
                    node_ranges = numpy.round(numpy.linspace(start_index, stop_index,
                                              self.num_procs + 1)).astype(numpy.int64)
                    data = self.data['cis_data'][node_ranges[self.rank]:node_ranges[self.rank + 1], :]
                    _binning.binning_bin_cis_observed(data,
                                                      filt,
                                                      mids,
                                                      bin_counts,
                                                      all_indices,
                                                      distance_cutoffs,
                                                      num_bins,
                                                      bin_divs,
                                                      distance_div,
                                                      distance_bins,
                                                      mindistance,
                                                      maxdistance)
                # Find number of possible interactions in each bin
                maxfend = chr_indices[chrint + 1]
                if self.num_procs == 1:
                    startfend = chr_indices[chrint]
                    stopfend = maxfend
                else:
                    numfends = maxfend - chr_indices[chrint]
                    fend_sizes = numpy.r_[0, numfends - numpy.arange(1, numfends)]
                    for i in range(1, fend_sizes.shape[0]):
                        fend_sizes[i] += fend_sizes[i - 1]
                    node_targets = numpy.round(numpy.linspace(0, numfends * (numfends - 1) / 2,
                                               self.num_procs + 1)).astype(numpy.int64)
                    node_ranges = numpy.searchsorted(fend_sizes[1:], node_targets).astype(numpy.int64)
                    if (node_ranges[self.rank] < fend_sizes.shape[0] - 1 and
                        node_targets[self.rank] - fend_sizes[node_ranges[self.rank]] > 
                        fend_sizes[node_ranges[self.rank] + 1] - node_targets[self.rank]):
                        node_ranges[self.rank] += 1
                    if (node_ranges[self.rank + 1] < fend_sizes.shape[0] - 1 and
                        node_targets[self.rank + 1] - fend_sizes[node_ranges[self.rank + 1]] > 
                        fend_sizes[node_ranges[self.rank + 1] + 1] - node_targets[self.rank + 1]):
                        node_ranges[self.rank + 1] += 1
                    startfend = node_ranges[self.rank] + chr_indices[chrint]
                    stopfend = node_ranges[self.rank + 1] + chr_indices[chrint]
                _binning.binning_bin_cis_expected(filt,
                                                  mids,
                                                  bin_counts,
                                                  all_indices,
                                                  distance_cutoffs,
                                                  num_bins,
                                                  bin_divs,
                                                  distance_div,
                                                  distance_bins,
                                                  mindistance,
                                                  maxdistance,
                                                  startfend,
                                                  stopfend,
                                                  maxfend)
            if usereads in ['trans', 'all']:
                # Find number of observations in each bin
                start_index = self.data['trans_indices'][chr_indices[chrint]]
                stop_index = self.data['trans_indices'][chr_indices[chrint + 1]]
                if start_index < stop_index:
                    node_ranges = numpy.round(numpy.linspace(start_index, stop_index,
                                              self.num_procs + 1)).astype(numpy.int64)
                    data = self.data['trans_data'][node_ranges[self.rank]:node_ranges[self.rank + 1], :]
                    _binning.binning_bin_trans_observed(data,
                                                        filt,
                                                        bin_counts,
                                                        all_indices,
                                                        num_bins,
                                                        bin_indices[1:] - bin_indices[:-1],
                                                        bin_divs,
                                                        distance_div,
                                                        distance_bins)
                valid = numpy.where(filt[chr_indices[chrint]:chr_indices[chrint + 1]] == 1)[0] + chr_indices[chrint]
                if valid.shape[0] > 0:
                    node_ranges = numpy.round(numpy.linspace(0, valid.shape[0],
                                              self.num_procs + 1)).astype(numpy.int64)
                    startfend = valid[node_ranges[self.rank]]
                    stopfend = valid[node_ranges[self.rank + 1] - 1] + 1
                    for chrom2 in chroms[(h + 1):]:
                        chrint2 = self.chr2int[chrom2]
                        _binning.binning_bin_trans_expected(filt,
                                                            bin_counts,
                                                            all_indices,
                                                            num_bins,
                                                            bin_divs,       
                                                            distance_div,
                                                            distance_bins,
                                                            startfend,
                                                            stopfend,
                                                            chr_indices[chrint2],
                                                            chr_indices[chrint2 + 1])
        # adjust fend indices to position them along correction array
        for i in range(correction_indices.shape[0] - 1):
            all_indices[:, i, 1] += correction_indices[i]
        self.binning_fend_indices = all_indices
        self.binning_correction_indices = correction_indices
        # Exchange bin_counts
        if self.rank == 0:
            for i in range(1, self.num_procs):
                self.comm.Recv(temp_counts, source=i, tag=13)
                bin_counts += temp_counts
            for i in range(1, self.num_procs):
                self.comm.Send(bin_counts, dest=i, tag=13)
        else:
            self.comm.Send(bin_counts, dest=0, tag=13)
            self.comm.Recv(bin_counts, source=0, tag=13)
        # if using pseudo counts, find how many should be added
        # when bin1 != bin2, the number should be doubled to represent both combinations
        if not pseudocounts is None and pseudocounts > 0:
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding psuedo-counts...") % (' ' * 80),
            pcounts = numpy.zeros(bin_counts.shape[0], dtype=numpy.int32) + pseudocounts * 2 ** len(model)
            for i in range(len(model)):
                indices = ((numpy.arange(num_bins[i]) * (num_bins[i]) - num_bins[i] * (num_bins[i] - 1) / 2) *
                            bin_divs[i])
                other_counts = bin_counts.shape[0] / (num_bins[i] * (num_bins[i] + 1) / 2)
                for j in range(other_counts):
                    pcounts[indices * other_counts / bin_divs[i]] /= 2
            bin_counts += pcounts.reshape(-1, 1)
        # Find seed values
        if not self.silent:
            print >> sys.stderr, ("\r%s\rFinding seed values...") % (' ' * 80),
        prior = numpy.sum(bin_counts[:, 0]) / numpy.sum(bin_counts[:, 1]).astype(numpy.float64)
        log_prior = numpy.log2(prior)
        log_2 = numpy.log(2.0)
        min_p = 1.0 / numpy.sum(bin_counts[:, 1])
        all_indices = numpy.zeros((bin_counts.shape[0], len(model)), dtype=numpy.int32)
        for i in range(correction_indices.shape[0] - 1):
            all_indices[:, i] = ((numpy.arange(bin_counts.shape[0], dtype=numpy.int32) / bin_divs[i]) %
                                 (correction_indices[i + 1] - correction_indices[i]))
            temp0 = numpy.bincount(all_indices[:, i], weights=bin_counts[:, 0], minlength=num_bins[i])
            temp1 = numpy.bincount(all_indices[:, i], weights=bin_counts[:, 1], minlength=num_bins[i])
            all_corrections[correction_indices[i]:correction_indices[i + 1]] = (
                   numpy.maximum(1, temp0) / numpy.maximum(1, temp1).astype(numpy.float64) / prior)
        if distance_bins > 0:
            distance_indices = ((numpy.arange(bin_counts.shape[0], dtype=numpy.int32) / distance_div) %
                                distance_corrections.shape[0])
            distance_corrections = (
                numpy.bincount(distance_indices, weights=bin_counts[:, 0],
                minlength=distance_corrections.shape[0]).astype(numpy.float64) /
                (numpy.maximum(1, numpy.bincount(distance_indices, weights=bin_counts[:, 1],
                minlength=distance_corrections.shape[0])) * prior)).astype(numpy.float64)
        else:
            distance_indices = None
        bin_counts[:, 1] -= bin_counts[:, 0]

        def find_ll(indices, corrections, correction_indices, distance_c, distance_i, counts, prior, min_p):
            prod = numpy.ones(counts.shape[0], dtype=numpy.float64)
            for i in range(indices.shape[1]):
                prod *= corrections[correction_indices[i] + indices[:, i]]
            if not distance_c is None:
                prod *= distance_c[distance_i]
            prod *= prior
            prod = numpy.minimum(1.0 - min_p, numpy.maximum(prod, min_p))
            sum_log = numpy.log2(prod)
            return (-numpy.sum(counts[:, 0] * sum_log + counts[:, 1] * numpy.log2(1.0 - prod)))

        def find_sum_log_prod(index, indices, corrections, correction_indices, distance_c, distance_i, sum_log, prod):
            prod.fill(1.0)
            for i in range(indices.shape[1]):
                if i == index:
                    continue
                prod *= corrections[correction_indices[i] + indices[:, i]]
            if not distance_c is None:
                prod *= distance_c[distance_i]
            sum_log[:] = numpy.log2(prod)
            return

        def temp_ll(x, *args):
            counts, sum_log, prod, prior, log_prior, log_2, min_p = args[:7]
            return -numpy.sum(counts[:, 0] * (log_prior + sum_log + numpy.log2(x[0])) +
                             counts[:, 1] * numpy.log2(numpy.maximum(min_p, 1.0 - prior * prod * x[0])))

        def temp_ll_grad(x, *args):
            counts, sum_log, prod, prior, log_prior, log_2, min_p = args[:7]
            grad = numpy.array(-numpy.sum(counts[:, 0] / (log_2 * x[0]) - counts[:, 1] *
                               prior * prod / (log_2 * numpy.maximum(min_p, 1.0 - prior * prod * x[0]))), dtype=numpy.float64)
            return grad

        ll = find_ll(all_indices, all_corrections, correction_indices, distance_corrections, distance_indices,
                     bin_counts, prior, min_p)
        if not self.silent:
            print >> sys.stderr, ("\r%s\rLearning binning corrections... iteration:00  ll:%f ") % (' ' * 80, ll),
        iteration = 0
        delta = numpy.inf
        pgtol = 1e-8

        old_settings = numpy.seterr(invalid='ignore', divide='ignore')
        while iteration < max_iterations and delta >= learning_threshold:
            new_corrections = numpy.copy(all_corrections)
            for h in range(len(model)):
                # don't learn constant parameters
                if parameters[h].count('const') > 0:
                    continue
                # divide bins among nodes
                num_cor = correction_indices[h + 1] - correction_indices[h]
                node_ranges = numpy.round(numpy.linspace(0, num_cor, self.num_procs + 1)).astype(numpy.int64)
                bins = bin_counts.shape[0] / num_cor
                temp_sum_log = numpy.zeros(bins, dtype=numpy.float64)
                temp_prod = numpy.ones(bins, dtype=numpy.float64)
                for i in range(node_ranges[self.rank], node_ranges[self.rank + 1]):
                    where = numpy.where(all_indices[:, h] == i)[0]
                    temp_bin_counts = bin_counts[where, :]
                    if not distance_corrections is None:
                        temp_dist_indices = distance_indices[where]
                    else:
                        temp_dist_indices = None
                    x0 = all_corrections[(correction_indices[h] + i):(correction_indices[h] + i + 1)]
                    find_sum_log_prod(h, all_indices[where, :], all_corrections, correction_indices,
                                      distance_corrections, temp_dist_indices, temp_sum_log,
                                      temp_prod)
                    x, f, d = bfgs(func=temp_ll, x0=x0, fprime=temp_ll_grad, pgtol=pgtol,
                                   args=(temp_bin_counts, temp_sum_log, temp_prod, prior, log_prior, log_2, min_p))
                    new_corrections[correction_indices[h] + i] = x[0]
                if self.rank == 0:
                    for i in range(1, self.num_procs):
                        if node_ranges[i + 1] > node_ranges[i]:
                            self.comm.Recv(new_corrections[(node_ranges[i] + correction_indices[h]):(
                                           node_ranges[i + 1] + correction_indices[h])], source=i, tag=13)
                else:
                    if node_ranges[self.rank + 1] > node_ranges[self.rank]:
                        self.comm.Send(new_corrections[(node_ranges[self.rank] + correction_indices[h]):(
                                       node_ranges[self.rank + 1] + correction_indices[h])], dest=0, tag=13)
            if self.rank == 0:
                all_corrections = new_corrections
                for i in range(1, self.num_procs):
                    self.comm.Send(all_corrections, dest=i, tag=13)
            else:
                self.comm.Recv(all_corrections, source=0, tag=13)
            iteration += 1
            new_ll = find_ll(all_indices, all_corrections, correction_indices, distance_corrections,
                             distance_indices, bin_counts, prior, min_p)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rLearning binning corrections... iteration:%02i  ll:%f\n") % (' ' * 80,
                                      iteration, new_ll),
            delta = ll - new_ll
            if delta < 0.0:
                delta = numpy.inf
            ll = new_ll
        numpy.seterr(**old_settings)
        if not self.silent:
            print >> sys.stderr, ("\r%s\rLearning binning corrections... Final ll:%f\n") % (' ' * 80, ll),
        self.normalization = 'binning'
        self.binning_corrections = all_corrections.astype(numpy.float32)
        self.model_parameters = numpy.array(model)
        self.history += "Success\n"
        return None

    def find_trans_means(self):
        """
        Calculate the mean signals across all valid fend-pair trans interactions for each chromosome pair.

        :returns: None

        :Attributes: * **trans_mean** (*float*) - A float corresponding to the mean signal of inter-chromosome interactions.
        """
        self.history += "HiC.find_trans_mean() - "
        if not self.silent:
            print >> sys.stderr, ("Finding mean signals across trans interactions..."),
        chr_indices = self.fends['chr_indices'][...]
        num_chroms = chr_indices.shape[0] - 1
        possible = numpy.zeros(num_chroms * (num_chroms - 1) / 2, dtype=numpy.int32)
        pos = 0
        for i in range(chr_indices.shape[0] - 2):
            valid1 = numpy.sum(self.filter[chr_indices[i]:chr_indices[i + 1]])
            for j in range(i + 1, chr_indices.shape[0] - 1):
                valid2 = numpy.sum(self.filter[chr_indices[j]:chr_indices[j + 1]])
                possible[pos] = valid1 * valid2
                pos += 1
        trans_data = self.data['trans_data'][...]
        valid = numpy.where(self.filter[trans_data[:, 0]] * self.filter[trans_data[:, 1]])[0]
        trans_data = trans_data[valid, :]
        del valid
        chrom = self.fends['fends']['chr'][trans_data[:, 0]]
        indices = chrom * (num_chroms - 1) - (chrom * (chrom + 1) / 2) - 1
        indices += self.fends['fends']['chr'][trans_data[:, 1]]
        del chrom
        counts = trans_data[:, 2] / (self.corrections[trans_data[:, 0]] * self.corrections[trans_data[:, 1]])
        del trans_data
        actual = numpy.bincount(indices, weights=counts, minlength=possible.shape[0])
        self.trans_means = actual / numpy.maximum(1.0, possible.astype(numpy.float32))
        if not self.silent:
            print >> sys.stderr, ('Done\n'),
        self.history += "Success\n"
        return None

    def learn_fend_3D_model(self, chrom, minobservations=10):
        """
        Learn coordinates for a 3D model of data using an approximate PCA dimensional reduction.

        This function makes use of the :mod:`mlpy` function :func:`PCAFast` to reduce the data to a set of three coordinates per fend. Cis data for all unfiltered fends for the specified chromosome are dynamically binned to yield a complete distance matrix. The diagonal is set equal to the highest valid enrichment value after dynamic binning. This N x N matrix is passed to :func:`PCAFast` and reduced to an N x 3 matrix.

        :param chrom: The chromosome to learn the model for.
        :type chrom: str.
        :param minobservations: The minimum number of observed reads needed to cease bin expansion in the dynamic binning phase.
        :type minobservations: int.
        :returns: Array containing a row for each valid fend and columns containing X coordinate, Y coordinate, Z coordinate, and sequence coordinate (fend midpoint).
        """
        if 'mlpy' not in sys.modules.keys():
            if not self.silent:
                print >> sys.stderr, ("The mlpy module must be installed to use this function.")
            return None
        if not self.silent:
            print >> sys.stderr, ("Learning fend-resolution 3D model..."),
        unbinned, mapping = hic_binning.unbinned_cis_signal(self, chrom, datatype='fend', arraytype='upper',
                                                            skipfiltered=True, returnmapping=True)
        mids = (mapping[:, 0] + mapping[:, 1]) / 2
        if not self.silent:
            print >> sys.stderr, ("Dynamically binning data..."),
        dynamic = numpy.copy(unbinned)
        hic_binning.dynamically_bin_unbinned_upper(unbinned, mids, dynamic, minobservations)
        dynamic = numpy.log(dynamic[:, 0], dynamic[:, 1])
        data_mat = numpy.zeros((mids.shape[0], mids.shape[0]), dtype=numpy.float32)
        indices = numpy.triu_indices(mids.shape[0], 1)
        data_mat[indices] = dynamic
        data_mat[indices[1], indices[0]] = dynamic
        del dynamic
        del indices
        if not self.silent:
            print >> sys.stderr, ("Done\nModeling chromosome..."),
        pca_fast = mlpy.PCAFast(k=3)
        pca_fast.learn(data_mat)
        coordinates = pca_fast.transform(data_mat)
        if not self.silent:
            print >> sys.stderr, ("Done\n"),
        return numpy.hstack((coordinates, mids.reshape(-1, 1)))

    def cis_heatmap(self, chrom, start=None, stop=None, startfend=None, stopfend=None, binsize=0, binbounds=None,
                    datatype='enrichment', arraytype='compact', maxdistance=0, skipfiltered=False, returnmapping=False,
                    dynamically_binned=False, minobservations=0, searchdistance=0, expansion_binsize=0,
                    removefailed=False, image_file=None, **kwargs):
        """
        Return a heatmap of cis data of the type and shape specified by the passed arguments.

        This function returns a heatmap for a single chromosome region, bounded by either 'start' and 'stop' or 'startfend' and 'stopfend' ('start' and 'stop' take precedence), or if given, the outer coordinates of the array passed by 'binbounds'. If none of these are specified, data for the complete chromosome is used. The data in the array is determined by the 'datatype', being raw, fend-corrected, distance-corrected, enrichment, or expected data. The array shape is given by 'arraytype' and can be compact, upper, or full. See :mod:`hic_binning <hifive.hic_binning>` for further explanation of 'datatype' and 'arraytype'. The returned data will include interactions ranging from zero to 'maxdistance' apart. If maxdistance is zero, all interactions within the requested bounds are returned. If using dynamic binning ('dynamically_binned' is set to True), 'minobservations', 'searchdistance', 'expansion_binsize', and 'removefailed' are used to control the dynamic binning process. Otherwise these arguments are ignored.

        :param chrom: The name of a chromosome to obtain data from.
        :type chrom: str.
        :param start: The smallest coordinate to include in the array, measured from fend midpoints. If both 'start' and 'startfend' are given, 'start' will override 'startfend'. If unspecified, this will be set to the midpoint of the first fend for 'chrom'. Optional.
        :type start: int.
        :param stop: The largest coordinate to include in the array, measured from fend midpoints. If both 'stop' and 'stopfend' are given, 'stop' will override 'stopfend'. If unspecified, this will be set to the midpoint of the last fend plus one for 'chrom'. Optional.
        :type stop: int.
        :param startfend: The first fend to include in the array. If unspecified and 'start' is not given, this is set to the first fend in 'chrom'. In cases where 'start' is specified and conflicts with 'startfend', 'start' is given preference. Optional
        :type startfend: int.
        :param stopfend: The first fend not to include in the array. If unspecified and 'stop' is not given, this is set to the last fend in 'chrom' plus one. In cases where 'stop' is specified and conflicts with 'stopfend', 'stop' is given preference. Optional.
        :type stopfend: str.
        :param binsize: This is the coordinate width of each bin. If 'binsize' is zero, unbinned data is returned. If 'binbounds' is not None, this value is ignored.
        :type binsize: int.
        :param binbounds: An array containing start and stop coordinates for a set of user-defined bins. Any fend not falling in a bin is ignored. Optional.
        :type binbounds: numpy array
        :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend', 'enrichment', and 'expected'. Observed values are always in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, unfiltered fends return value of one. Expected values are returned for 'distance', 'fend', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fend' uses only fend correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
        :type datatype: str.
        :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'compact', 'full', and 'upper'. 'compact' means data are arranged in a N x M x 2 array where N is the number of fends or bins, M is the maximum number of steps between included fend pairs or bin pairs and data are stored such that bin n,m contains the interaction values between n and n + m + 1. 'full' returns a square, symmetric array of size N x N x 2. 'upper' returns only the flattened upper triangle of a full array, excluding the diagonal of size (N * (N - 1) / 2) x 2.
        :type arraytype: str.
        :param maxdistance: This specifies the maximum coordinate distance between bins that will be included in the array. If set to zero, all distances are included.
        :type maxdistance: str.
        :param skipfiltered: If 'True', all interaction bins for filtered out fends are removed and a reduced-size array is returned.
        :type skipfiltered: bool.
        :param returnmapping: If 'True', a list containing the data array and a 1d array containing fend numbers included in the data array if unbinned or a 2d array of N x 4 containing the first fend and last fend plus one included in each bin and first and last coordinates if binned is return. Otherwise only the data array is returned.
        :type returnmapping: bool.
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
        :param image_file: If a filename is specified, a PNG image file is written containing the heatmap data. Arguments for the appearance of the image can be passed as additional keyword arguments.
        :type image_file: str.
        :returns: Array in format requested with 'arraytype' containing data requested with 'datatype'. If returnmapping is True, a list is returned containined the requested data array and an array of associated positions (dependent on the binning options selected).

        """
        # check that all values are acceptable
        if datatype not in ['raw', 'fend', 'distance', 'enrichment', 'expected']:
            if not self.silent:
                print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
            return None
        if arraytype not in ['full', 'compact', 'upper']:
            if not self.silent:
                print >> sys.stderr, ("Unrecognized array type. No data returned.\n"),
            return None
        if dynamically_binned and arraytype == 'full':
            dynamic_arraytype = 'upper'
        else:
            dynamic_arraytype = arraytype
        # determine if data is to be dynamically binned
        if not dynamically_binned:
            data = hic_binning.find_cis_signal(self, chrom, binsize=binsize, binbounds=binbounds, start=start,
                                               stop=stop, startfend=startfend, stopfend=stopfend, datatype=datatype,
                                               arraytype=arraytype, maxdistance=maxdistance, skipfiltered=skipfiltered,
                                               returnmapping=returnmapping, silent=self.silent)
        else:
            if not binbounds is None:
                estart = binbounds[0, 0]
                estop = binbounds[-1, 1]
            else:
                estart = start
                estop = stop
            expansion, exp_mapping = hic_binning.find_cis_signal(self, chrom, start=estart, stop=estop,
                                                                 startfend=startfend, stopfend=stopfend,
                                                                 binsize=expansion_binsize, binbounds=None,
                                                                 datatype=datatype, arraytype=dynamic_arraytype,
                                                                 maxdistance=maxdistance, returnmapping=True,
                                                                 silent=self.silent)
            binned, mapping = hic_binning.find_cis_signal(self, chrom, start=start, stop=stop, startfend=startfend,
                                                          stopfend=stopfend, binsize=binsize, binbounds=binbounds,
                                                          datatype=datatype, arraytype=dynamic_arraytype,
                                                          maxdistance=maxdistance, returnmapping=True,
                                                          silent=self.silent)
            hic_binning.dynamically_bin_cis_array(expansion, exp_mapping, binned, mapping,
                                                  minobservations=minobservations, searchdistance=searchdistance,
                                                  removefailed=removefailed, silent=self.silent)
            if arraytype == 'full':
                indices = numpy.triu_indices(mapping.shape[0], 1)
                full_binned = numpy.zeros((mapping.shape[0], mapping.shape[0], 2), dtype=numpy.float32)
                full_binned[indices[0], indices[1], :] = binned
                full_binned[indices[1], indices[0], :] = binned
                binned = full_binned
            if returnmapping:
                data = [binned, mapping]
            else:
                data = binned
        if not image_file is None:
            if 'symmetricscaling' not in kwargs:
                if datatype == 'enrichment':
                    kwargs['symmetricscaling'] = True
                else:
                    kwargs['symmetricscaling'] = False
            if isinstance(data, list):
                binned = data[0]
            else:
                binned = data
            if arraytype == 'compact':
                img = plotting.plot_compact_array(binned, silent=self.silent, **kwargs)
            elif arraytype == 'full':
                img = plotting.plot_full_array(binned, silent=self.silent, **kwargs)
            else:
                img = plotting.plot_upper_array(binned, silent=self.silent, **kwargs)
            img.save(image_file, format='png')
        return data

    def trans_heatmap(self, chrom1, chrom2, start1=None, stop1=None, startfend1=None, stopfend1=None, binbounds1=None,
                      start2=None, stop2=None, startfend2=None, stopfend2=None, binbounds2=None, binsize=1000000,
                      skipfiltered=False, datatype='enrichment', returnmapping=False, dynamically_binned=False,
                      minobservations=0, searchdistance=0, expansion_binsize=0, removefailed=False, image_file=None,
                      **kwargs):
        """
        Return a heatmap of trans data of the type and shape specified by the passed arguments.

        This function returns a heatmap for trans interactions between two chromosomes within a region, bounded by either 'start1', 'stop1', 'start2' and 'stop2' or 'startfend1', 'stopfend1', 'startfend2', and 'stopfend2' ('start' and 'stop' take precedence), or if given, the outer coordinates of the arrays passed by 'binbounds1' and 'binbounds2'. The data in the array is determined by the 'datatype', being raw, fend-corrected, distance-corrected, enrichment, or expected data. The array shape is always rectangular. See :mod:`hic_binning <hifive.hic_binning>` for further explanation of 'datatype'. If using dynamic binning ('dynamically_binned' is set to True), 'minobservations', 'searchdistance', 'expansion_binsize', and 'removefailed' are used to control the dynamic binning process. Otherwise these arguments are ignored.

        :param chrom1: The name of the first chromosome to obtain data from.
        :type chrom1: str.
        :param chrom2: The name of the second chromosome to obtain data from.
        :type chrom2: str.
        :param start1: The coordinate at the beginning of the smallest bin from 'chrom1'. If unspecified, 'start1' will be the first multiple of 'binsize' below the 'startfend1' mid. If there is a conflict between 'start1' and 'startfend1', 'start1' is given preference. Optional.
        :type start1: int.
        :param stop1: The largest coordinate to include in the array from 'chrom1', measured from fend midpoints. If both 'stop1' and 'stopfend1' are given, 'stop1' will override 'stopfend1'. 'stop1' will be shifted higher as needed to make the last bin of size 'binsize'. Optional.
        :type stop1: int.
        :param startfend1: The first fend from 'chrom1' to include in the array. If unspecified and 'start1' is not given, this is set to the first valid fend in 'chrom1'. In cases where 'start1' is specified and conflicts with 'startfend1', 'start1' is given preference. Optional
        :type startfend1: int.
        :param stopfend1: The first fend not to include in the array from 'chrom1'. If unspecified and 'stop1' is not given, this is set to the last valid fend in 'chrom1' + 1. In cases where 'stop1' is specified and conflicts with 'stopfend1', 'stop1' is given preference. Optional.
        :param binbounds1: An array containing start and stop coordinates for a set of user-defined bins to use for partitioning 'chrom1'. Any fend not falling in a bin is ignored.
        :type binbounds1: numpy array
        :param start2: The coordinate at the beginning of the smallest bin from 'chrom2'. If unspecified, 'start2' will be the first multiple of 'binsize' below the 'startfend2' mid. If there is a conflict between 'start2' and 'startfend2', 'start2' is given preference. Optional.
        :type start2: int.
        :param stop2: The largest coordinate to include in the array from 'chrom2', measured from fend midpoints. If both 'stop2' and 'stopfend2' are given, 'stop2' will override 'stopfend2'. 'stop2' will be shifted higher as needed to make the last bin of size 'binsize'. Optional.
        :type stop2: int.
        :param startfend2: The first fend from 'chrom2' to include in the array. If unspecified and 'start2' is not given, this is set to the first valid fend in 'chrom2'. In cases where 'start2' is specified and conflicts with 'startfend2', 'start2' is given preference. Optional
        :type startfend2: int.
        :param stopfend2: The first fend not to include in the array from 'chrom2'. If unspecified and 'stop2' is not given, this is set to the last valid fend in 'chrom2' + 1. In cases where 'stop2' is specified and conflicts with 'stopfend2', 'stop1' is given preference. Optional.
        :type stopfend2: str.
        :param binbounds2: An array containing start and stop coordinates for a set of user-defined bins to use for partitioning 'chrom2'. Any fend not falling in a bin is ignored.
        :type binbounds2: numpy array
        :param binsize: This is the coordinate width of each bin. If binbounds is not None, this value is ignored.
        :type binsize: int.
        :param skipfiltered: If 'True', all interaction bins for filtered out fends are removed and a reduced-size array is returned.
        :type skipfiltered: bool.
        :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend', 'enrichment', and 'expected'. Observed values are always in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, unfiltered fends return value of one. Expected values are returned for 'distance', 'fend', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fend' uses only fend correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
        :type datatype: str.
        :param returnmapping: If 'True', a list containing the data array and two 2d arrays of N x 4 containing the first fend and last fend plus one included in each bin and first and last coordinates for the first and second chromosomes is returned. Otherwise only the data array is returned.
        :type returnmapping: bool.
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
        :param image_file: If a filename is specified, a PNG image file is written containing the heatmap data. Arguments for the appearance of the image can be passed as additional keyword arguments.
        :type image_file: str.
        :returns: Array in format requested with 'arraytype' containing data requested with 'datatype'. If returnmapping is True, a list is returned containined the requested data array and an array of associated positions (dependent on the binning options selected).
        """
        # check that all values are acceptable
        if datatype not in ['raw', 'fend', 'distance', 'enrichment', 'expected']:
            if not self.silent:
                print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
            return None
        # determine if data is to be dynamically binned
        if not dynamically_binned:
            data = hic_binning.find_trans_signal(self, chrom1, chrom2, start1=start1, stop1=stop1,
                                                startfend1=startfend1, stopfend1=stopfend1, binbounds1=binbounds1,
                                                start2=start2, stop2=stop2, startfend2=startfend2, stopfend2=stopfend2,
                                                binbounds2=binbounds2, binsize=binsize, datatype=datatype,
                                                returnmapping=returnmapping, skipfiltered=skipfiltered,
                                                silent=self.silent)
        else:
            if not binbounds1 is None:
                estart1 = binbounds1[0, 0]
                estop1 = binbounds1[-1, 1]
            else:
                estart1 = start1
                estop1 = stop1
            if not binbounds2 is None:
                estart2 = binbounds2[0, 0]
                estop2 = binbounds2[-1, 1]
            else:
                estart2 = start2
                estop2 = stop2
            expansion, exp_mapping1, exp_mapping2 = hic_binning.find_trans_signal(self, chrom1, chrom2, start1=estart1,
                                                                        stop1=estop1, startfend1=startfend1,
                                                                        stopfend1=stopfend1, binbounds1=binbounds1,
                                                                        start2=estart2, stop2=estop2,
                                                                        startfend2=startfend2, stopfend2=stopfend2,
                                                                        binbounds2=binbounds2,
                                                                        binsize=expansion_binsize, datatype=datatype,
                                                                        returnmapping=True, silent=self.silent)
            binned, mapping1, mapping2 = hic_binning.find_trans_signal(self, chrom1, chrom2, start1=start1,
                                                                       stop1=stop1, startfend1=startfend1,
                                                                       stopfend1=stopfend1, binbounds1=binbounds1,
                                                                       start2=start2, stop2=stop2,
                                                                       startfend2=startfend2, stopfend2=stopfend2,
                                                                       binbounds2=binbounds2, binsize=binsize,
                                                                       datatype=datatype, returnmapping=True,
                                                                       silent=self.silent)
            hic_binning.dynamically_bin_trans_array(expansion, exp_mapping1, exp_mapping2, binned, mapping1,
                                                    mapping2, minobservations=minobservations,
                                                    searchdistance=searchdistance, removefailed=removefailed,
                                                    silent=self.silent)
            if returnmapping:
                data = [binned, mapping1, mapping2]
            else:
                data = binned
        if not image_file is None:
            if 'symmetricscaling' not in kwargs:
                if datatype == 'enrichment':
                    kwargs['symmetricscaling'] = True
                else:
                    kwargs['symmetricscaling'] = False
            if isinstance(data, list):
                binned = data[0]
            else:
                binned = data
            img = plotting.plot_full_array(binned, silent=self.silent, **kwargs)
            img.save(image_file, format='png')
        return data

    def write_heatmap(self, filename, binsize, includetrans=True, datatype='enrichment', chroms=[], 
                      dynamically_binned=False, minobservations=0, searchdistance=0, expansion_binsize=0,
                      removefailed=False):
        """
        Create an h5dict file containing binned interaction arrays, bin positions, and an index of included chromosomes. This function is MPI compatible.

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

        The following attributes are created within the hdf5 dictionary file. Arrays are accessible as datasets while the resolution is held as an attribute.

        :Attributes: * **resolution** (*int.*) - The bin size that data are accumulated in.
                     * **chromosomes** (*ndarray*) - A numpy array of strings listing all of the chromosomes included in the heatmaps.
                     * **N.positions** (*ndarray*) - A series of numpy arrays of type int32, one for each chromosome where N is the chromosome name, containing one row for each bin and four columns denoting the start and stop coordinates and first fend and last fend plus one for each bin.
                     * **N.counts** (*ndarray*) - A series of numpy arrays of type int32, one for each chromosome where N is the chromosome name, containing the observed counts for valid fend combinations. Arrays are in an upper-triangle format such that they have N * (N - 1) / 2 entries where N is the number of fends or bins in the chromosome.
                     * **N.expected** (*ndarray*) - A series of numpy arrays of type float32, one for each chromosome where N is the chromosome name, containing the expected counts for valid fend combinations. Arrays are in an upper-triangle format such that they have N * (N - 1) / 2 entries where N is the number of fends in the chromosome.
                     * **N_by_M.counts** (*ndarray*) - A series of numpy arrays of type int32, one for each chromosome pair N and M if trans data are included, containing the observed counts for valid fend combinations. The chromosome name order specifies which axis corresponds to which chromosome.
                     * **N_by_M.expected** (*ndarray*) - A series of numpy arrays of type float32, one for each chromosome pair N and M if trans data are included, containing the expected counts for valid fend combinations. The chromosome name order specifies which axis corresponds to which chromosome.
        """
        history = self.history
        history += "HiC.write_heatmap(filename='%s', binsize=%i, includetrans=%s, datatype='%s', chroms=%s, dynamically_binned=%s, minobservations=%i, searchdistance=%i, expansion_binsize=%i, removefailed=%s)" % (filename, binsize, includetrans, datatype, str(chroms), dynamically_binned, minobservations, searchdistance, expansion_binsize, removefailed)
        if (chroms is None or
                (isinstance(chroms, list) and
                (len(chroms) == 0 or
                (len(chroms) == 1 and chroms[0] == ''))) or
                chroms == ''):
            chroms = self.chr2int.keys()
            chroms.sort()
        hic_binning.write_heatmap_dict(self, filename, binsize, includetrans=includetrans, datatype=datatype,
                                       chroms=chroms, dynamically_binned=dynamically_binned,
                                       minobservations=minobservations, searchdistance=searchdistance,
                                       expansion_binsize=expansion_binsize, removefailed=removefailed,
                                       silent=self.silent, history=history)
        return None
