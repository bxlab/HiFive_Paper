#!/usr/bin/env python

"""A class for handling HiC read data."""

import os
import sys

import numpy
import h5py
try:
    import pysam
except:
    pass


class HiCData(object):

    """This class handles interaction count data for HiC experiments.

    This class stores mapped paired-end reads, indexing them by fend-end (fend) number, in an h5dict.

    .. note::
      This class is also available as hifive.HiCData

    When initialized, this class creates an h5dict in which to store all data associated with this object.
    
    :param filename: The file name of the h5dict. This should end with the suffix '.hdf5'
    :type filename: str.
    :param mode: The mode to open the h5dict with. This should be 'w' for creating or overwriting an h5dict with name given in filename.
    :type mode: str.
    :param silent: Indicates whether to print information about function execution for this object.
    :type silent: bool.
    :returns: :class:`HiCData` class object.

    :Attributes: * **file** (*str.*) A string containing the name of the file passed during object creation for saving the object to.
                 * **silent** (*bool.*) - A boolean indicating whether to suppress all of the output messages.
                 * **history** (*str.*) - A string containing all of the commands executed on this object and their outcomes.
    """

    def __init__(self, filename, mode='r', silent=False):
        """Create a :class:`HiCData` object."""
        self.file = os.path.abspath(filename)
        self.silent = silent
        self.history = ''
        if mode != "w":
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

    def save(self):
        """
        Save analysis parameters to h5dict.

        :returns: None
        """
        self.history.replace("'None'", "None")
        datafile = h5py.File(self.file, 'w')
        for key in self.__dict__.keys():
            if key in ['file', 'chr2int', 'fends', 'silent', 'cuts']:
                continue
            elif self[key] is None:
                continue
            elif isinstance(self[key], numpy.ndarray):
                datafile.create_dataset(key, data=self[key])
            elif not isinstance(self[key], dict):
                datafile.attrs[key] = self[key]
        datafile.close()
        return None

    def load(self):
        """
        Load data from h5dict specified at object creation.

        Any call of this function will overwrite current object data with values from the last :func:`save` call.

        :returns: None
        """
        datafile = h5py.File(self.file, 'r')
        for key in datafile.keys():
            self[key] = numpy.copy(datafile[key])
        for key in datafile['/'].attrs.keys():
            self[key] = datafile['/'].attrs[key]
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

    def _find_cut_sites(self):
        self.cuts = []
        chroms = self.fends['chromosomes'][...]
        chr_indices = self.fends['chr_indices'][...]
        for i in range(len(chroms)):
            self.cuts.append(numpy.r_[
                self.fends['fends']['start'][chr_indices[i]:(chr_indices[i + 1])][::2],
                self.fends['fends']['stop'][chr_indices[i + 1] - 1]])
        return None

    def load_data_from_raw(self, fendfilename, filelist, maxinsert):
        """
        Read interaction counts from a text file(s) and place in h5dict.

        Files should contain both mapped ends of a read, one read per line, separated by tabs. Each line should be in the following format::

          chromosome1    coordinate1  strand1   chromosome2    coordinate2  strand2

        where strands are given by the characters '+' and '-'.

        :param fendfilename: This specifies the file name of the :class:`Fend` object to associate with the dataset.
        :type fendfilename: str.
        :param filelist: A list containing all of the file names of mapped read text files to be included in the dataset. If only one file is needed, this may be passed as a string.
        :type filelist: list
        :param maxinsert: A cutoff for filtering paired end reads whose total distance to their respective restriction sites exceeds this value.
        :type maxinsert: int.
        :returns: None

        :Attributes: * **fendfilename** (*str.*) - A string containing the relative path of the fend file.
                     * **cis_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero intra-chromosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend, the second column contains the idnex of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **cis_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'cis_data'. For example, all of the downstream cis interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :]. 
                     * **trans_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero inter-chroosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend (upstream also refers to the lower indexed chromosome in this context), the second column contains the index of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **trans_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'trans_data'. For example, all of the downstream trans interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :].
                     * **frags** (*ndarray*) - A filestream to the hdf5 fend file such that all saved fend attributes can be accessed through this class attribute.
                     * **maxinsert** (*int.*) - An interger denoting the maximum included distance sum between both read ends and their downstream RE site.

        When data is loaded the 'history' attribute is updated to include the history of the fend file that becomes associated with it.
        """
        self.history += "HiCData.load_data_from_raw(fendfilename='%s', filelist=%s, maxinsert=%i) - " % (fendfilename, str(filelist), maxinsert)
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            if not self.silent:
                print >> sys.stderr, \
                ("The fend file %s was not found. No data was loaded.\n") % (fendfilename),
            self.history += "Error: '%s' not found\n" % fendfilename
            return None
        self.fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                       os.path.dirname(self.file)), os.path.basename(fendfilename))
        self.maxinsert = maxinsert
        self.fends = h5py.File(fendfilename, 'r')
        self.history = self.fends['/'].attrs['history'] + self.history
        self.chr2int = {}
        chroms = self.fends['chromosomes'][...]
        for i, j in enumerate(chroms):
            self.chr2int[j] = i
        # load data from all files, skipping if chromosome not in the fend file.
        if isinstance(filelist, str):
            filelist = [filelist]
        fend_pairs = []
        for i in range(len(chroms)):
            fend_pairs.append([])
            for j in range(i + 1):
                fend_pairs[i].append({})
        total_reads = 0
        for fname in filelist:
            data = []
            for i in range(len(chroms)):
                data.append([])
                for j in range(i + 1):
                    data[i].append({})
            if not os.path.exists(fname):
                if not self.silent:
                    print >> sys.stderr, ("The file %s was not found...skipped.\n") % (fname.split('/')[-1]),
                self.history += "'%s' not found, " % fname
                continue
            if not self.silent:
                print >> sys.stderr, ("Loading data from %s...") % (fname.split('/')[-1]),
            input = open(fname, 'r', 1)
            for line in input:
                temp = line.strip('\n').split('\t')
                temp[0].strip('chr')
                temp[3].strip('chr')
                if temp[0] not in self.chr2int or temp[3] not in self.chr2int:
                    continue
                chrint1 = self.chr2int[temp[0]]
                chrint2 = self.chr2int[temp[3]]
                start1 = int(temp[1])
                start2 = int(temp[4])
                if chrint1 == chrint2:
                    if start1 < start2:
                        if temp[2] == '-':
                            start1 = -start1
                        if temp[5] == '-':
                            start2 = -start2
                        data[chrint2][chrint1][(start1, start2)] = None
                    else:
                        if temp[2] == '-':
                            start1 = -start1
                        if temp[5] == '-':
                            start2 = -start2
                        data[chrint1][chrint2][(start2, start1)] = None
                elif chrint1 < chrint2:
                    if temp[2] == '-':
                        start1 = -start1
                    if temp[5] == '-':
                        start2 = -start2
                    data[chrint2][chrint1][(start1, start2)] = None
                else:
                    if temp[2] == '-':
                        start1 = -start1
                    if temp[5] == '-':
                        start2 = -start2
                    data[chrint1][chrint2][(start2, start1)] = None
            input.close()
            new_reads = 0
            for i in range(len(data)):
                for j in range(len(data[i])):
                    data[i][j] = numpy.array(data[i][j].keys(), dtype=numpy.int32)
                    new_reads += data[i][j].shape[0]
            if not self.silent:
                print >> sys.stderr, ("%i validly-mapped reads pairs loaded.\n") % (new_reads),
            total_reads += new_reads
            # map data to fends, filtering as needed
            if new_reads > 0:
                self._find_fend_pairs(data, fend_pairs)
        self._clean_fend_pairs(fend_pairs)
        total_fend_pairs = 0
        for i in range(len(fend_pairs)):
            for j in range(len(fend_pairs[i])):
                total_fend_pairs += len(fend_pairs[i][j])
        if total_fend_pairs == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            self.history += "Error: no valid data loaded\n"
            return None
        if not self.silent:
            print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i valid fend pairs\n") %\
                             (total_reads, total_fend_pairs),
        # write fend pairs to h5dict
        self._parse_fend_pairs(fend_pairs)
        self.history += 'Success\n'
        return None

    def load_data_from_bam(self, fendfilename, filelist, maxinsert):
        """
        Read interaction counts from pairs of BAM-formatted alignment file(s) and place in h5dict.

        :param fendfilename: This specifies the file name of the :class:`Fend` object to associate with the dataset.
        :type fendfilename: str.
        :param filelist: A list containing lists of paired end bam files. If only one pair of files is needed, the list may contain both file path strings.
        :type filelist: list of mapped sequencing runs. Each run should be a list of the first and second read end bam files ([[run1_1, run1_2], [run2_1, run2_2]...])
        :param maxinsert: A cutoff for filtering paired end reads whose total distance to their respective restriction sites exceeds this value.
        :type maxinsert: int.
        :returns: None

        :Attributes: * **fendfilename** (*str.*) - A string containing the relative path of the fend file.
                     * **cis_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero intra-chromosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend, the second column contains the idnex of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **cis_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'cis_data'. For example, all of the downstream cis interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :]. 
                     * **trans_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero inter-chroosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend (upstream also refers to the lower indexed chromosome in this context), the second column contains the index of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **trans_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'trans_data'. For example, all of the downstream trans interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :].
                     * **frags** (*ndarray*) - A filestream to the hdf5 fend file such that all saved fend attributes can be accessed through this class attribute.
                     * **maxinsert** (*int.*) - An interger denoting the maximum included distance sum between both read ends and their downstream RE site.

        When data is loaded the 'history' attribute is updated to include the history of the fend file that becomes associated with it.
        """
        self.history += "HiCData.load_data_from_bam(fendfilename='%s', filelist=%s, maxinsert=%i) - " % (fendfilename, str(filelist), maxinsert)
        if 'pysam' not in sys.modules.keys():
            if not self.silent:
                print >> sys.stderr, ("The pysam module must be installed to use this function.")
            self.history += 'Error: pysam module missing\n'
            return None
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            if not self.silent:
                print >> sys.stderr, ("The fend file %s was not found. No data was loaded.\n") % (fendfilename),
            self.history += "Error: '%s' not found\n" % fendfilename
            return None
        self.fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                       os.path.dirname(self.file)), os.path.basename(fendfilename))
        self.maxinsert = maxinsert
        self.fends = h5py.File(fendfilename, 'r')
        self.history = self.fends['/'].attrs['history'] + self.history
        self.chr2int = {}
        chroms = self.fends['chromosomes'][...]
        for i, j in enumerate(chroms):
            self.chr2int[j] = i
        # load data from all files, skipping if chromosome not in the fend file.
        if isinstance(filelist[0], str):
            filelist = [[filelist[0], filelist[1]]]
        total_reads = 0
        fend_pairs = []
        for i in range(len(chroms)):
            fend_pairs.append([])
            for j in range(i + 1):
                fend_pairs[i].append({})
        for filepair in filelist:
            # determine which files have both mapped ends present
            present = True
            if not os.path.exists(filepair[0]):
                if not self.silent:
                    print >> sys.stderr, ("%s could not be located.") % (filepair[0]),
                self.history += "'%s' not found, " % filepair[0]
                present = False
            if not os.path.exists(filepair[1]):
                if not self.silent:
                    print >> sys.stderr, ("%s could not be located.") % (filepair[1]),
                self.history += "'%s' not found, " % filepair[1]
                present = False
            if not present:
                if not self.silent:
                    print >> sys.stderr, ("No data for one or both ends could be located. Skipping this run.\n")
                continue
            unpaired = {}
            # load first half of paired ends
            if not self.silent:
                print >> sys.stderr, ("Loading data from %s...") % (filepair[0].split('/')[-1]),
            input = pysam.Samfile(filepair[0], 'rb')
            idx2int = {}
            for i in range(len(input.header['SQ'])):
                chrom = input.header['SQ'][i]['SN'].strip('chr')
                if chrom in self.chr2int:
                    idx2int[i] = self.chr2int[chrom]
            for read in input.fetch(until_eof=True):
                # Only consider reads with an alignment
                if read.is_unmapped:
                    continue
                # if chromosome not in chr2int, skip
                if read.tid not in idx2int:
                    continue
                if read.is_reverse:
                    end = -(read.pos + len(read.seq))
                else:
                    end = read.pos

                unpaired[read.qname] = (idx2int[read.tid], end)
            input.close()
            if not self.silent:
                print >> sys.stderr, ("Done\n"),
            # load second half of paired ends
            if not self.silent:
                print >> sys.stderr, ("Loading data from %s...") % (filepair[1].split('/')[-1]),
            data = []
            for i in range(len(chroms)):
                data.append([])
                for j in range(i + 1):
                    data[i].append({})
            input = pysam.Samfile(filepair[1], 'rb')
            idx2int = {}
            for i in range(len(input.header['SQ'])):
                chrom = input.header['SQ'][i]['SN'].strip('chr')
                if chrom in self.chr2int:
                    idx2int[i] = self.chr2int[chrom]
            for read in input.fetch(until_eof=True):
                # Only consider reads with an alignment
                if read.is_unmapped:
                    continue
                # if chromosome not in chr2int, skip
                if read.tid not in idx2int:
                    continue
                if read.is_reverse:
                    start2 = -(read.pos + len(read.seq))
                else:
                    start2 = read.pos
                if read.qname not in unpaired:
                    continue
                chr1, start1 = unpaired[read.qname]
                chr2 = idx2int[read.tid]
                if chr1 == chr2:
                    if abs(start1) < abs(start2):
                        data[chr2][chr1][(start1, start2)] = None
                    else:
                        data[chr1][chr2][(start2, start1)] = None
                elif chr1 < chr2:
                    data[chr2][chr1][(start1, start2)] = None
                else:
                    data[chr1][chr2][(start2, start1)] = None
            input.close()
            del unpaired
            new_reads = 0
            for i in range(len(data)):
                for j in range(len(data[i])):
                    data[i][j] = numpy.array(data[i][j].keys(), dtype=numpy.int32)
                    new_reads += data[i][j].shape[0]
            if not self.silent:
                print >> sys.stderr, ("Done\nRead %i validly-mapped read pairs.\n") % (new_reads),
            total_reads += new_reads
            # map data to fends, filtering as needed
            if new_reads > 0:
                self._find_fend_pairs(data, fend_pairs)
            del data
        self._clean_fend_pairs(fend_pairs)
        total_fend_pairs = 0
        for i in range(len(fend_pairs)):
            for j in range(len(fend_pairs[i])):
                total_fend_pairs += len(fend_pairs[i][j])
        if total_fend_pairs == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            self.history += "Error: no valid data loaded\n"
            return None
        if not self.silent:
            print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i valid fend pairs\n") %\
                                 (total_reads, total_fend_pairs),
        self._parse_fend_pairs(fend_pairs)
        self.history += "Success\n"
        return None

    def load_data_from_mat(self, fendfilename, filename, maxinsert=0):
        """
        Read interaction counts from a :mod:`HiCPipe`-compatible 'mat' text file and place in h5dict.

        :param fendfilename: This specifies the file name of the :class:`Fend` object to associate with the dataset.
        :type fendfilename: str.
        :param filename: File name of a 'mat' file containing fend pair and interaction count data.
        :type filename: str.
        :param maxinsert: A cutoff for filtering paired end reads whose total distance to their respective restriction sites exceeds this value.
        :type maxinsert: int.
        :returns: None

        :Attributes: * **fendfilename** (*str.*) - A string containing the relative path of the fend file.
                     * **cis_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero intra-chromosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend, the second column contains the idnex of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **cis_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'cis_data'. For example, all of the downstream cis interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :]. 
                     * **trans_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero inter-chroosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend (upstream also refers to the lower indexed chromosome in this context), the second column contains the index of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **trans_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'trans_data'. For example, all of the downstream trans interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :].
                     * **frags** (*ndarray*) - A filestream to the hdf5 fend file such that all saved fend attributes can be accessed through this class attribute.
                     * **maxinsert** (*int.*) - An interger denoting the maximum included distance sum between both read ends and their downstream RE site.

        When data is loaded the 'history' attribute is updated to include the history of the fend file that becomes associated with it.
        """
        self.history += "HiCData.load_data_from_mat(fendfilename='%s', filename='%s', maxinsert=%i) - " % (fendfilename, filename, maxinsert)
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            if not self.silent:
                print >> sys.stderr, ("The fend file %s was not found. No data was loaded.\n") % (fendfilename),
            self.history += "Error: '%s' not found\n" % fendfilename
            return None
        self.fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                       os.path.dirname(self.file)), os.path.basename(fendfilename))
        self.maxinsert = maxinsert
        self.fends = h5py.File(fendfilename, 'r')
        self.history = self.fends['/'].attrs['history'] + self.history
        # load data from mat file. This assumes that the mat data was mapped
        # using the same fend numbering as in the fend file.
        chr_indices = self.fends['chr_indices'][...]
        fend_pairs = []
        for i in range(chr_indices.shape[0] - 1):
            fend_pairs.append([])
            for j in range(i + 1): 
                fend_pairs[i].append({})
        if not os.path.exists(filename):
            if not self.silent:
                print >> sys.stderr, ("%s not found... no data loaded.\n") % (filename.split('/')[-1]),
            self.history += "Error: '%s' not found\n" % fendfilename
            return None
        if not self.silent:
            print >> sys.stderr, ("Loading data from mat file..."),
        input = open(filename, 'r')
        for line in input:
            temp = line.strip('\n').split('\t')
            try:
                fend1 = int(temp[0]) - 1
            except ValueError:
                continue
            fend2 = int(temp[1]) - 1
            chr1 = numpy.searchsorted(chr_indices, fend1, side='right') - 1
            chr2 = numpy.searchsorted(chr_indices, fend2, side='right') - 1
            if chr2 < chr1:
                fend_pairs[chr1][chr2][(fend2 - chr_indices[chr2], fend1 - chr_indices[chr1])] = int(temp[2])
            elif chr1 < chr2:
                fend_pairs[chr2][chr1][(fend1 - chr_indices[chr1], fend2 - chr_indices[chr2])] = int(temp[2])
            elif fend1 < fend2:
                fend_pairs[chr2][chr1][(fend1 - chr_indices[chr1], fend2 - chr_indices[chr2])] = int(temp[2])
            else:
                fend_pairs[chr1][chr2][(fend2 - chr_indices[chr2], fend1 - chr_indices[chr1])] = int(temp[2])
        input.close()
        self._clean_fend_pairs(fend_pairs)
        total_fend_pairs = 0
        for i in range(len(fend_pairs)):
            for j in range(len(fend_pairs[i])):
                total_fend_pairs += len(fend_pairs[i][j])
        if total_fend_pairs == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            self.history += "Error: no valid data loaded\n"
            return None
        if not self.silent:
            print >> sys.stderr, ("%i valid fend pairs loaded.\n") % (total_fend_pairs),
        # write fend pairs to h5dict
        self._parse_fend_pairs(fend_pairs)
        self.history += "Success\n"
        return None

    def _find_fend_pairs(self, data, fend_pairs):
        """Return array with lower fend, upper fend, and count for pair."""
        if 'cuts' not in self.__dict__.keys():
            self._find_cut_sites()
        chroms = self.fends['chromosomes'][...]
        # assign fends on a per-chromosome pair basis
        for i in range(len(data)):
            for j in range(len(data[i])):
                if data[i][j].shape[0] == 0:
                    continue
                mapped_fends = numpy.empty((data[i][j].shape[0], 2), dtype=numpy.int32)
                distances = numpy.zeros(data[i][j].shape[0], dtype=numpy.int32)
                distances.fill(self.maxinsert + 1)
                signs = numpy.minimum(0, numpy.sign(data[i][j]))
                data[i][j] = numpy.abs(data[i][j])
                if not self.silent:
                    print >> sys.stderr, ("\rMapping fends for %s by %s") % (chroms[i].ljust(10), chroms[j].ljust(10)),
                mapped_fends[:, 0] = numpy.searchsorted(self.cuts[j], data[i][j][:, 0])
                mapped_fends[:, 1] = numpy.searchsorted(self.cuts[i], data[i][j][:, 1])
                valid = numpy.where((mapped_fends[:, 0] > 0) * (mapped_fends[:, 0] < self.cuts[j].shape[0]) *
                                    (mapped_fends[:, 1] > 0) * (mapped_fends[:, 1] < self.cuts[i].shape[0]))[0]
                # find distance from cutsite to mapped coordinate
                distances[valid] = (numpy.abs(data[i][j][valid, 0] - self.cuts[j][mapped_fends[valid, 0] +
                                    signs[valid, 0]]) + numpy.abs(data[i][j][valid, 1] -
                                    self.cuts[i][mapped_fends[valid, 1] + signs[valid, 1]]))
                # remove fends with too great an insert distance
                valid = numpy.where(distances <= self.maxinsert)[0]
                # convert from fragments to fends
                mapped_fends[valid, :2] = mapped_fends[valid, :2] * 2 - 1 + signs[valid, :]
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rCounting fend pairs...") % (' ' * 50),
                for k in valid:
                    key = (mapped_fends[k, 0], mapped_fends[k, 1])
                    if key not in fend_pairs[i][j]:
                        fend_pairs[i][j][key] = 1
                    else:
                        fend_pairs[i][j][key] += 1
        if not self.silent:
            print >> sys.stderr, ("Done\n"),
        return None

    def _clean_fend_pairs(self, fend_pairs):
        chr_indices = self.fends['chr_indices'][...]
        # remove fend pairs from same fend or opposite strand adjacents
        for i in range(len(fend_pairs)):
            for j in range(0, chr_indices[i + 1] - chr_indices[i] - 2, 2):
                # same fend
                name = (j, j)
                if name in fend_pairs[i][i]:
                    del fend_pairs[i][i][name]
                name = (j + 1, j + 1)
                if name in fend_pairs[i][i]:
                    del fend_pairs[i][i][name]
                # same fend
                name = (j, j + 1)
                if name in fend_pairs[i][i]:
                    del fend_pairs[i][i][name]
                # same fend
                name = (j + 1, j)
                if name in fend_pairs[i][i]:
                    del fend_pairs[i][i][name]
                # adjacent fends, opposite strands
                name = (j, j + 3)
                if name in fend_pairs[i][i]:
                    del fend_pairs[i][i][name]
                name = (j + 1, j + 2)
                if name in fend_pairs[i][i]:
                    del fend_pairs[i][i][name]
            j = chr_indices[i + 1] - chr_indices[i] - 2
            # same fend
            name = (j, j)
            if name in fend_pairs[i][i]:
                del fend_pairs[i][i][name]
            name = (j + 1, j + 1)
            if name in fend_pairs[i][i]:
                del fend_pairs[i][i][name]
            name = (j, j + 1)
            if name in fend_pairs[i][i]:
                del fend_pairs[i][i][name]
            name = (j + 1, j)
            if name in fend_pairs[i][i]:
                del fend_pairs[i][i][name]
        return None

    def _parse_fend_pairs(self, fend_pairs):
        """Separate fend pairs into cis and trans interactions index."""
        if not self.silent:
            print >> sys.stderr, ("Parsing fend pairs..."),
        chr_indices = self.fends['chr_indices'][...]
        # determine number of cis pairs
        cis_count = 0
        for i in range(len(fend_pairs)):
            cis_count += len(fend_pairs[i][i])
        # create cis array
        self.cis_data = numpy.empty((cis_count, 3), dtype=numpy.int32)
        pos = 0
        # fill in each chromosome's cis interactions
        for i in range(len(fend_pairs)):
            if len(fend_pairs[i][i]) == 0:
                continue
            chr_start = pos
            for key, count in fend_pairs[i][i].iteritems():
                self.cis_data[pos, :2] = key
                self.cis_data[pos, 2] = count
                pos += 1
            fend_pairs[i][i] = None
            # sort interactions
            order = numpy.lexsort((self.cis_data[chr_start:pos, 1], self.cis_data[chr_start:pos, 0]))
            self.cis_data[chr_start:pos, :] = self.cis_data[order + chr_start, :]
            self.cis_data[chr_start:pos, :2] += chr_indices[i]
            del order
        # determine number of trans pairs
        trans_count = 0
        for i in range(len(fend_pairs)):
            for j in range(i + 1, len(fend_pairs)):
                trans_count += len(fend_pairs[j][i])
        # create cis array
        self.trans_data = numpy.empty((trans_count, 3), dtype=numpy.int32)
        pos = 0
        # fill in each chromosome's trans interactions
        for i in range(len(fend_pairs) - 1):
            for j in range(i + 1, len(fend_pairs)):
                if len(fend_pairs[j][i]) == 0:
                    continue
                chr_start = pos
                for key, count in fend_pairs[j][i].iteritems():
                    self.trans_data[pos, :2] = key
                    self.trans_data[pos, 2] = count
                    pos += 1
                fend_pairs[j][i] = None
                # sort interactions
                order = numpy.lexsort((self.trans_data[chr_start:pos, 1], self.trans_data[chr_start:pos, 0]))
                self.trans_data[chr_start:pos, :] = self.trans_data[order + chr_start, :]
                self.trans_data[chr_start:pos, 0] += chr_indices[i]
                self.trans_data[chr_start:pos, 1] += chr_indices[j]
                del order
        # create data indices
        if self.cis_data.shape[0] > 0:
            cis_indices = numpy.r_[0, numpy.bincount(self.cis_data[:, 0],
                                   minlength=self.fends['fends'].shape[0])].astype(numpy.int64)
            for i in range(1, cis_indices.shape[0]):
                cis_indices[i] += cis_indices[i - 1]
            self.cis_indices = cis_indices
        else:
            self.cis_data = None
        if self.trans_data.shape[0] > 0:
            trans_indices = numpy.r_[0, numpy.bincount(self.trans_data[:, 0],
                                     minlength=self.fends['fends'].shape[0])].astype(numpy.int64)
            for i in range(1, trans_indices.shape[0]):
                trans_indices[i] += trans_indices[i - 1]
            self.trans_indices = trans_indices
        else:
            self.trans_data = None
        if not self.silent:
            print >> sys.stderr, ("Done  %i cis reads, %i trans reads\n") % (numpy.sum(self.cis_data[:, 2]),
                                                                             numpy.sum(self.trans_data[:, 2])),
        return None

    def export_to_mat(self, outfilename):
        """
        Write reads loaded in data object to text file in :mod:`HiCPipe`-compatible 'mat' format.

        :param outfilename: Specifies the file to save data in.
        :type outfilename: str.
        :returns: None
        """
        if not self.silent:
            print >> sys.stderr, ("Writing data to mat file..."),
        output = open(outfilename, 'w')
        if not self.silent:
            print >> output, "fend1\tfend2\tcount"
        for i in range(self.cis_indices.shape[0] - 1):
            for j in range(self.cis_indices[i], self.cis_indices[i + 1]):
                # One is added to indices so numbering starts from one.
                print >> output, "%i\t%i\t%i" % (i + 1, self.cis_data[j, 1] + 1, self.cis_data[j, 2])
            for j in range(self.trans_indices[i], self.trans_indices[i + 1]):
                print >> output, "%i\t%i\t%i" % (i + 1, self.trans_data[j, 1] + 1, self.trans_data[j, 2])
        output.close()
        if not self.silent:
            print >> sys.stderr, ("Done\n"),
        return None
