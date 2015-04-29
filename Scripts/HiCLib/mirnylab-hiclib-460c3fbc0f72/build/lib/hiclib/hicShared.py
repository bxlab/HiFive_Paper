import re
import warnings 
import numpy as np 
import logging
from scipy import stats
import os 
import mirnylib.h5dict
import mirnylib.numutils
import mirnylib.numutils_new 
from scipy.stats.stats import spearmanr
log = logging.getLogger(__name__)


try:
    from fastBinSearch import binarySearch as bs  # @UnresolvedImport
    fastBS = True
except:
    print "Please compile binarySearch!!! It is in the main folder of the library"
    fastBS = False

def fileIsFragment(filename):
    "Checks whether a file is a fragment-level h5dict"
    try:
        a = mirnylib.h5dict.h5dict(filename, "r")
    except:
        return False
    if "chrms1" in a.keys():
        return True
    return False


def fileIsHeatmap(filename):
    """Checks whether a file is a heatmap
    Returns a file type as a first argument, and resolution as a second
    """
    try:
        a = mirnylib.h5dict.h5dict(filename, "r")
    except:
        return False

    keys = a.keys()
    if ("heatmap" in keys) and ("resolution" in keys):
        resolution = getResolution(filename)
        assert resolution == a["resolution"]
        return ("heatmap", a["resolution"])
    elif "heatmap" in keys:
        warnings.warn("Resolution not found in {0}".format(filename))
        resolution = getResolution(filename)
        return ("heatmap", resolution)

    if "0 0" in keys:
        resolution = getResolution(filename)
        if "0 1" in keys:
            return "byChr", resolution
        else:
            return "byChrCis", resolution
 
    for key in keys:
        if len(key.split("_")) == 3:
            resolution = getResolution(filename)
            return "byChrSuper",resolution
    return False


def getResolution(fname):
    "Extracts resolution from a Hi-C filename"
    matches = re.findall(r"-\d{1,5}[k,M]", fname)
    if len(matches) == 0:
        raise ValueError("resolution not found")
    if len(matches) > 1:
        warnings.warn("undetermined resolution")
        print "found resolutions", matches
        print "using the last found"
    match = matches[-1]
    if match[-1] == "k":
        mult = 1000
    elif match[-1] == "M":
        mult = 1000000
    else:
        raise ValueError("Bad things happened")

    num = int(match[1:-1])
    return num * mult



def binarySearch(x, y):
    """
    ultrafast np.int64 binary search
    """
    if len(x) < 3000000:
        log.debug("Using regular searchsorted because dataset is short: {0}".format(len(x)))
        a =  np.searchsorted(y, x)
    else:
        if fastBS:
            log.debug("using fast binary search")
            a =  bs(x, y)            
            
        else:
            log.debug("Using searchsorted because fast binary search not found")            
            a = np.searchsorted(y, x)
    log.debug("Binary search finished")
    return a

# binarySearch = lambda x, y:np.searchsorted(y, x)
r_ = np.r_


def cleanFile(filename):
    if os.path.exists(filename):
        os.remove(filename)


class sliceableDataset(object):
    """
    Implements slicing for a function getFunction which has a syntax
    data[start:end] = getFunction(name,start,end), where name is a "key" for data
    """
    def __init__(self, getFunction, name, length):
        self.getFunction = getFunction
        self.name = name
        self.length = length
    def __getitem__(self, val):

        if issubclass(type(val), int):
            val = int(val)
            data = self.getFunction(self.name, val, val + 1)
            print "fetched one element from vectors2. This is bad!"
            return data[0]

        start = val.start
        stop = val.stop
        step = val.step
        if (step == None) or (step == 1):
            return self.getFunction(self.name, start, stop)
        else:
            ar = self.getFunction(self.name, start, stop)
            return ar[::step]
    def __array__(self):
        if self.length > 25000000:
            log.info("fetched numpy array as a whole. This is memory inefficient.")
        return self.getFunction(self.name)

    def __len__(self):
        return self.length


def corr(x, y):
    return stats.spearmanr(x, y)[0]


"""
Functions to work with fragmentHiC related numpy record arrays and with hdf5 datasets
"""

#defining data types for building heatmaps, and sorter functions for externalMergeSort 
mydtype = np.dtype("i1,i4,i1,i4,b,b")
mydtype.names = ("chrms1","pos1","chrms2","pos2","strands1","strands2")

def mydtypeSorter(x):
    """fuction which sorts mydtype type datasets over first two fields
    Parameters
    ----------
    x : numpy array of mydtype
        array to be sorted
    """
    inds = np.lexsort((x["pos1"],x["chrms1"]))
    toret = x.view(np.dtype((str, x.dtype.itemsize)))[inds].view(x.dtype)
    assert len(toret) == len(x)
    assert toret.dtype == x.dtype     
    return toret

def searchsorted(array, element):
    "matching searchsorted"
    c1 = array["chrms1"]
    p1 = array["pos1"]
    val1 = element[0]
    val2 = element[1]
    low = np.searchsorted(c1, val1,"left")
    high = np.searchsorted(c1,val1,"right")
    toret =  low + np.searchsorted(p1[low:high], val2, "right")     
    return toret 


def mycmp(i,j):
    "A cmp function for a tuple (chr, pos)"
    if i[0] > j[0]:
        return 1
    if i[0] < j[0]:
        return -1
    if i[1] > j[1]:
        return 1
    if i[1] < j[1]:
        return -1
    return 0
        

def h5dictBinarySearch(chrms1, pos1, value, side="left", mycmp = mycmp):
    """perform binary search in a sorted h5dict array of fragmentHiC data
    Parameters
    ----------
    chrms1 : h5py dataset
        h5py array of chromosomes (obtained with h5dict.getDataset)
    pos1 : h5py dataset
        h5py array of positions (cuts) (obtained h5dict.getDataset)
    values : tuple (chrms, pos) 
        values to search in h5dict dataset        
        
    """
    
    low = 0
    high = len(chrms1)-1
    if side == "left":
        more = [0,1]        
    elif side == "right":
        more = [1]        
    else:
        raise ValueError("side should be left or right")
    
        
    if mycmp((chrms1[0],pos1[0]), value) == 1:
        return 0 
    if mycmp((chrms1[-1],pos1[-1]), value) == -1:
        return len(chrms1)          
    
    while high - low > 1:                 
        mid = (low+high)//2        
        cmpResult = mycmp((chrms1[mid],pos1[mid]),value)
        if cmpResult  in more: 
            high = mid
        else: 
            low = mid         
    if side == "left":
        return low
    else:
        return high
    
""" Tools to manipulate with Hi-C datasets""" 
    
def saveFile(datasetFilename, outFolder, IC = True, smooth = True):
    if not os.path.exists(outFolder):
        os.mkdir(outFolder)
    mydict = mirnylib.h5dict.h5dict(datasetFilename)
    raw = os.path.join(outFolder,"raw")
    if not os.path.exists(raw):
        os.mkdir(raw)
    if IC:
        ICedFol = os.path.join(outFolder, "corrected")
        if not os.path.exists(ICedFol):
            os.mkdir(ICedFol)
    if smooth:
        smoothedFol = os.path.join(outFolder,"smoothed") 
        if not os.path.exists(smoothedFol):
            os.mkdir(smoothedFol)
    for key in mydict: 
        print key
        data = mydict[key] 
        key = key.replace(" ","_")
        if type(data) == np.ndarray and len(data.shape) == 2:
            
            mirnylib.numutils_new.matrixToGzippedFile(data, os.path.join(raw,key+".txt.gz"))  # @UndefinedVariable
            print "saved raw"
            if IC:
                ICed = mirnylib.numutils.completeIC(data)
                mirnylib.numutils_new.matrixToGzippedFile(ICed, os.path.join(ICedFol,key+".txt.gz"))  # @UndefinedVariable
                print "saved ICed"
            if smooth:
                smoothed = mirnylib.numutils.adaptiveSmoothing(ICed, 20, originalCounts = data, maxSmooth = 10)
                mirnylib.numutils_new.matrixToGzippedFile(smoothed, os.path.join(smoothedFol, key+".txt.gz"))  # @UndefinedVariable
                print "saved smoothed"
        elif type(data) == np.ndarray:            
            np.savetxt(os.path.join(outFolder,key), data)
        else:
            data_str = repr(data)
            open(os.path.join(outFolder,key),'w').write(data_str) 
            
            
def completeEig(data, GC = None,doSmooth = False ):
    """
    Performs eigenvector decomposition of cis data with all proper filters
    Parameters
    ----------
    data : 2D np.array 
        heatmap to analyze
    GC : 1D np.array
        vector of GC content to flip eigenvector (optional but recommended)
    doSmooth : bool (default=False)
        Perform gentle iterative smoothing
    """    
        
    data = np.clip(data, 0, np.percentile(data, 99.99))    
    mirnylib.numutils.removeDiagonals(data, 1)
    ICdata = mirnylib.numutils.completeIC(data, 20)
    mask = np.sum(ICdata, axis=0) == 0
    if (-mask).sum() < 3:
        return np.zeros(len(data), dtype = np.float)
        
    for dia in [-1, 0, 1]:
        mirnylib.numutils.fillDiagonal(ICdata, np.mean(np.diag(ICdata, 2)))
    ICdata[mask] = 0
    ICdata[:, mask] = 0
    if doSmooth:
        ICdata = mirnylib.numutils.adaptiveSmoothing(ICdata, 3, originalCounts=data, maxSmooth = 5)
    ICdata = mirnylib.numutils.observedOverExpected(ICdata)
    ICdata = mirnylib.numutils.iterativeCorrection(ICdata)[0]
    PCs = mirnylib.numutils.zeroEIG(ICdata)[0]
    PC1 = PCs[0]
    PC2 = PCs[1]
    if not GC is None:
        if np.abs(spearmanr(GC, PC2)[0]) > np.abs(spearmanr(GC, PC1)[0]):
            PC2, PC1 = PC1, PC2
        if spearmanr(GC, PC1)[0] < 0:
            PC1 = -PC1
    return PC1

def byChrEig(filename, genome, chromosomes = "all",resolution = "auto", byArm = True, doSmooth = False):    
    from mirnylib.genome import Genome        
    if resolution == "auto":
        resolution = getResolution(filename)            
    if type(genome) == str: 
        genome = Genome(genome)         
    assert isinstance(genome, Genome)
    genome.setResolution(resolution)
    mydict = mirnylib.h5dict.h5dict(filename)    
    if chromosomes == "all":
        chromosomes = range(genome.chrmCount)
        chromosomes = [i for i in chromosomes if "{0} {0}".format(i) in mydict]
        if len(chromosomes) == 0:
            raise ValueError("No chromosomes left. Check h5dict file.")
            
    result = []    
    for chrom in chromosomes:
        data = mydict["{0} {0}".format(chrom)]
        if not byArm:
            result.append(completeEig(data, genome.GCBin[chrom],doSmooth=doSmooth))
        else:
            GC = genome.GCBin[chrom]            
            result.append(np.zeros(len(GC), dtype = float))
            cent = genome.cntrMids[chrom] / resolution 
            result[-1][:cent] =  completeEig(data[:cent,:cent], genome.GCBin[chrom][:cent],doSmooth=doSmooth)
            result[-1][cent:] =  completeEig(data[cent:,cent:], genome.GCBin[chrom][cent:],doSmooth=doSmooth)
    return result 
                          
        
    
