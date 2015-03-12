# (c) 2013 Massachusetts Institute of Technology. All Rights Reserved
# Code written by: Maksim Imakaev (imakaev@mit.edu)
"""
This file contains a bunch of method to work on contact maps of a Hi-C data.
It uses a lot of methods from mirnylib repository.




Find contacts of a conformation
-------------------------------

To be filled in later


Find average contact maps
-------------------------

"""


import subprocess
import numpy
from mirnylib.h5dict import h5dict
import os
import traceback
from mirnylib.numutils import sumByArray
import exceptions
np = numpy
from scipy import weave
from math import sqrt
from mirnylib.systemutils import fmapred, fmap, deprecate, setExceptionHook
import sys
import mirnylib.numutils
from polymerutils import load, save
import warnings
import polymerutils




try:
    import simtk.openmm
    openmm = True
except:
    openmm = False
    CPU = False

if openmm == True:
    CPU = True

    ver = simtk.openmm.__version__
    nums = tuple([int(i) for i in ver.split(".")])
    if nums < (6, 0):
        print "OpenMM CPU contactlist would only work with OpenMM >= 6.1"
        print 'Switching to default contact list finder'
        CPU = False
    if nums > (6, 2):
        print"You might need to recompile OpenMM neighbors list"
        print "please do this if something fails"
    if ((nums < (6, 1)) and (nums >= (6, 0))):
        CPUfile = "getCpuNeighborList6.0"
    else:
        CPUfile = "getCpuNeighborList6.1"
    folderName = os.path.split(__file__)[0]
    CPUfile = os.path.join(folderName, CPUfile)
    if not os.path.exists(CPUfile):
        raise ValueError("OPenMM contactmap file missing! This bug should be reported.")

try:
    simtk.openmm.Platform_getPlatformByName("CPU")
except:
    print "Not using OpenMM contact map finder"
    CPU = False



# a set of checks for OpenMM contact list, and selector of OpenMM contact list finder



def intload(filename, center="N/A"):
    """
    .. warning:: About to be deprecated
    """
    "dedicated weave.inline for loading int polymers"
    f = open(filename, 'r')
    line = f.readline()
    f.close()
    N = int(line)
    ret = numpy.zeros((3, N), order="C", dtype=numpy.int)
    code = """
    #line 85 "binary_search.py"
    using namespace std;
    FILE *in;
    const char * myfile = filename.c_str();
    in = fopen(myfile,"r");
    int dummy;
    dummy = fscanf(in,"%d",&dummy);
    for (int i = 0; i < N; i++)
    {
    dummy = fscanf(in,"%ld %ld %ld ",&ret[i],&ret[i + N],&ret[i + 2 * N]);
    if (dummy < 3){printf("Error in the file!!!");throw -1;}
    }
    fclose(in);
    """
    support = """
    #include <math.h>
    """
    weave.inline(code, ['filename', 'N', 'ret'], extra_compile_args=[
        '-march=native -malign-double'], support_code=support)

    return ret

class h5dictLoad(object):
    """
    An experimental class to fetch h5dict values based on a fake filename.

    It accepts filenames in a form
    /path-to-h5dict/h5dictKey

    It automatically caches h5dicts that are already open

    use simpleFetch method for concurrent fetching using fmap
    """
    def __init__(self):
        self.baseDict = {}

    def fetch(self, filename, dummy=True):
        base, num = os.path.split(filename)

        if base not in self.baseDict.keys():
            self.baseDict[base] = h5dict(base, mode='r')
        return self.baseDict[base][num]

    def simpleFetch(self, filename, dummy=True):
        base, num = os.path.split(filename)
        h5 = h5dict(base, mode="r")
        toret = h5[num]
        del h5
        return toret


def rad2(data):
    """
    Returns Rg(N^(2/3)
    This is an important point, as for a dense equilibrium globule at N^{2/3}
    both P(s) and Rg(s) plots reach a plateau.
    """

    if len(data) != 3:
        data = numpy.transpose(data)
    if len(data) != 3:
        raise ValueError("Wrong dimensions of data")

    def give_radius_scaling(data):
        N = len(data[0])
        target = int(N ** (2 / 3.))
        coms = numpy.cumsum(data, 1)
        coms2 = numpy.cumsum(data ** 2, 1)

        def radius_gyration(len2):
            coms2d = (-coms2[:, :-len2] + coms2[:, len2:]) / len2
            comsd = ((coms[:, :-len2] - coms[:, len2:]) / len2) ** 2
            diffs = coms2d - comsd
            sums = numpy.sqrt(numpy.sum(diffs, 0))
            return numpy.mean(sums)
        return radius_gyration(target)
    return give_radius_scaling(data)


def giveIntContacts(data):
    """give all contacts of a polymer on a cubic lattice
    Intersections are not counted as contacts. Sorry :(

    Parameters
    ----------
    data : Nx3 or 3xN array of ints
    """

    data = numpy.asarray(data, dtype=int)

    if len(data.shape) != 2:
        raise ValueError("Wrong dimensions of data")
    if 3 not in data.shape:
        raise ValueError("Wrong size of data: %s,%s" % data.shape)
    if data.shape[0] == 3:
        data = data.T

    data -= numpy.min(data, axis=0)[None, :]

    M = numpy.max(data) + 1
    if M > 1500:
        raise ValueError("Polymer is to big, can't create bounding box!")

    N = len(data)
    tocheck = numpy.zeros(M * M * M, dtype=numpy.int32) - 1
    tocheck[data[:, 0] + data[:, 1] * M + data[:, 2] * M *
            M] = numpy.arange(N, dtype=numpy.int32)
    tocheck.shape = (M, M, M)
    contacts1 = numpy.concatenate([tocheck[1:, :, :].ravel(), tocheck[
        :, 1:, :].ravel(), tocheck[:, :, 1:].ravel()])
    contacts2 = numpy.concatenate([tocheck[:-1, :, :].ravel(), tocheck[
        :, :-1, :].ravel(), tocheck[:, :, :-1].ravel()])
    mask = (contacts1 != -1) * (contacts2 != -1)
    contacts1 = contacts1[mask]
    contacts2 = contacts2[mask]
    contacts3 = numpy.minimum(contacts1, contacts2)
    contacts4 = numpy.maximum(contacts1, contacts2)
    return numpy.concatenate([contacts3[:, None], contacts4[:, None]], 1)


def giveContactsOpenMM(data, cutoff=1.7):
    """
    A wrapper to an ultra-fast openmm filter by VJ Pande
    If you use it, you may want to cite OpenMM
    """
    assert CPU == True
    assert os.path.exists(CPUfile)
    data = numpy.asarray(data)
    if len(data.shape) != 2:
        print "bad data"
        raise ValueError("Wrong dimensions of data")
    if 3 not in data.shape:
        print "bad data"
        raise ValueError("Wrong size of data: %s,%s" % data.shape)
    if data.shape[0] == 3:
        data = data.T

    data = np.asarray(data, order="C", dtype=np.float32)


    newProcess = subprocess.Popen([CPUfile, str(cutoff), str(len(data))], stdin=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=-1)

    try:
        output, err = newProcess.communicate(data.tostring(order="C"))
    except:
        print "Process failed!"
        return None

    returncode = newProcess.returncode
    if returncode != 0:
        print "Bad return code"
        return None

    array = np.fromstring(output, dtype=np.int32)
    array = array.reshape((-1, 2))

    return array


testData = np.random.random((10, 3))
try:
    giveContactsOpenMM(testData)
except:
    CPU = False
    print '---Failed to execute OpenMM contact finder-----'
    print "Error code below"
    traceback.print_exc(file=sys.stdout)
    print "---Continuing withont OpenMM contact finder---"
    print "You may try to recompile OpenMM contact finder to make it work"





def giveContactsAny(data, cutoff=1.7, maxContacts=300):
    """Returns contacts of any sets of molecules with a given cutoff.
    Is slower than give_contacts, but can tolerate multiple chains.

    Parameters
    ----------
    data : Nx3 or 3xN array
        Polymer configuration. One chaon only.
    cutoff : float , optional
        Cutoff distance that defines contact
    maxContacts : int
        Maximum number of contacts per monomer.
        If total number of contacts exceeds maxContacts*N, program becomes uncontrollable.

    Returns
    -------

    k by 2 array of contacts. Each row corresponds to a contact.
    """
    data = numpy.asarray(data)
    if len(data.shape) != 2:
        raise ValueError("Wrong dimensions of data")
    if 3 not in data.shape:
        raise ValueError("Wrong size of data: %s,%s" % data.shape)
    if data.shape[0] == 3:
        data = data.T
    data = numpy.asarray(data, float, order="C")
    N = len(data)
    points = numpy.zeros((maxContacts * N, 2), int, order="C")
    N  # Eclipse warning remover

    code = """
    #line 50 "binary_search.py"
        using namespace std;
        double CUTOFF2 = cutoff * cutoff;

        int counter  = 0;

        double  d;
        for (int i=0;i<N;i++)
            {
            for (int j=i+1;j<N;j++)
                {

                if (dist(data,N,i,j) <= CUTOFF2)
                    {
                    int t = j - i;
                    points[2*counter] = i;
                    points[2*counter + 1 ] = j;
                    counter ++ ;
                    //cout << i << " "   << j  << " "  << d <<  " " << counter << endl;
                    }
            }
            }
       """
    support = """
        #include <math.h>
        double sq(double x) {return x*x;}
        double  dist(double* data, int N, int i,int j)
        {
        return (pow((data[3 * i]-data[3 * j]),2.) + pow((data[3 * i + 1] - data[3 * j + 1]),2.) + pow((data[3 * i + 2] - data[3 * j + 2]),2.));
        }
       """
    #
    weave.inline(code, ['data', 'N', 'cutoff', 'points'], extra_compile_args=[
        '-march=native -malign-double'], support_code=support)
    k = numpy.max(numpy.nonzero(points[:, 1]))
    return points[:k + 1, :]





def giveContacts(data, cutoff=1.7, maxContacts=300, method="auto", tryOpenMM=True):
    """Returns contacts of a single polymer with a given cutoff

    .. warning:: Use this only to find contacts of a single polymer chain
    with distance between monomers of 1.
    Multiple chains will lead to silent bugs.

    Parameters
    ----------
    data : Nx3 or 3xN array
        Polymer configuration. One chaon only.
    cutoff : float , optional
        Cutoff distance that defines contact
    maxContacts : int
        Maximum number of contacts per monomer.
        If total number of contacts exceeds maxContacts*N,
        program becomes uncontrollable.
    method : "auto", "continuous", "any"
        Selects between a method suitable for a continuous polymer, and an auto method

    Returns
    -------

    k by 2 array of contacts. Each row corresponds to a contact.
    """

    data = numpy.asarray(data)
    if np.isnan(data).any():
        raise RuntimeError("Data contains NANs")

    if max(data.shape) < 2000:
        return giveContactsAny(data, cutoff, maxContacts)

    if len(data.shape) != 2:
        raise ValueError("Wrong dimensions of data")
    if 3 not in data.shape:
        raise ValueError("Wrong size of data: %s,%s" % data.shape)
    if data.shape[0] == 3:
        data = data.T


    if (CPU == True) and (tryOpenMM == True) and (len(data) > 2000):

        conts = giveContactsOpenMM(data, cutoff)
        if conts == None:
            warnings.warn(RuntimeWarning("OpenMM did not work. Please use OpenMM 6.1 or recompile contact searcher. "))
        else:
            return conts


    dists2 = numpy.sqrt(numpy.sum(numpy.diff(data, axis=0) ** 2, axis=1))
    maxRatio = dists2.max() / numpy.median(dists2)

    if method == "any":
        return giveContactsAny(data, cutoff, maxContacts)
    if (method == "auto") and (maxRatio > 5):
        warnings.warn("\nPolymer does not seem continuous, falling back"
                      "to arbitrary contact finger 'give_contacts_any'"
                      "\n This is ok, just be aware of this! ")
        print "ratio of smaller to larger bonds is {0}".format(maxRatio)
        return giveContactsAny(data, cutoff, maxContacts)
        safeDistance = numpy.percentile(dists2, 99)
        cutoff = cutoff / safeDistance
        data = data / safeDistance

    data = numpy.asarray(data, float, order="C")
    cutoff = float(cutoff)
    N = len(data)
    points = numpy.zeros((maxContacts * N, 2), int, order="C")
    N  # Eclipse warning remover
    code = r"""
    #line 196 "binary_search.py"
        using namespace std;
        int counter  = 0;
          double  d;
        for (int i=0;i<N;i++)
            {
            for (int j=i+1;j<N;j++)
                {
                d = dist(data,N,i,j);

                if (d<= cutoff)
                    {
                    int t = j - i;
                    points[2*counter] = i;
                    points[2*counter + 1 ] = j;
                    counter ++ ;
                    //cout << i << " "   << j  << " "  << d;
                    //cout <<  " " << counter << endl;

                    }
                else if (d>4)
                    {
                    j +=  d-4;
                    }
            }
            }
       """
    support = """
        #include <math.h>
        double  dist(double* data, int N, int i,int j)
        {
        return sqrt(pow((data[3 * i]-data[3 * j]),2.) + pow((data[3 * i + 1] - data[3 * j + 1]),2.) + pow((data[3 * i + 2] - data[3 * j + 2]),2.));
        }
       """
    weave.inline(code, ['data', 'N', 'cutoff', 'points'], extra_compile_args=[
        '-march=native -malign-double -O3'], support_code=support)

    k = numpy.max(numpy.nonzero(points[:, 1]))
    return points[:k + 1, :]


def giveDistanceMap(data, size=1000):
    """returns  a distance map of a polymer with a given size"""
    toret = numpy.zeros((size, size), float, order="C")

    data = numpy.asarray(data)
    if len(data.shape) != 2:
        raise ValueError("Wrong dimensions of data")
    if 3 not in data.shape:
        raise ValueError("Wrong size of data: %s,%s" % data.shape)
    if data.shape[0] == 3:
        data = data.T
    data = numpy.asarray(data, float, order="C")

    N = len(data)
    N  # Eclipse warning remover
    code = """
    #line 50 "binary_search.py"
        using namespace std;
        int counter  = 0;
          double  d;
        for (int i=0;i<size;i++)
            {
            for (int j=i;j<size;j++)
                {
                d = dist(data,N,N*i/size,N*j/size);
                toret[size*i+j] = d;
                toret[size*j+i] = d;
            }
            }
       """
    support = """
        #include <math.h>
        double  dist(double* data, int N, int i,int j)
        {
        return sqrt(pow((data[3 * i]-data[3 * j]),2.) + pow((data[3 * i + 1] - data[3 * j + 1]),2.) + pow((data[3 * i + 2] - data[3 * j + 2]),2.));
        }
       """
    weave.inline(code, ['data', 'N', 'size', 'toret'], extra_compile_args=[
        '-march=native -malign-double'], support_code=support)
    return toret




def rescalePoints(points, bins):
    "converts array of contacts to the reduced resolution contact map"
    a = numpy.histogram2d(points[:, 0], points[:, 1], bins)[0]
    a = a + numpy.transpose(a)

    return a


def rescaledMap(data, bins, cutoff=1.7, contactMap=None):
    # print data.sum(), bins.sum(), cutoff
    """calculates a rescaled contact map of a structure
    Parameters
    ----------
    data : Nx3 or 3xN array
        polymer conformation

    bins : Lx1 array
        bin starts

    cutoff : float, optional
        cutoff for contacts

    Returns
    -------
        resXres array with the contact map
    """

    t = giveContacts(data, cutoff)
    x = np.searchsorted(bins, t[:, 0]) - 1
    y = np.searchsorted(bins, t[:, 1]) - 1
    assert x.min() >= 0
    assert y.min() >= 0
    assert x.max() < len(bins) - 1
    assert y.max() < len(bins) - 1
    matrixSize = len(bins) - 1
    index = matrixSize * x + y
    unique = np.unique(index)
    inds = sumByArray(index, unique)
    uniquex = unique / matrixSize
    uniquey = unique % matrixSize

    if contactMap is None:
        contactMap = np.zeros((matrixSize, matrixSize), dtype=int)
    contactMap[uniquex, uniquey] += inds

    return contactMap








def pureMap(data, cutoff=1.7, contactMap=None):
    """calculates an all-by-all contact map of a single polymer chain.
    Doesn't work for multi-chain polymers!!!
    If contact map is supplied, it just updates it

    Parameters
    ----------
    data : Nx3 or 3xN array
        polymer conformation
    cutoff : float
        cutoff for contacts
    contactMap : NxN array, optional
        contact map to update, if averaging is used
    """
    data = numpy.asarray(data)
    if len(data.shape) != 2:
        raise ValueError("Wrong dimensions of data")
    if 3 not in data.shape:
        raise ValueError("Wrong size of data: %s,%s" % data.shape)
    if data.shape[0] == 3:
        data = data.T
    data = numpy.asarray(data, float, order="C")

    t = giveContacts(data, cutoff)
    N = data.shape[0]
    if contactMap is None:
        contactMap = numpy.zeros((N, N), "int32")
    contactMap[t[:, 0], t[:, 1]] += 1
    contactMap[t[:, 1], t[:, 0]] += 1
    return contactMap


observedOverExpected = deprecate(mirnylib.numutils.observedOverExpected)


# print observedOverExpected(numpy.arange(100).reshape((10,10)))

def cool_trunk(data):
    "somehow trunkates the globule so that the ends are on the surface"
    CUTOFF = 0.7
    CUTOFF2 = 3
    datax, datay, dataz = data[0], data[1], data[2]
    N = len(data[0])
    found1, found2 = 0, 0

    def dist(i):
        return sqrt(datax[i] ** 2 + datay[i] ** 2 + dataz[i] ** 2)

    def sqdist(i, j):
        return sqrt((datax[i] - datax[j]) ** 2 + (datay[i] - datay[j]) ** 2
                    + (dataz[i] - dataz[j]) ** 2)

    def sqdist2(i, j, scale):
        return sqrt((scale * datax[i] - datax[j]) ** 2 +
                    (scale * datay[i] - datay[j]) ** 2 +
                    (scale * dataz[i] - dataz[j]) ** 2)
    breakflag = 0
    escapeflag = 0
    for i in xrange(N):
        exitflag = 0
        pace = 0.25 / dist(i)
        for j in numpy.arange(1, 10, pace):
            # print j
            escapeflag = 1
            for k in xrange(N):
                if k == i:
                    continue
                escapeflag, breakflag

                if 0.001 < sqdist2(i, k, j) < CUTOFF:
                    breakflag = 1
                    print "breaking at", k
                    break
                if 0.001 < sqdist2(i, k, j) < CUTOFF2:
                    escapeflag = 0

            if breakflag == 1:
                breakflag = 0
                print i, dist(i), j
                break
            if escapeflag == 1:
                print i, dist(i)
                found1 = i
                exitflag
                exitflag = 1
                break
        if exitflag == 1:
            break
    for i in xrange(N - 1, 0, -1):
        exitflag = 0
        pace = 0.25 / dist(i)
        for j in numpy.arange(1, 10, pace):
            # print j
            escapeflag = 1
            for k in xrange(N - 1, 0, -1):
                if k == i:
                    continue
                escapeflag, breakflag

                if 0.001 < sqdist2(i, k, j) < CUTOFF:
                    breakflag = 1
                    print "breaking at", k
                    break
                if 0.001 < sqdist2(i, k, j) < CUTOFF2:
                    escapeflag = 0

            if breakflag == 1:
                breakflag = 0
                print i, dist(i), j
                break
            if escapeflag == 1:
                print i, dist(i)
                found2 = i
                exitflag
                exitflag = 1
                break
        if exitflag == 1:
            break

    f1 = found1
    f2 = found2
    datax[f1 - 1], datay[f1 - 1], dataz[f1 - 1] = datax[f1] * 2, datay[
        f1] * 2, dataz[f1] * 2
    datax[f2 + 1], datay[f2 + 1], dataz[f2 + 1] = datax[f2] * 2, datay[
        f2] * 2, dataz[f2] * 2
    a = (datax[f1 - 1], datay[f1 - 1], dataz[f1 - 1])
    b = (datax[f2 + 1], datay[f2 + 1], dataz[f2 + 1])
    c = (a[0] + b[0], a[1] + b[1], a[2] + b[2])
    absc = sqrt(c[0] ** 2 + c[1] ** 2 + c[2] ** 2)
    c = (c[0] * 40 / absc, c[1] * 40 / absc, c[2] * 40 / absc)
    datax[f1 - 2], datay[f1 - 2], dataz[f1 - 2] = c
    datax[f2 + 2], datay[f2 + 2], dataz[f2 + 2] = c
    return [datax[f1 - 2:f2 + 3], datay[f1 - 2:f2 + 3], dataz[f1 - 2:f2 + 3]]


def averageContactMap(*args, **kvargs):
    print 'please use averageBinnedContactMap or averagePureContactMap'
    raise Exception('deprecated function')


def averageBinnedContactMapOld(filenames, chains=None, binSize=None, cutoff=1.7,
                            n=4,  # Num threads
                            loadFunction=load,
                            exceptionsToIgnore=None):
    """
    Returns an average contact map of a set of conformations.
    Non-existing files are ignored if exceptionsToIgnore is set to IOError.
    example:\n

    An example:

    .. code-block:: python
        >>> filenames = ["myfolder/blockd%d.dat" % i for i in xrange(1000)]
        >>> cmap = averageBinnedContactMap(filenames) + 1  #getting cmap
        #either showing a log of a map (+1 for zeros)
        >>> plt.imshow(numpy.log(cmap +1))
        #or truncating a map
        >>> vmax = np.percentile(cmap, 99.9)
        >>> plt.imshow(cmap, vmax=vmax)
        >>> plt.show()

    Parameters
    ----------
    filenames : list of strings
        Filenames to average map over
    chains : list of tuples or Nx2 array
        (start,end+1) of each chain
    binSize : int
        size of each bin in monomers
    cutoff : float, optional
        Cutoff to calculate contacts
    n : int, optional
        Number of threads to use.
        By default 4 to minimize RAM consumption.
    exceptionsToIgnore : list of Exceptions
        List of exceptions to ignore when finding the contact map.
        Put IOError there if you want it to ignore missing files.

    Returns
    -------
    tuple of two values:
    (i) MxM numpy array with the conntact map binned to binSize resolution.
    (ii) chromosomeStarts a list of start sites for binned map.

    """

    getResolution = 0
    fileInd = 0
    while getResolution == 0:
        try:
            data = loadFunction(filenames[fileInd])  # load filename
            getResolution = 1
        except:
            fileInd = fileInd + 1
        if fileInd >= len(filenames):
            print "no valid files found in filenames"
            raise ValueError("no valid files found in filenames")

    if chains is None:
        chains = [[0, len(data)]]
    if binSize is None:
        binSize = int(numpy.floor(len(data) / 500))

    bins = []
    chains = numpy.asarray(chains)
    chainBinNums = (
        numpy.ceil((chains[:, 1] - chains[:, 0]) / (0.0 + binSize)))
    for i in xrange(len(chainBinNums)):
        bins.append(binSize * (numpy.arange(int(chainBinNums[i])))
                    + chains[i, 0])
    bins.append(numpy.array([chains[-1, 1] + 1]))
    bins = numpy.concatenate(bins)
    bins = bins - .5
    Nbase = len(bins) - 1

    if Nbase > 10000:
        warnings.warn(UserWarning('very large contact map'
                                  ' may be difficult to visualize'))

    chromosomeStarts = numpy.cumsum(chainBinNums)
    chromosomeStarts = numpy.hstack((0, chromosomeStarts))

    def action(i):  # Fetch rescale map from filename.
        print i
        try:
            data = loadFunction(i)  # load filename

        except exceptionsToIgnore:
            # if file faled to load, return empty map
            print "file not found"
            return numpy.zeros((Nbase, Nbase), "float")
        value = rescaledMap(data, bins, cutoff=cutoff)
        return value

    return fmapred(action, filenames, n=n, exceptionList=exceptionsToIgnore), chromosomeStarts[0:-1]


def averageBinnedContactMap(filenames, chains=None, binSize=None, cutoff=1.7,
                            n=4,  # Num threads
                            loadFunction=load,
                            exceptionsToIgnore=None):
    """
    Returns an average contact map of a set of conformations.
    Non-existing files are ignored if exceptionsToIgnore is set to IOError.
    example:\n

    An example:

    .. code-block:: python
        >>> filenames = ["myfolder/blockd%d.dat" % i for i in xrange(1000)]
        >>> cmap = averageBinnedContactMap(filenames) + 1  #getting cmap
        #either showing a log of a map (+1 for zeros)
        >>> plt.imshow(numpy.log(cmap +1))
        #or truncating a map
        >>> vmax = np.percentile(cmap, 99.9)
        >>> plt.imshow(cmap, vmax=vmax)
        >>> plt.show()

    Parameters
    ----------
    filenames : list of strings
        Filenames to average map over
    chains : list of tuples or Nx2 array
        (start,end+1) of each chain
    binSize : int
        size of each bin in monomers
    cutoff : float, optional
        Cutoff to calculate contacts
    n : int, optional
        Number of threads to use.
        By default 4 to minimize RAM consumption.
    exceptionsToIgnore : list of Exceptions
        List of exceptions to ignore when finding the contact map.
        Put IOError there if you want it to ignore missing files.

    Returns
    -------
    tuple of two values:
    (i) MxM numpy array with the conntact map binned to binSize resolution.
    (ii) chromosomeStarts a list of start sites for binned map.

    """
    n = min(n, len(filenames))
    subvalues = [filenames[i::n] for i in xrange(n)]

    getResolution = 0
    fileInd = 0
    while getResolution == 0:
        try:
            data = loadFunction(filenames[fileInd])  # load filename
            getResolution = 1
        except:
            fileInd = fileInd + 1
        if fileInd >= len(filenames):
            print "no valid files found in filenames"
            raise ValueError("no valid files found in filenames")

    if chains is None:
        chains = [[0, len(data)]]
    if binSize is None:
        binSize = int(numpy.floor(len(data) / 500))

    bins = []
    chains = numpy.asarray(chains)
    chainBinNums = (
        numpy.ceil((chains[:, 1] - chains[:, 0]) / (0.0 + binSize)))
    for i in xrange(len(chainBinNums)):
        bins.append(binSize * (numpy.arange(int(chainBinNums[i])))
                    + chains[i, 0])
    bins.append(numpy.array([chains[-1, 1] + 1]))
    bins = numpy.concatenate(bins)
    bins = bins - .5
    Nbase = len(bins) - 1

    if Nbase > 10000:
        warnings.warn(UserWarning('very large contact map'
                                  ' may be difficult to visualize'))

    chromosomeStarts = numpy.cumsum(chainBinNums)
    chromosomeStarts = numpy.hstack((0, chromosomeStarts))

    def myaction(values):  # our worker receives some filenames
        mysum = None  # future contact map.
        for i in values:
            try:
                data = loadFunction(i)
                print i
            except tuple(exceptionsToIgnore):
                print "file not found", i
                continue

            if data.shape[0] == 3:
                data = data.T
            if mysum is None:  # if it's the first filename,


                mysum = rescaledMap(data, bins, cutoff)  # create a map

            else:  # if not
                rescaledMap(data, bins, cutoff, mysum)
                    # use existing map and fill in contacts

        return mysum
    blocks = fmap(myaction, subvalues)
    blocks = [i for i in blocks if i is not None]
    a = blocks[0]
    for i in blocks[1:]:
        a = a + i
    a = a + a.T

    return a, chromosomeStarts


def averagePureContactMap(filenames,
                          cutoff=1.7,
                          n=4,  # Num threads
                          loadFunction=load,
                          exceptionsToIgnore=None):
    """
        Parameters
    ----------
    filenames : list of strings
        Filenames to average map over
    cutoff : float, optional
        Cutoff to calculate contacts
    n : int, optional
        Number of threads to use.
        By default 4 to minimize RAM consumption with pure maps.
    exceptionsToIgnore : list of Exceptions
        List of exceptions to ignore when finding the contact map.
        Put IOError there if you want it to ignore missing files.

    Returns
    -------

    An NxN (for pure map) numpy array with the contact map.
    """

    """
    Now we actually need to modify our contact map by adding
    contacts from each new file to the contact map.
    We do it this way because our contact map is huge (maybe a gigabyte!),
    so we can't just add many gigabyte-sized arrays together.
    Instead of this each worker creates an empty "average contact map",
    and then loads files one by one and adds contacts from each file to a contact map.
    Maps from different workers are then added together manually.
    """

    n = min(n, len(filenames))
    subvalues = [filenames[i::n] for i in xrange(n)]

    def myaction(values):  # our worker receives some filenames
        mysum = None  # future contact map.
        for i in values:
            try:
                data = loadFunction(i)
                print i
            except tuple(exceptionsToIgnore):
                print "file not found", i
                continue
            except:
                print "Unexpected error:", sys.exc_info()[0]
                print "File is: ", i
                return -1

            if data.shape[0] == 3:
                data = data.T
            if mysum is None:  # if it's the first filename,

                if len(data) > 6000:
                    warnings.warn(UserWarning('very large contact map'
                                              ' may cause errors. these may be fixed with n=1 threads.'))
                if len(data) > 20000:
                    warnings.warn(UserWarning('very large contact map'
                                              ' may be difficult to visualize.'))

                mysum = pureMap(data, cutoff)  # create a map

            else:  # if not
                pureMap(data, cutoff, mysum)
                    # use existing map and fill in contacts

        return mysum

    blocks = fmap(myaction, subvalues)
    blocks = [i for i in blocks if i is not None]
    a = blocks[0]
    for i in blocks[1:]:
        a = a + i
    return a


def _test():
    setExceptionHook()
    print "----> Testing giveContacts subroutine"
    try:
        c = polymerutils.load("/net/evolution/home/magus/trajectories/"
                              "globule_creation/32000_RW/crumpled1.dat")
        c = c[:16000]
    except:
        N = 16000
        a = numpy.random.random((N, 3))
        b = numpy.sqrt(numpy.sum(a ** 2, axis=1))
        a = a / b[:, None]
        b = (numpy.random.randn(N) * 0.1 + 1)
        a *= b[:, None]

        c = numpy.cumsum(a, axis=0)
        c = numpy.asarray(c)

    def transformContacts(contacts):
        a = np.sort(contacts, axis=0)
        a = [tuple(i) for i in a]
        return set(a)

    c2 = giveContactsAny(c, cutoff=2.2)
    c2  # Simply initialize weave.inline first

    from time import time
    a = time()
    c2 = giveContacts(c, cutoff=2.2, tryOpenMM=False)
    print "time for giveContacts is: ",
    print time() - a

    if CPU == True:
        a = time()
        c3 = giveContactsOpenMM(c, cutoff=2.2)
        print "time for giveContacts with openmm is: ",
        print time() - a
    else:
        c3 = c2

    a = time()
    c1 = giveContactsAny(c, cutoff=2.2)
    print "time for contactsAny is:", time() - a
    assert transformContacts(c1) == transformContacts(c2)
    assert transformContacts(c1) == transformContacts(c3)

    print "Test completed successfully"

    print "----> Testing averageContactMap subroutine"

    def funnyLoad(x):
        if x % 2 == 1:
            raise IOError("Blah Error")
        else:
            return c

    onemap = rescaledMap(c, -0.5 + np.arange(0, len(c) + 10, 10), cutoff=1.7)

    manyMap = averageBinnedContactMap(range(10),
                                      binSize=10,
                                      n=1,
                                      cutoff=1.7,
                                      loadFunction=funnyLoad,
                                      exceptionsToIgnore=[IOError, exceptions.IOError])
    manyMap = manyMap[0]
    fmapMap = averageBinnedContactMap(range(50),
                                      binSize=10,
                                      n=8,
                                      cutoff=1.7,
                                      loadFunction=funnyLoad,
                                      exceptionsToIgnore=[IOError, exceptions.IOError])[0]

    assert  np.abs(manyMap - 5 * onemap).sum() < 1
    assert  np.abs(fmapMap - 25 * onemap).sum() < 1
    print "Finished"

    print "----> Testing concurrent h5dict load"

    a = h5dict("bla", "w")
    for i in xrange(40):
        a[str(i)] = c + numpy.random.random(c.shape) * 0.2
    del a
    loader = h5dictLoad()
    robustLoad = loader.simpleFetch
    filenames = ["bla/%d" % i for i in xrange(50)]
    oneThreadMap = averageBinnedContactMap(filenames,
                                           binSize=10,
                                           n=1,
                                           loadFunction=robustLoad,
                                           exceptionsToIgnore=[IOError, KeyError])[0]
    eightThreadMap = averageBinnedContactMap(filenames,
                                             binSize=10,
                                             n=8,
                                             loadFunction=robustLoad,
                                             exceptionsToIgnore=[IOError, KeyError])[0]

    assert np.abs(oneThreadMap - eightThreadMap).sum() < 1
    print "Concurrent h5dict load test successful!"
    os.remove("bla")



if __name__ == "__main__":
    _test()
