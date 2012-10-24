import os.path
import sys

import numpy as np
import cPickle
#from polymerutils import save, load
from tempfile import NamedTemporaryFile
from polymerutils import save
from scipy.interpolate.interpolate import interp1d
from scipy.interpolate.fitpack2 import InterpolatedUnivariateSpline
from mirnylib.systemutils import deprecate





def interpolateData(data, targetN=90000, colorArrays=[]):
    """
    Converts a polymer of any length to a smoothed chain with (hopefully)
    fixed distance between neighboring monomers. Does it by cubic spline
    interpolation as following.

    1. Interpolate the data using cubic spline \n
    2. Evaluate cubic spline at targetN*10 values \n
    3. Rescale the evaluated spline such that total distance is targetN \n
    4. Select targetN points along the path with distance between
    neighboring points _along the chain_ equal to 1. 
    
    Parameters
    ----------
    data : Nx3 array
        Input xyz coordinates
    targetN : int
        Length of output polymer. 
        It is not adviced to make it many times less than N
    
    Returns
    -------
    (about targetN) x 3 array
    """

    fineGrain = 10

    N = len(data)
    numDim = len(data[0])
    targetDataSize = targetN * fineGrain

    evaluateRange = np.arange(N)
    targetRange = np.arange(0, N - 1, N / float(targetDataSize))

    splined = np.zeros((len(targetRange), numDim), float)
    colorsSplined = []
    for coor in xrange(numDim):
        spline = InterpolatedUnivariateSpline(evaluateRange,
                                        data[:, coor], k=3)
        evaled = spline(targetRange)
        splined[:, coor] = evaled

    for color in colorArrays:
        spline = InterpolatedUnivariateSpline(evaluateRange,
                                        color, k=1)
        evaled = spline(targetRange)
        colorsSplined.append(evaled)

    dists = np.sqrt(np.sum(np.diff(splined, 1, axis=0) ** 2, axis=1))
    totalDist = np.sum(dists)
    mult = totalDist / targetN
    splined /= mult
    dists /= mult
    cumDists = np.cumsum(dists)
    searched = np.searchsorted(cumDists, np.arange(1, targetN))
    v1 = cumDists[searched]
    v2 = cumDists[searched - 1]
    vals = np.floor(v1)
    p1 = (v1 - vals) / (v1 - v2)
    p2 = 1 - p1

    colorReturn = [i[searched] for i in colorsSplined]

    evaled = p2[:, None] * splined[searched] + \
    p1[:, None] * splined[searched - 1]
    return evaled, colorReturn
















def convert_xyz(data, out_file):
    "Converts an XYZ data to a fake pdb file"

    data = data - np.min(data, axis=0)[None, :]

    N = len(data)
    retret = ""
    for i, line in enumerate(data):
        def add(st, n):
            if len(st) > n: return st[:n]
            else:

                return st + " "*(n - len(st))


        line = [1. * (float(j) + 1) for j in line]
        ret = add("ATOM", 7)
        ret = add(ret + "%i" % (i + 1), 13)
        ret = add(ret + "CA", 17)
        ret = add(ret + "ALA", 22)

        ret = add(ret + "%i" % (i), 30)
        ret = add(ret + ("%8.3f" % line[0]), 37)
        ret = add(ret + ("%8.3f" % line[1]), 45)
        ret = add(ret + ("%8.3f" % line[2]), 53)
        ret = add(ret + (" 1.00"), 61)
        ret = add(ret + str(float(i % 8 > 4)), 67)
        retret += (ret + "\n")

    f = out_file
    f.write(retret)
    f.flush()
    #print retret



def create_regions(a):
    """
    Creates array of non-zero regions of a.
    if a is 0 1 1 0 1 0 
    result will be (1,3), (4,5), because elements 1,2 and 4 are non-zero.     
     
    """

    a = a > 0
    a = np.r_[False, a, False]
    a1 = np.nonzero(a[1:] * (1 - a[:-1]))[0]
    a2 = np.nonzero(a[:-1] * (1 - a[1:]))[0]

    return np.transpose(np.array([a1, a2]))


def do_coloring(data, regions, colors, transparencies,
                chain_radius=0.02, subchain_radius=0.04,
                chain_transparency=0.5, support="",
                multiplier=.4,
                spherePositions=[],
                sphereRadius=.3,
                force=False,
                misc_arguments=""):

    """
    !!! Please read this completely. Otherwise you'll suck :( !!!
    
    Creates a PDB file and a rasmol script that shows an XYZ polymer
    using pymol. Polymer consists of two parts: chain and subchain.
    A chain is a grey polymer, that is meant to resemble the main chain.
    It is meant to be thin, gray and transparent (overall, here transparency 1
    means transparent, transparency 0 means fully visible). A subchain
    consists of a certain number of regions, each has it's own color, trans
    parency, etc.
    
    Parameters
    ----------
    
    data : an Nx3 array of XYZ coordinates
    
    regions : a list of tuples (start, end)
        Note that rasmol acceps subchains in a format 
        (first monomer, last monomer), not the usual python 
        convention (first, last+1)!!! An overlap check will watch this. 
    
    colors : a list of colors ("red", "green", "blue", etc.)  for each region
    
    transparencies : a list of floats between 0 and 1. 0 is fully visible
    
    chain_radius : radius of a main chain in arbitraty units
    
    subchain_radius : radius of a subchain in arbitrary units
    
    chain_transparency : transparency of a main chain
    
    support : code to put at the end of the script
        put all the "save" or "ray" commands here if you want automation
    multiplier : a number, probably between .1 and 3
        Increasing it makes chain more smooth
        Decreasing it makes it more kinky, but may cause bugs
        or even missing chain regions
    misc_arguments : str
        Misc arguments to pymol command at the very end (mainly >/dev/null)
    
    .. warning :: Do not call this scripy "pymol.py!"
    
    .. warning ::
        Please resize the window to the needed size and run 
        "ray" command (press "ray" button) to get a nice image.
        Then DO NOT MOVE the image and find "export" in the menu.
        Otherwise your image will be not that high quality
        
    .. note ::
        performance of "ray" command depends on two things.
        First is resolution : it is more than quadratic in that
        Second is chain complexity. Tens thousand of monomers
        at high resolution may take up to an hour to ray.
        Though it actually looks awesome then!

    Run an example method below to see how the code works.
    See full automation examples below. 
    """
    data = np.array(data)
    data *= multiplier
    chain_radius *= multiplier
    subchain_radius *= multiplier
    sphereRadius *= multiplier

    #starting background check
    N = len(data)
    nregions = np.array(regions)
    if len(nregions) > 0:
        if nregions.min() < 0 or nregions.max() >= N:
            raise ValueError("region boundaries should be between 0 and N-1")
    covered = np.zeros(len(data), int)
    for i in regions:
        covered[i[0]:i[1] + 1] += 1
    if (covered.max() > 1) and (force == True):

        raise ValueError("Overlapped regions detected! Rasmol will not work"\
                         " Note that regions is (first,last), not last+1!")
    bgcolor = "grey"
    letters = [i for i in "1234567890abcdefghijklmnopqrstuvwxyz"]
    names = [i + j + k for i in letters for j in letters for k in letters]

    out = NamedTemporaryFile()
    pdbFile = NamedTemporaryFile()
    pdbname = os.path.split(pdbFile.name)[-1]

    out.write("hide all\n")
    out.write("bg white\n")

    for i in xrange(len(regions)):

        out.write("select %s, resi %d-%d\n" % (names[i], regions[i][0], regions[i][1]))
        out.write("create subchain%s,%s\n" % (names[i], names[i]))
        #out.write("remove subchain%s in %s\n"%(names[i],pdbname))

    out.write("set cartoon_trace_atoms,1,%s\n" % pdbname)
    out.write("cartoon tube,%s\n" % pdbname)
    out.write("set cartoon_tube_radius,%f,%s\n" % (chain_radius, pdbname))
    out.write("set cartoon_transparency,%f,%s\n" % (chain_transparency, pdbname))
    out.write("color %s,%s\n" % (bgcolor, pdbname))
    for i in xrange(len(regions)):

        name = "subchain%s" % names[i]
        out.write("set cartoon_trace_atoms,1,%s\n" % name)
        out.write("cartoon tube,%s\n" % name)
        out.write("set cartoon_tube_radius,%f,%s\n" % (subchain_radius, name))
        out.write("color %s,subchain%s\n" % (colors[i], names[i]))
        out.write("set cartoon_transparency,%f,%s\n" % (transparencies[i], name))
    for i  in spherePositions:
        out.write("show spheres, i. {0}-{0}\n".format(i))
        out.write("set sphere_color, grey60 \n")

    out.write("alter all, vdw={0} \n".format(sphereRadius))
    out.write("show cartoon,name ca\n")
    out.write("zoom %s" % pdbname)
    out.write(support)
    out.flush()

    #saving data

    convert_xyz(data, pdbFile)

    from time import sleep
    sleep(0.5)

    print os.system("pymol {1} -u {0} {2}".format(out.name, pdbFile.name,
                                                  misc_arguments))




def example_pymol():
    #Creating a random walk
    rw = .4 * np.cumsum(np.random.random((1000, 3)) - 0.5, axis=0)

    #Highlighting first 100 monomers and then 200-400
    regions = [(0, 100), (200, 400)]

    #Coloring them red and blue
    colors = ["red", "green"]

    #Making red semi-transparent
    transp = [0.7, 0]

    #Running the script with default chain radiuses
    do_coloring(
                data=rw,
                regions=regions,
                colors=colors,
                transparencies=transp,
                spherePositions=[500, 600])



def show_chain(data, chain_radius=0.3, dataMult=1, support="",
                spherePositions=[],
                sphereRadius=.3,
               ):
    """This was meant to show rainbow colored worms. 
    Not sure if it still works, but you can try
    """
    data *= dataMult
    data -= np.min(data, axis=0)[None, :]
    print data.min()

    #regions = [(10,20),(120,140),(180,250)]
    dataFile = NamedTemporaryFile()
    out = NamedTemporaryFile()
    convert_xyz(data, dataFile)
    bgcolor = "grey"
    pdbname = dataFile.name.split("/")[-1]
    out.write("hide all\n")
    out.write("bg white\n")

    out.write("hide all\n")
    out.write("set cartoon_trace_atoms,1,%s\n" % pdbname)
    out.write("cartoon tube,%s\n" % pdbname)
    out.write("set cartoon_tube_radius,%f,%s\n" % (chain_radius, pdbname))
    out.write("spectrum\n")
    out.write("show cartoon,name ca\n")
    out.write("zoom %s" % pdbname)
    for i  in spherePositions:
        out.write("show spheres, i. {0}-{0}\n".format(i))
        out.write("set sphere_color, grey60 \n")
    out.write(support)
    out.flush()

    os.system("pymol {0} -u {1}".format(dataFile.name, out.name))
    #out.close()



