import os.path
import sys

import numpy as np
import cPickle
#from polymerutils import save, load
from tempfile import NamedTemporaryFile
from polymerutils import save

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
    "creates array of uniform regions"

    a = np.array(a, int)
    a = np.concatenate([np.array([0], int), a, np.array([0], int)])
    a1 = np.nonzero(a[1:] * (1 - a[:-1]))[0]
    a2 = np.nonzero(a[:-1] * (1 - a[1:]))[0]

    return np.transpose(np.array([a1, a2]))


def do_coloring(data, regions, colors, transparencies,
                chain_radius=0.02, subchain_radius=0.04,
                chain_transparency=0.5, support="",
                multiplier=.4,
                spherePositions=[],
                sphereRadius=.3):

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
    if nregions.min() < 0 or nregions.max() >= N:
        raise ValueError("region boundaries should be between 0 and N-1")
    covered = np.zeros(len(data), int)
    for i in regions:
        covered[i[0]:i[1] + 1] += 1
    if covered.max() > 1:
        raise ValueError("Overlapped regions detected! Rasmol will not work"\
                         " Note that regions is (first,last), not last+1!")
    bgcolor = "grey"
    letters = [i for i in "abcdefghigklmnopqrstuvwxyz"]
    names = [i + j for i in letters for j in letters]

    out = NamedTemporaryFile()
    pdbFile = NamedTemporaryFile()
    pdbname = os.path.split(pdbFile.name)[1]

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

    print os.system("pymol {1} -u {0}".format(out.name, pdbFile.name))




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
example_pymol()


def show_chain(data, chain_radius=0.003, dataMult=1, support=""):
    """This was meant to show rainbow colored worms. 
    Not sure if it still works, but you can try
    """
    data *= dataMult

    #regions = [(10,20),(120,140),(180,250)]
    dataFile = NamedTemporaryFile()
    out = NamedTemporaryFile()
    convert_xyz(data, dataFile)
    bgcolor = "grey"
    pdbname = "1pdb"
    out.write("hide all\n")
    out.write("bg white\n")

    out.write("hide all\n")
    out.write("set cartoon_trace_atoms,1,%s\n" % pdbname)
    out.write("cartoon tube,%s\n" % pdbname)
    out.write("set cartoon_tube_radius,%f,%s\n" % (chain_radius, pdbname))
    out.write("spectrum\n")
    out.write("show cartoon,name ca\n")
    out.write("zoom %s" % pdbname)
    out.write(support)
    out.flush()

    os.system("pymol {0} -u {1}".format(dataFile.name, out.name))
    #out.close()



