#(c) 2013 Massachusetts Institute of Technology. All Rights Reserved
# Code written by: Maksim Imakaev (imakaev@mit.edu)

import numpy
from scipy import weave
import os
import os.path
from tempfile import NamedTemporaryFile
from polymerutils import grow_rw, getLinkingNumber, create_random_walk
from polymerutils import findSimplifiedPolymer

folderName = os.path.split(__file__)[0]
reduceKnotFilename = os.path.join(folderName, "Reduce_knot20")



def getKnotNumber(data, evalAt= -1.1):
    """A wrapper to code which gets knotting number of a polymer

    Parameters
    ----------
    data : an (nx3) or (3xn) array of points

    evalAt : float
        A number to evaluate Alexanders(x) * Alexandrs(1/x)
        You can use -1 or -1.1 - these are classics
    """
    data = numpy.array(data)
    if len(data) == 3:
        data = data.T

    with  NamedTemporaryFile() as newfile:
        newfile.write("t=0\n\n%d\n" % len(data))
        for j, i in enumerate(data):
            newfile.write("%d %lf %lf %lf\n" % tuple([j + 1] + list(i)))

        name = newfile.name
        newfile.flush()

        os.system("{0} {1} -p {4}  > {2}_{3}".format(reduceKnotFilename, name,
                                     name, "_output", evalAt))
        lines = open("%s_%s" % (name, "_output")).readlines()
        os.remove("%s_%s" % (name, "_output"))
        return lines


def expandPolymerRing(data, mode="auto", steps=20):
    """
    Expands polymer ring or chain using OpenMM.

    Parameters
    ----------
    data : Nx3 or 3xN array of coordinates
        Input coordinates of the polymer
    mode : str, optional
        "ring", or "chain", default - autodetect
    """

    from openmmlib import Simulation
    from time import sleep
    sim = Simulation(
        timestep=70, thermostat=0.002, velocityReinitialize=True)
    sim.setup(platform="cuda")
    sim.load(data)
    sim.randomizeData()
    if mode == "auto":
        if sim.dist(0, sim.N - 1) < 2:
            mode = "ring"
        else:
            mode = "chain"
    sim.setLayout(mode=mode)
    sim.addHarmonicPolymerBonds(wiggleDist=0.06)
    sim.addGrosbergRepulsiveForce(trunc=60)
    sim.addGrosbergStiffness(k=3)
    #sim.energyMinimization(stepsPerIteration=50)
    #sim.localEnergyMinimization()
    sim.doBlock(40)
    for _ in xrange(steps):
        sim.doBlock(2000)
    data = sim.getData()
    del sim
    sleep(0.5)
    return data



def analyzeKnot(data, useOpenmm=False, evalAt= -1.1, lock=None):
    """
    Takes a polymer ring or chain, and analyzes knot number

    Parameters
    ----------

    data : (nx3) or (3xn) array
        Input data to analyze
    useOpenmm : bool
        If True, first simplify the polymer
    evalAt : float
        Evaluate A(x) * A(1/x). Use x=-1, or x=-1.1
    lock : lock object
        Lock object to prevent concurrent use of OpenMM
        Use this if you use multithreading
        I prefer multiprocessing.Pool.map, and multiprocessing.Manager.Lock()

    """

    data = numpy.asarray(data)
    if len(data) == 3:
        data = data.T

    t = findSimplifiedPolymer(data)
    if useOpenmm == True:
        if len(t) > 250:
            ll = len(t)
            if ll < 300:
                steps = 2
            elif ll < 400:
                steps = 4
            elif ll < 450:
                steps = 10
            elif ll < 500:
                steps = 15
            elif ll < 550:
                steps = 25
            elif ll < 600:
                steps = 35
            elif ll < 800:
                steps = 60
            elif ll < 1200:
                steps = 100
            else:
                steps = 150
            if lock != None:
                lock.acquire()
                data = expandPolymerRing(data, steps=steps)
                lock.release()
            else:
                data = expandPolymerRing(data, steps=steps)
            t = findSimplifiedPolymer(data)



    #t = data
    print "simplified from {0} to {1} monomers".format(len(data), len(t))

    try:
        print "OpenMM helped: %d to %d" % (ll, len(t))
    except:
        pass
    number = getKnotNumber(t, evalAt=evalAt)
    print number
    num = float(number[0].split()[1])
    return num


