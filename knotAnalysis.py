# (c) 2013 Massachusetts Institute of Technology. All Rights Reserved
# Code written by: Maksim Imakaev (imakaev@mit.edu)

import numpy
from scipy import weave
import os
import os.path
from tempfile import NamedTemporaryFile
from polymerutils import grow_rw, getLinkingNumber, create_random_walk
from polymerutils import findSimplifiedPolymer
import platform
import polymerutils
arch = platform.architecture()

folderName = os.path.split(__file__)[0]

if arch == "32bit":
    reduceKnotFilename = os.path.join(folderName, "Reduce_knot20_x86")
else:
    reduceKnotFilename = os.path.join(folderName, "Reduce_knot20")



def getKnotNumber(data, evalAt=-1.1):
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


        print "runnung command {0} {1} -p {4}  > {2}_{3}".format(reduceKnotFilename, name,
                                     name, "_output", evalAt)
        os.system("{0} {1} -p {4}  > {2}_{3}".format(reduceKnotFilename, name,
                                     name, "_output", evalAt))
        lines = open("%s_%s" % (name, "_output")).readlines()
        print "Contents of the output: -----"
        print lines
        print "End of the output-----"
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
    sim.setup(platform="cuda",integrator = "variableLangevin", errorTol = 0.01)
    sim.load(data)
    sim.randomizeData()
    if mode == "auto":
        if sim.dist(0, sim.N - 1) < 2:
            mode = "ring"
        else:
            mode = "chain"
            sim.tetherParticles([0, sim.N - 1], 5)
            # sim.addGravity()
    sim.setChains()
    sim.addHarmonicPolymerBonds(wiggleDist=0.06)
    sim.addGrosbergRepulsiveForce(trunc=60)
    sim.addGrosbergStiffness(k=3)
    # sim.localEnergyMinimization(tolerance = 0.001)
    # sim.localEnergyMinimization()
    sim.doBlock(40)
    for _ in xrange(steps):
        sim.doBlock(2000)
    data = sim.getData()
    del sim
    sleep(0.5)
    return data



def analyzeKnot(data, useOpenmm=False, simplify=True, evalAt=-1.1, lock=None, offset=0, stepMult=1):
    """
    Takes a polymer ring or chain, and analyzes knot number

    Parameters
    ----------

    data : (nx3) or (3xn) array
        Input data to analyze
    useOpenmm : bool
        If True, first simplify the polymer with OpenMM.
        Note that this may unfold chains a little bit.
    evalAt : float
        Evaluate A(x) * A(1/x). Use x=-1, or x=-1.1
    lock : lock object
        Lock object to prevent concurrent use of OpenMM
        Use this if you use multithreading
        I prefer multiprocessing.Pool.map, and multiprocessing.Manager.Lock()

    offset : int (optional)
        OpenMM is usually adjusted to kick in only for sufficiently complex knots.
        OpenMM would not start for knots which simplify to less than 250 monomers.
        However, sometimes knot calculation may take longer than expected.
        This is, for example, relevant for very knotted short polymers.
        Then you can set offset to a negative number to let OpenMM kick in earlier.
        Or set it to a positive number to let it kick in later.

        If offset is -100, then OpenMM will start working when a polymer is
        simplified only to 150 monomers, not to 250.

    stepMult : float (optional)
        Multiplies the number of steps which OpenMm will do.
        OpenMM is run for some time. If you want to reduce or increase this time,
        You can use this flag. The main purpose would be to trigger OpenMM early, but
        let it run for less.
        Setting stepMult to 0.25 will do four times less OpenMM.
        Setting it to 2 would make twice the amount of OpenMM


    """

    data = numpy.asarray(data)
    if len(data) == 3:
        data = data.T

    if simplify:
        t = findSimplifiedPolymer(data)
        if useOpenmm == True:
            if len(t) > 250 + offset:
                ll = len(t)
                if ll < 300 + offset:
                    steps = 2
                elif ll < 400 + offset:
                    steps = 4
                elif ll < 450 + offset:
                    steps = 10
                elif ll < 500 + offset:
                    steps = 15
                elif ll < 550 + offset:
                    steps = 25
                else:
                    steps = 30
                if lock != None:
                    lock.acquire()
                    data = expandPolymerRing(data, steps=int((steps - 1) * stepMult) + 1)
                    lock.release()
                else:
                    data = expandPolymerRing(data, steps=steps)
                t = findSimplifiedPolymer(data)
    else:
        t = data

    print "simplified from {0} to {1} monomers".format(len(data), len(t))

    try:
        print "OpenMM helped: %d to %d" % (ll, len(t))
    except:
        pass
    output = getKnotNumber(t, evalAt=evalAt)
    word = output[0].split()[1]
    if word == "0_0":
        return 1
    return float(word)


def _testAnalyzeKnot():
    np = numpy

    for _ in xrange(100):
        a = polymerutils.grow_rw(1000, 12, method="standard")
        kn = analyzeKnot(a, simplify=False)
        print kn
        assert kn == 1


    for _ in xrange(50):
        a = np.random.random((40, 3))
        ka = analyzeKnot(a, simplify=False)
        mat = np.random.random((3, 3))
        kb = analyzeKnot(np.dot(a, mat), simplify=False)
        print
        print ka, kb
        assert np.abs(ka / kb - 1) < 0.0001
        print

# _testAnalyzeKnot()

def _testSimplify():
    np = numpy

    for _ in xrange(200):
        s = 20
        # a = np.cumsum(np.random.randn(s, 3), axis=0) + np.random.randn(s, 3) * 2
        a = np.random.randn(s, 3) * 2

        ka = analyzeKnot(a, simplify=False)
        kb = analyzeKnot(a, simplify=True)
        print
        print ka, kb
        assert np.abs(ka / kb - 1) < 0.0001
        print


# _testSimplify()
