import os
import numpy as np
from openmmlib import SimulationWithCrosslinks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mirnylib.numutils import logbins, create_regions, continuousRegions
import polymerutils
from mirnylib.plotting import mat_img
from polymerutils import create_spiral, load
import cPickle
from mirnylib.systemutils import setExceptionHook
from brushManaging import makeBondsForBrush
import pymol_show
setExceptionHook()
import sys

if len(sys.argv) < 2:
    print "Set GPU number as first command line argument"
    print "Default value of 0 used"
    sys.argv = ["dummy", "0"]




def exampleOpenmm(supercoilCompaction=4, plectonemeLength=10, plectonemeGap=0, stifness=2, ind=0, domainNum=30, domainLen=5):
    """
    A method that performs one simulation.
    """

    a = SimulationWithCrosslinks(timestep=130, thermostat=0.01)

    a.setup(platform="cuda", verbose=True, GPU=sys.argv[1])

    foldername = "newSweep_lambda{0}_L{1}_gap{2}_stiff{3}_ind{4}".format(supercoilCompaction, plectonemeLength,
                                                         plectonemeGap, stifness, ind)
    a.saveFolder(foldername)  # folder where to save trajectory


    """
    Defining all parameters of the simulation
    """
    genomeLengthBp = 4000000
    lengthNm = 1500
    diamNm = 450
    supercoilBpPerNm = 4.5 * supercoilCompaction

    supercoilDiameterNm = 3.45 * np.sqrt(supercoilBpPerNm)
    supercoilBpPerBall = supercoilBpPerNm * supercoilDiameterNm * 2  #2 strands
    print supercoilBpPerBall
    print supercoilDiameterNm

    ballNumber = 2 * int(0.5 * (genomeLengthBp / supercoilBpPerBall))
    radiusMon = 0.5 * (diamNm / float(supercoilDiameterNm))
    halfLengthMon = 0.5 * (lengthNm / float(supercoilDiameterNm))

    N = ballNumber
    print N
    avLength = plectonemeLength

    #Growing a circular polymer
    from polymerutils import grow_rw
    bacteria = grow_rw(N, int(halfLengthMon * 2), "line")
    bacteria = bacteria - np.mean(bacteria, axis=0)
    a.load(bacteria)  # filename to load
    print bacteria[0]

    a.tetherParticles([0, N - 1])

    a.addCylindricalConfinement(r=radiusMon, bottom= -halfLengthMon, top=halfLengthMon, k=1.5)


    BD = makeBondsForBrush(chainLength=a.N)

    #Coordinates of top highly expressed genes
    geneCoord = [1162773, 3509071, 1180887, 543099, 1953250, 2522439, 3328524, 1503879, 900483, 242693, 3677144, 3931680, 3677704, 3762707, 3480870, 3829656, 1424678, 901855, 1439056, 3678537]
    particleCoord = [(i / 4042929.) * a.N for i in geneCoord]
    particleCoord = [int(i) for i in particleCoord]
    particleCoord = sorted(particleCoord)

    #Below making sure that if two PFRs overlap, we just create a PFR which is twice longer
    gapShift = 7
    gaps = []
    gapStart = particleCoord[0]
    gapEnd = particleCoord[0]
    for i in particleCoord:
        if i > gapEnd:
            gaps.append((gapStart, gapEnd))
            gapStart = i
            gapEnd = i + gapShift
        else:
            gapEnd += gapShift
    gaps.append((gapStart, gapEnd))

    #Adding PFRs at highly expressed genes
    M = domainNum
    particles = []
    for st, end in gaps:
        BD.addGap(st, end)

    #Making bonds
    BD.addBristles(3, plectonemeLength, 0, plectonemeGap)
    BD.sortSegments()
    print BD.segments
    BD.createBonds()
    BD.checkConnectivity()

    BD.save(os.path.join(a.folder, "chains"))
    a.setLayout(mode="chain", chains=BD.getChains())


    a._initHarmonicBondForce()
    for i in BD.bonds:
        a.addBond(i[0], i[1], bondWiggleDistance=0.15, bondType="Harmonic")


    a.addGrosbergStiffness(k=stifness)
    a.addGrosbergRepulsiveForce(trunc=1.)
    #a.addSoftLennardJonesForce(epsilon=0.46, trunc=2.5, cutoff=2.3)
    a.save(os.path.join(a.folder, "start"))
    a.localEnergyMinimization()
    a.save()
    counter = 0

    for step in xrange(500):
        a.doBlock(37000)
        a.save()
    a.printStats()


"""
stifnesses = [5, 6, 8]
gaps = [0, 2, 4]
lengths = [35, 40, 45]
compactions = [2.5, 3, 3.5]
domainNum = [20, 30, 40]
domainLen = [3, 5, 8]
"""
stifnesses = [5, 4, 6]
gaps = [2]
lengths = [30, 35]
compactions = [4, 3.5, 4.5]
domainNum = [30]
domainLen = [5]

stifnesses = [6]
gaps = [0, 2, 4]
lengths = [30, 35]
compactions = [4, 3.5, 4.5]
domainNum = [30]
domainLen = [2, 5, 8]

stifnesses = [6]
gaps = [0, 2, 4]
lengths = [30, 35]
compactions = [4]
domainNum = [20, 30, 10]
domainLen = [10, 5, 15]

stifnesses = [6]
gaps = [2]
lengths = [35]
compactions = [4]
domainNum = [30]
domainLen = [15]

stifnesses = [6, 5, 7]
gaps = [2, 0]
lengths = [35, 32, 38]
compactions = [4, 3.5, 5.5]
domainNum = [30, 25, 35]
domainLen = [15, 12, 18]


stifnesses = [ 6, ]
gaps = [2, ]
lengths = [35]
compactions = [3.5]
domainNum = [15]
domainLen = [6]

stifnesses = [2, 3.5, 5, 8]
gaps = [2, 4, 8, 16]
lengths = [15, 25, 40, 55]
compactions = [2.5, 3, 3.7, 4.5]
domainNum = [15]
domainLen = [6]


stifnesses = [8]
gaps = [8]
lengths = [25]
compactions = [3.7]
inds = range(40)

stifnesses = [7, 8, 9]
gaps = [1, 2, 3]
lengths = [35, 40, 45]
compactions = [3.5, 3.7, 3.9]
inds = range(3)

stifnesses = [6]
gaps = [2]
lengths = [15, 35, 60, 80]
compactions = [3.5]
inds = [1]




combinations = [(i, j, k, l, m) for i in compactions for j in lengths for k in gaps for l in stifnesses for m in inds]
import random
random.shuffle(combinations)
for i in combinations:

    foldername = "newSweep_nokick_lambda{0}_L{1}_gap{2}_stiff{3}_ind{4}".format(*i)
    print foldername
    if os.path.exists(foldername):
        continue
    exampleOpenmm(*i)



