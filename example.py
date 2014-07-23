# (c) 2013 Massachusetts Institute of Technology. All Rights Reserved
# Code written by: Maksim Imakaev (imakaev@mit.edu)

from openmmlib import Simulation
import polymerutils
import os

def exampleOpenmm():
    """
    An example script which generates a compact polymer and expands it to occupy the entire volume

    You need to have a OpenMM-compatible GPU and
    OpenMM installed to run this script.
    Otherwise you can switch to "reference" platform
    a.setup(platform = "reference")

    To install OpenMM please go to their website:
    https://simtk.org/home/openmm

    """

    # a = Simulation(timestep=80, thermostat=0.002)
    # Old version to use with integrator="langevin"

    a = Simulation(thermostat=0.002)
    # timestep not necessary for variableLangevin

    # a.setup(platform="cuda", verbose=True)

    a.setup(platform="cuda", integrator="variableLangevin", errorTol=0.01, verbose=True)
    # If simulation blows up, decrease errorTol twice and try again

    # We use CUDA on 680 GTX, and OpenCL on 580 GTS
    # MirnyLab: Now use "Cuda" on quill, proteome, kulibin, zubr and wiz

    a.saveFolder("trajectory")  # folder where to save trajectory
    # Folder to save trajectory

    # polymer = polymerutils.load("globule")
    # loads compact polymer conformation

    # polymer = polymerutils.grow_rw(8000, 50, method="standard")
    # grows a compact polymer ring of a length 8000 in a 50x50x50 box

    # polymer = polymerutils.create_spiral(r1=4, r2=20, N=8000)
    # Creates a polymer arranged in a cylinder of diameter 20, 8000 monomers long

    polymer = polymerutils.create_random_walk(1, 12000)

    a.load(polymer, center=True)  # loads a polymer, puts a center of mass at zero

    a.save(os.path.join(a.folder, "original"))
    # saves the original file in the same folder

    a.setLayout(mode="chain")
    # This line initializes the fact that we have one chain

    a.addSphericalConfinement(density=0.85, k=1)
    # Specifying density is more intuitive
    # k is the slope of confinement potential, measured in kt/mon
    # set k=5 for harsh confinement
    # and k = 0.2 or something for collapse simulation

    a.addHarmonicPolymerBonds(wiggleDist=0.05)
    # Bonds will fluctuate +- 0.05 on average

    a.addGrosbergRepulsiveForce(trunc=50)
    # this will resolve chain crossings and will not let chain cross anymore

    # a.addGrosbergRepulsiveForce(trunc=5)
    # this will let chains cross sometimes

    a.addStiffness(k=4)
    # K is more or less arbitrary, k=4 corresponds to presistence length of 4,
    # k=1.5 is recommended to make polymer realistically flexible; k=8 is very stiff

    # If your simulation does not start, consider using energy minimization below

    # a.localEnergyMinimization()
    # A very efficient algorithm to reach local energy minimum
    # Use it to minimize energy if you're doing diffusion simulations
    # If you're simulating dynamics of collapse or expansion, please do not use it

    # a.energyMinimization(stepsPerIteration=10)
    # An algorithm to start a simulation
    # Works only with langevin integrator
    # Decreases a timestep for some time and lets simulation settle down

    a.save()
    for _ in xrange(10):
        a.doBlock(2000)
        a.save()
    a.printStats()
    a.show()


exampleOpenmm()
exit()
