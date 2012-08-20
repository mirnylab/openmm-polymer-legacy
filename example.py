from math import sin, cos, sqrt
import numpy
from openmmlib import Simulation
from polymerutils import create_spiral as createSpiral

from polymerutils import load


def exampleOpenmm():
    """
    You need to have a OpenMM-compatible GPU and
    OpenMM installed to run this script.
    Otherwise you can switch to "reference" platform
    a.setup(platform = "reference")
    But this will be extremely slow...

    Installing OpenMM may be not easy too... but you can try
    """

    a = Simulation(timestep=80, thermostat=0.002)
    assert isinstance(a, Simulation)
    a.initStorage("mystorage", mode="w")
    a.setup(platform="OpenCL", verbose=True)
    a.saveFolder("trajectory")  # folder where to save trajectory

    a.load("globule")
    a.setLayout(mode="chain")  # default = chain
    a.addSphericalConfinement(density=0.55)
    a.addHarmonicPolymerBonds()
    a.addGrosbergRepulsiveForce()  # Fastest pure repulsive force
    #a.addGrosbergStiffness()
    a.addStiffness()

    a.addBond(0, 5999, 0.05, 1)

    a.energyMinimization()

    for _ in xrange(10):
        a.doBlock(3000)
        a.save()
    a.printStats()
    a.show()


exampleOpenmm()
exit()

"""Below are scripts that can be used to create starting conformation.
create_spiral will create a spiral.
create_sausage will create a sausage with a certain
length/width ratio of a certain N.
"""
