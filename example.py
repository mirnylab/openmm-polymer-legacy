#(c) 2013 Massachusetts Institute of Technology. All Rights Reserved
# Code written by: Maksim Imakaev (imakaev@mit.edu)

from openmmlib import Simulation


def exampleOpenmm():
    """
    You need to have a OpenMM-compatible GPU and
    OpenMM installed to run this script.
    Otherwise you can switch to "reference" platform
    a.setup(platform = "reference")
    But this will be extremely slow...

    To install OpenMM please go to their website:
    https://simtk.org/home/openmm

    """

    a = Simulation(timestep=80, thermostat=0.002)
    # Consider increasing timestep if your system is very relaxed

    a.setup(platform="cuda", verbose=True)
    #We use CUDA on 680 GTX, and OpenCL on 580 GTS
    #MirnyLab: Now use "Cuda" on quill, proteome, kulibin, zubr and wiz

    a.saveFolder("trajectory")  # folder where to save trajectory
    #Folder to save trajectory

    a.load("globule")
    a.setLayout(mode="chain")
    #This line initializes the fact that we have one chain

    a.addSphericalConfinement(density=0.85)
    #You can specify radius as well

    a.addHarmonicPolymerBonds(wiggleDist=0.05)
    #Bonds will fluctuate +- 0.05 on average

    a.addGrosbergRepulsiveForce(trunc=50)
    # this will resolve chain crossings and will not let chain cross anymore

    #a.addGrosbergRepulsiveForce(trunc=5)
    # this will let chains cross sometimes

    a.addStiffness(k=4)
    #K is more or less arbitrary, k=4 corresponds to presistence length of 4,

    a.localEnergyMinimization()
    #New fancy algorithm to minimize energy

    for _ in xrange(10):
        a.doBlock(3000)
        a.save()
    a.printStats()
    a.show()


exampleOpenmm()
exit()
