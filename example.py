from openmmlib import Simulation


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
    # Consider increasing timestep if your system is very relaxed
    assert isinstance(a, Simulation)
    #This line is for Eclipse to know the type of a

    a.setup(platform="OpenCL", verbose=True)
    #Now use "Cuda" on dau, quill, proteome, kulibin, and on wiz CPU #1 (will be changed soon)
    a.saveFolder("trajectory")  # folder where to save trajectory
    #Folder to save trajectory

    a.load("globule")
    a.setLayout(mode="chain")
    #This line initializes the fact that we have one chain

    a.addSphericalConfinement(density=0.85)
    #You can specify radius as well

    a.addHarmonicPolymerBonds(wiggleDist=0.05)
    #Bonds will fluctuate +- 0.05 on average

    #a.addGrosbergRepulsiveForce(trunc=5)  # Fastest pure repulsive force
    #a.addPolynomialRepulsiveForce(trunc=5)  # Fastest pure repulsive force
    #a.addPolynomialAttractionForce(
    #    repulsionEnergy=5.0,
    #    repulsionRadius=1.0,
    #    attractionEnergy=1.0,
    #    attractionRadius=1.5,
    #    )  # Fastest pure repulsive force
    a.addSoftLennardJonesForce(
        epsilon=4.0,
        trunc=5.0,
        cutoff=3.0,
        )  # Fastest pure repulsive force
    #truncation at 50 kT to resolve chain overlaps in the original conformation
    #Optional with this starting conformation, but may be useful in general

    a.addStiffness(k=4)
    #K is more or less arbitrary, k=4 corresponds to presistence length of 4

    a.energyMinimization()
    #Adaptive algorithm to minimize energy

    for _ in xrange(10):
        a.doBlock(3000)
        a.save()
    a.printStats()
    a.show()


exampleOpenmm()
exit()
