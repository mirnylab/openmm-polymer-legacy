#(c) 2013 Massachusetts Institute of Technology. All Rights Reserved
# Code written by: Maksim Imakaev (imakaev@mit.edu)
#                  Anton Goloborodko (golobor@mit.edu)

"""
Openmm-lib - a wrapper above Openmm to use with polymer simulations
===================================================================

Summary
-------

This is a wrapper above a GPU-assisted molecular dynamics package Openmm.

You can find extensive description of openmm classes here:
https://simtk.org/api_docs/openmm/api10/annotated.html

Input/Output file format
------------------------

Polymer configuration is represented as a Nx3 numpy array of coordinates.
Start/end of chains/rings are not directly specified in the file,
and have to be added through method :py:func:`setLayout <Simulation.setLayout>`

Input file may have the simplistic format, described in txtToJoblib.py (first line with
number of particles, then N lines with three floats corresponding to x,y,z coordinates each).
Input file can also be any of the output files.

Output file format is a dictionary, saved with joblib.dump.
Nx3 data array is stored under the key "data".
The rest of the dictionary consists of metadata, describing details of the simulation.
This metadata is for informative purpose only, and is ignored by the code.

New Input/Output format
-----------------------

A new format of the input data was recently introduced.
Now the data is sotred in a single h5dict file,
where each saved conformation is represented as a single dictionary record.

Keys are strings like "1", "2", etc.

This was done to save some space, and avoid creating millions of files.
This should be actually very fast, and is a recommended storage method.
Show command was modified to incorporate new storage format.

Now if show is runned with two arguments, it assumes that the first is
filename of h5dict, and the second is a number to show.
"-1" stands for the last number.

Implemented forces
------------------

All forces of the system are contained in the self.ForceDict dictionary.
After the force is added, different methods are free to modify parameters
of the force.
Once the system is started, all forces are automatically applied
and cannot be modified.

Two types of bond forces are harmonic bond force
and FENE-type bond as described in Grosberg papers.
Individual bonds can be added using :py:func:`addBond <Simulation.addBond>`,
while polymer bonds can be added using
:py:func:`addHarmonicPolymerBonds <Simulation.addHarmonicPolymerBonds>`, etc.

Repulsive force can be of 3 types: simple repulsife force U = 1/r^12;
Grosberg repulsive force - a faster and better implementation
of repulsive force; and LennardJones Force, that can be attractive and
allows to specify extra attraction/repulsion between any pairs of particles.

Stiffness force can be harmonic, or "special" Grosberg force, kept only
for compatibility with the systems used in Grosberg forces

External forces include spherical confinement, cylindrical confinement,
attraction to "lamina"- surface of a cylinder, gravity, etc.

Information, printed on current step
------------------------------------

A sample line of information printed on each step looks like this:

minim  bl=5 . . . . (i)  pos[1]=[99.2 52.1 52.4]  shift=0.54
4.58 kin, 49.18 pot, 53.76 tot,  Rg=107.312 SPS=144:

Let's go over it step by step.

minim - simulation name (sim -default, minim - energy minimization.
Other name can be provided in self.name).

bl=5 - name of a current block

. . . . . : each dot indicates current subblock. If dots stopped,
maybe the system have crushed

(i) indicate that velocity reinitialization was done at this step.
You will simultaneously see that Ek is more than 2.4

pos[1] is a position of a first monomer

shift is sqrt(mean square displacement) of monomers, i.e.
how much did a monomer shift on average.

4.58 kin, 49.18 pot, 53.76 tot - energies: kinetic, potential, total

Rg=107.312 - current radius of gyration (size of the system)

SPS - steps per second


Functions
---------
Functions depend on each other, and have to be applied in certain groups

1.


:py:func:`load <Simulation.load>`  ---    Mandatory

:py:func:`setup <Simulation.setup>`  ---   Mandatory

:py:func:`setLayout <Simulation.setLayout>` --- Mandatory, after self.load()

:py:func:`saveFolder <Simulation.saveFolder>`  ---  Optional
(default is folder with the code)

2.

self.add___Force()  --- Use any set of forces

Then go all the addBond, and other modifications and tweaks of forces.

3.

Before running actual simulation, it is advised to resolve all possible
conflict by doing :py:func:`energyMinimization <Simulation.energyMinimization>`

4.


:py:func:`doBlock <Simulation.doBlock>`  --- the actual simulation

:py:func:`save <Simulation.save>` --- saves conformation


Frequently-used settings - where to specify them?
-------------------------------------------------

Select GPU ("0" or "1") - :py:func:`setup <Simulation.setup>`

Select chain/ring - :py:func:`setLayout <Simulation.setLayout>`

Select timestep or collision rate - :py:class:`Simulation`


-------------------------------------------------------------------------------

"""
# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php


import numpy
import numpy as np
from polymerutils import getLinkingNumber
import cPickle
import sys
import os
import time
import joblib
import tempfile
import warnings
import polymerutils

os.environ["LD_LIBRARY_PATH"] = "/usr/local/cuda/lib64:/usr/local/openmm/lib"

import simtk.openmm as openmm
import simtk.unit as units

nm = units.meter * 1e-9
fs = units.second * 1e-15
ps = units.second * 1e-12


class Simulation():
    """Base class for openmm simulations

    """
    def __init__(
        self, timestep=80, thermostat=0.001, temperature=300 * units.kelvin,
        verbose=False,
        velocityReinitialize=True,
        # reinitialize velocities at every block if E_kin is more than 2.4
        name="sim",
        length_scale=1.0,
        mass_scale=1.0):  # name to print out
        """

        Parameters
        ----------

        timestep : number
            timestep in femtoseconds. Default value is good.

        thermostat : number
            collision rate in inverse picoseconds. Default value is ok...

        temperature : simtk.units.quantity(units.kelvin), optional
            Temperature of the simulation. Devault value is 300 K.

        verbose : bool, optional
            If True, prints a lot of stuff in the command line.

        velocityReinitialize : bool, optional
            If true, velocities are reinitalized if Ek is more than 2.4 kT.
            Set to False for simulations where inertia is very important.

        name : string, optional
            Name to be printed out as a first line of each block.
            Use it if you run simulations one after another
            and want to see what's going on.

        length_scale : float, optional
            The geometric scaling factor of the system.
            By default, length_scale=1.0 and harmonic bonds and repulsive
            forces have the scale of 1 nm.

        mass_scale : float, optional
            The scaling factor of the mass of the system.
            By default, length_scale=1.0 and harmonic bonds and repulsive forces
            have the scale of 1 nm.


        """

        self.name = name
        self.timestep = timestep * fs
        self.collisionRate = thermostat * (1 / ps)
        self.temperature = temperature
        self.verbose = verbose
        self.velocityReinitialize = velocityReinitialize
        self.loaded = False  # check if the data is loaded
        self.forcesApplied = False
        self.folder = "."
        self.metadata = {}
        self.length_scale = length_scale
        self.mass_scale = mass_scale



    def setup(self, platform="CUDA", PBC=False, PBCbox=None, GPU="default",
              integrator="langevin", verbose=True, errorTol=None):
        """Sets up the important low-level parameters of the platform.
        Mandatory to run.

        Parameters
        ----------

        platform : string, optional
            Platform to use

        PBC : bool, optional
            Use periodic boundary conditions, default:False

        PBCbox : (float,float,float), optional
            Define size of the bounding box for PBC

        GPU : "0" or "1", optional
            Switch to another GPU. Mostly unnecessary.
            Machines with 1 GPU automatically select right GPU.
            Machines with 2 GPUs select GPU that is less used.

        integrator : "langevin", "variableLangevin", "brownian", optional
            Integrator to use (see Openmm class reference)

        verbose : bool, optional
            Shout out loud about every change.

        errorTol : float, optional
            Error tolerance parameter for variableLangevin integrator

        """

        self.step = 0
        if PBC == True:
            self.metadata["PBC"] = True

        self.kB = units.BOLTZMANN_CONSTANT_kB * \
            units.AVOGADRO_CONSTANT_NA  # Boltzmann constant
        self.kT = self.kB * self.temperature  # thermal energy
        self.mass = 100.0 * units.amu * self.mass_scale
        # All masses are the same,
        #changing them would be difficult in this formalism
        self.bondsForException = []
        self.mm = openmm
        self.conlen = 1. * nm * self.length_scale
        self.system = self.mm.System()
        self.PBC = PBC

        if self.PBC == True:  # if periodic boundary conditions
            if PBCbox is None:
                data = self.getData()
                data -= numpy.min(data, axis=0)

                datasize = 1.1 * (2 + (numpy.max(self.getData(), axis=0) - \
                                       numpy.min(self.getData(), axis=0)))
                # size of the system plus some overhead

                self.SolventGridSize = (datasize / 1.1) - 2
                print "density is ", self.N / (datasize[0]
                    * datasize[1] * datasize[2])
            else:
                PBCbox = numpy.array(PBCbox)
                datasize = PBCbox

            self.metadata["PBCbox"] = PBCbox
            self.system.setDefaultPeriodicBoxVectors([datasize[0], 0.,
                0.], [0., datasize[1], 0.], [0., 0., datasize[2]])
            self.BoxSizeReal = datasize

        self.GPU = GPU  # setting default GPU

        if platform.lower() == "opencl":
            platformObject = self.mm.Platform.getPlatformByName('OpenCL')
            if self.GPU.lower() != "default":
                platformObject.setPropertyDefaultValue(
                    'OpenCLDeviceIndex', self.GPU)
            platformObject.setPropertyDefaultValue('OpenCLPrecision', "mixed")

        elif platform.lower() == "reference":
            platformObject = self.mm.Platform.getPlatformByName('Reference')

        elif platform.lower() == "cuda":
            platformObject = self.mm.Platform.getPlatformByName('CUDA')
            if self.GPU.lower() != "default":
                platformObject.setPropertyDefaultValue('CudaDeviceIndex', self.GPU)
            platformObject.setPropertyDefaultValue('CudaPrecision', "mixed")
            platformObject.setPropertyDefaultValue('CudaUseBlockingSync', "true")

        else:
            self.exit("\n!!!!!!!!!!unknown platform!!!!!!!!!!!!!!!\n")
        self.platform = platformObject

        self.forceDict = {}  # Dictionary to store forces
        try:
            for _ in xrange(self.N):
                self.system.addParticle(self.mass)
            print self.N, "particles loaded"
        except:
            pass

        self.integrator_type = integrator
        if integrator.lower() == "langevin":
            self.integrator = self.mm.LangevinIntegrator(self.temperature,
                self.collisionRate, self.timestep)
        elif integrator.lower() == "variablelangevin":
            self.integrator = self.mm.VariableLangevinIntegrator(self.temperature,
                self.collisionRate, errorTol)
        elif integrator.lower() == 'brownian':
            self.integrator = self.mm.BrownianIntegrator(self.temperature,
                self.collisionRate, self.timestep)
        else:
            self.integrator = integrator

    def saveFolder(self, folder):
        """
        sets the folder where to save data.

        Parameters
        ----------
            folder : string
                folder to save the data

        """
        if os.path.exists(folder) == False:
            os.mkdir(folder)
        self.folder = folder

    def exitProgram(self, line):
        """Prints error line and exits program

        Parameters
        ----------

        line : str
            Line to print

        """
        print line
        print "--------------> Bye <---------------"
        exit()

    def setLayout(self, mode="chain", chains=None, Nchains=1):
        """sets layout of chains for chains or rings.
        By default makes one chain. You can change it to one ring (mode=ring).
        You can either have chains/rings of equal length (mode=, Nchains=).
        Or you can have chains or rings of different lengthes (mode=, chains=).
        You can't have a mix of rings and chains

        .. note :: If some monomers are unused in the chains,
        they become freely floating.

        Parameters
        ----------

        mode : "chain" or "ring"
            Does the system consist of rings or chains?

        chains : None or ((0,L1),(L1,L2),(L2,L3)...)
            Specifies exact chain/ring start/end particle numbers,
            if chains are of different lengths. Nchains is ignored.
            E.g. if you have 3 chains of length 5,10,15,
            chains should be [(0,5),(5,15),(15,30)]

        Nchains : int
            Number of chains, if they all are of the same lengths.
            Ignored if chains is specified exactly.


        """

        if mode in ["chain", "ring"]:

            if chains is not None:
                self.chains = chains
            else:
                if not hasattr(self, "N"):
                    raise ValueError("Load the chain first, or provide chain length")
                self.chains = []
                for i in xrange(Nchains):
                    self.chains.append(((self.N * i) /
                        Nchains, (self.N * (i + 1)) / Nchains))
        self.mode = mode
        #print self.N, chains
        layout = {"chains": chains, "mode": mode, "Nchains": Nchains}
        self.metadata["layout"] = repr(layout)

    def getLayout(self):
        "returns configuration of chains"
        return self.chains

    def load(self, filename,  # Input filename, or input data array
             center=False,  # Shift center of mass to zero?
             h5dictKey=None
             ):
        """loads data from file.
        Accepts text files, joblib files or pure data as Nx3 or 3xN array

        If h5dictKey is specified, uses new h5dict-based I/O

        Parameters
        ----------

        filename : joblib file, or text file name, or Nx3 or 3xN numpy array
            Input filename or array with data

        center : bool, optional
            Move center of mass to zero before starting the simulation

        h5dictKey : int or str, optional
            Indicates that you need to load from an h5dict
        """
        if h5dictKey is not None:
            from mirnylib.h5dict import h5dict
            mydict = h5dict(path=filename, mode="r")
            data = mydict[str(h5dictKey)]

        elif type(filename) == str:
            try:
                "checking for a text file here"
                line0 = open(filename).readline()
                try:
                    N = int(line0)
                except ValueError:
                    raise TypeError("Cannot read text file..."\
                                    " reading pickle file")

                lines = open(filename).readlines()[1:]
                data = [[float(i) for i in j.split(
                    )] for j in lines if len(j) > 3]

                if len(data) != N:
                    raise ValueError("N does not correspond"\
                                     " to the number of lines!")

            except TypeError:
                "loading from a joblib file here"
                mydict = dict(joblib.load(filename))
                data = mydict.pop("data")
                self.oldMetadata = data
        else:
            data = filename

        data = numpy.asarray(data, float)

        if len(data) == 3:
            data = numpy.transpose(data)
        if len(data[0]) != 3:
            self.exitProgram("strange data file")
        if numpy.isnan(data).any():
            self.exitProgram("\n!!!!!!!!!!file contains NANS!!!!!!!!!\n")

        if center is True:
            av = numpy.mean(data, 0)
            data -= av

        if center == "zero":
            minvalue = numpy.min(data, 0)
            data -= minvalue

        self.setData(data)
        self.randomizeData()

        if self.verbose == True:
            print "center of mass is", numpy.mean(self.data, 0)
            print "Radius of gyration is,", self.RG()

        try:
            for i in xrange(self.N):
                self.system.addParticle(self.mass)
            if self.verbose == True:
                print "%d particles loaded" % self.N
            self.loaded = True
        except:
            pass

    def save(self, filename=None, mode="auto"):
        """Saves conformation plus some metadata.
        Metadata is not interpreted by this library, and is for your reference

        If data is saved to the .vtf format,
        the same filename should be specified each time.
        .vtf format is then viewable by VMD.

        Parameters:
            mode : str
                "h5dict" : use build-in h5dict storage
                "joblib" : use joblib storage
                "txt" : use text file storage
                "vtf" : append data to the .vtf trajectory file
                "auto" : use h5dict if initialized, joblib otherwise


            filename : str or None
                Filename not needed for h5dict storage
                (use initStorage command) for joblib and txt,
                if filename not provided, it is created automatically.

        """
        mode = mode.lower()

        if mode == "vtf":
            if not hasattr(self, "vtf"):
                with open(filename, 'w') as ff:
                    ff.write("atom 0:%d  radius 1.2 name H \n" % (self.N - 1))
                self.vtf = True
            with open(filename, 'a') as ff:
                ff.write("timestep\n")
                if self.PBC == True:
                    ff.write("pbc %lf %lf %lf" % tuple([
                        2 * i for i in self.BoxSizeReal]))
                data = 2 * self.getScaledData()
                lines = ["%.2lf %.2lf %.2lf \n" % tuple(i) for i in data]
                ff.writelines(lines)
            return

        if mode == "auto":
            if hasattr(self, "storage"):
                mode = "h5dict"
            else:
                mode = "joblib"
        if mode == "h5dict":
            if not hasattr(self, "storage"):
                raise StandardError("Cannot save to h5dict!"\
                                    " Initialize storage first!")
            self.storage[str(self.step)] = self.getScaledData()
            return

        if filename is None:
            if mode == "xyz":
                filename = "block%d.xyz" % self.step
            else:
                filename = "block%d.dat" % self.step
            filename = os.path.join(self.folder, filename)

        if mode == "joblib":
            self.metadata["data"] = self.getScaledData()
            self.metadata["timestep"] = repr(self.timestep / fs)
            self.metadata["Collision rate"] = repr(self.collisionRate / ps)
            joblib.dump(self.metadata, filename=filename, compress=3)

        elif (mode == "txt") or (mode == "xyz"):
            data = self.getScaledData()
            lines = [str(len(data)) + "\n"]
            if mode == "xyz":
                lines.append("\n")
            for j, particle in enumerate(data):
                if mode == "txt":
                    lines.append("".join([str(j) +
                        " " for j in particle]) + "\n")
                else:
                    lines.append("CA " + "".join([str(j)
                        + " " for j in particle]) + "\n")
            with open(filename, 'w') as myfile:
                myfile.writelines(lines)
        else:
            raise ValueError("Unknown mode : %s, use h5dict, joblib or txt" %
                mode)

    def initStorage(self, filename, mode="w-"):
        """
        Initializes an HDF5-based storage for multiple conformation.

        When the storage is initialized,
        save() by default saves to the storage.

        If "r+" mode is chosen, simulation is automaticaly set up
        to continue from existing conformation.
        In that case, last step is automatically determined, and data
        from second-to-last file is loaded.

        .. note :: Continue mode does not copy forces and layouts!

        Parameters
        ----------

        filename : str
            Filename of an h5dict storage file

        mode :
            'w'  - Create file, truncate if exists
            'w-' - Create file, fail if exists         (default)
            'r+' - Continue simulation, file must exist.
        """
        from mirnylib.h5dict import h5dict

        if mode not in ['w', 'w-', 'r+']:
            raise ValueError("Wrong mode to open file."
                             " Only 'w','w-' and 'r+' are supported")
        if (mode == "w-") and os.path.exists(filename):
            raise IOError("Cannot create file... file already exists."\
                          " Use mode ='w' to override")
        self.storage = h5dict(path=filename, mode=mode)
        if mode == "r+":
            myKeys = []
            for i in self.storage.keys():
                try:
                    myKeys.append(int(i))
                except:
                    pass
            maxkey = max(myKeys) if myKeys else 1
            self.step = maxkey - 1
            self.setData(self.storage[str(maxkey - 1)])

    def getData(self):
        "Returns an Nx3 array of positions"
        return numpy.asarray(self.data / nm, dtype="float32")

    def getScaledData(self):
        """Returns data, scaled back to PBC box """
        if self.PBC != True:
            return self.getData()
        alldata = self.getData()
        boxsize = numpy.array(self.BoxSizeReal)
        mults = numpy.floor(alldata / boxsize[None, :])
        toRet = alldata - mults * boxsize[None, :]
        assert toRet.min() >= 0
        return toRet

    def setData(self, data):
        """Sets particle positions

        Parameters
        ----------

        data : Nx3 array-line
            Array of positions with distance ~1 between connected atoms.
        """
        data = numpy.asarray(data, dtype="float")
        self.data = units.Quantity(data, nm)
        self.N = len(self.data)
        if hasattr(self, "context"):
            self.initPositions()

    def randomizeData(self):
        """
        run this if your data is integer-based - adds small offsets
        """
        data = self.getData()
        data = data + numpy.random.randn(*data.shape) * 0.0001
        self.setData(data)

    def RG(self):
        """
        Returns
        -------

        Gyration ratius in units of length (bondlength).
        """
        data = self.getScaledData()
        return numpy.sqrt(numpy.sum(numpy.var(numpy.array(data), 0)))

    def RMAX(self, percentile=None):
        """
        Returns
        -------
        Distance to the furthest from the origin particle.

        """
        data = self.getScaledData()
        dists = numpy.sqrt(numpy.sum((numpy.array(data)) ** 2, 1))
        if percentile == None:
            return numpy.max(dists)
        else:
            return numpy.percentile(dists, percentile)

    def dist(self, i, j):
        """
        Calculates distance between particles i and j
        """
        data = self.getData()
        dif = data[i] - data[j]
        return numpy.sqrt(sum(dif ** 2))

    def useDomains(self, domains=None, filename=None):
        """
        Sets up domains for the simulation.
        Also, pickles domain vector to "domains.dat".

        Parameters
        ----------

        domains : boolean array or None
            N-long array with domain vector
        filename : str or None
            Filename with pickled domain vector

        """

        if domains is not None:
            self.domains = domains

        elif filename is not None:
            self.domains = cPickle.load(open(domains))
        else:
            self.exit("You have to specify at least some domains!")

        if len(self.domains) != self.N:
            self.exitProgram("Wrong domain lengths")

        cPickle.dump(self.domains, open(os.path.join(self.folder,
            "domains.dat"), 'wb'))
        if hasattr(self, "storage"):
            self.storage["domains"] = self.domains

    def _initHarmonicBondForce(self):
        "Internal, inits harmonic forse for polymer and non-polymer bonds"
        if "HarmonicBondForce" not in self.forceDict.keys():
            self.forceDict["HarmonicBondForce"] = self.mm.HarmonicBondForce()
        self.bondType = "Harmonic"

    def _initGrosbergBondForce(self):
        "inits Grosberg FENE bond force"
        if "GrosbergBondForce" not in self.forceDict.keys():
            force = (
                "- 0.5 * GROSk * GROSr0 * GROSr0 * log(1-(r/GROSr0)* (r / GROSr0))"
                " + (4 * GROSe * ((GROSs/r)^12 - (GROSs/r)^6) + GROSe) * step(GROScut - r)")
            bondforceGr = self.mm.CustomBondForce(force)
            bondforceGr.addGlobalParameter("GROSk", 30 *
                self.kT / (self.conlen * self.conlen))
            bondforceGr.addGlobalParameter("GROSr0", self.conlen * 1.5)
            bondforceGr.addGlobalParameter('GROSe', self.kT)
            bondforceGr.addGlobalParameter('GROSs', self.conlen)
            bondforceGr.addGlobalParameter(
                "GROScut", self.conlen * 2. ** (1. / 6.))
            self.forceDict["GrosbergBondForce"] = bondforceGr

    def _initAbsBondForce(self):
        "inits abs(x) FENE bond force"
        if "AbsBondForce" not in self.forceDict.keys():
            force = "(1. / ABSwiggle) * ABSunivK * "\
            "(sqrt((r-ABSr0 * ABSconlen)* "\
            " (r - ABSr0 * ABSconlen) + ABSa * ABSa) - ABSa)"

            bondforceAbs = self.mm.CustomBondForce(force)
            bondforceAbs.addPerBondParameter("ABSwiggle")
            bondforceAbs.addPerBondParameter("ABSr0")
            bondforceAbs.addGlobalParameter("ABSunivK", self.kT / self.conlen)
            bondforceAbs.addGlobalParameter("ABSa", 0.02 * self.conlen)
            bondforceAbs.addGlobalParameter("ABSconlen", self.conlen)
            self.forceDict["AbsBondForce"] = bondforceAbs

    def _initAbsDistanceLimitation(self):
        "inits abs(x) FENE bond force"
        if "AbsLimitation" not in self.forceDict.keys():
            force = (
                "(1. / ABSwiggle) * ABSunivK * step(r - ABSr0 * ABSconlen) "
                "* (sqrt((r-ABSr0 * ABSconlen)"
                "*(r - ABSr0 * ABSconlen) + ABSa * ABSa) - ABSa)")
            bondforceAbsLim = self.mm.CustomBondForce(force)
            bondforceAbsLim.addPerBondParameter("ABSwiggle")
            bondforceAbsLim.addPerBondParameter("ABSr0")
            bondforceAbsLim.addGlobalParameter(
                "ABSunivK", self.kT / self.conlen)
            bondforceAbsLim.addGlobalParameter("ABSa", 0.02 * self.conlen)
            bondforceAbsLim.addGlobalParameter("ABSconlen", self.conlen)
            self.forceDict["AbsLimitation"] = bondforceAbsLim

    def addCenterOfMassRemover(self):
        remover = self.mm.CMMotionRemover(10)
        self.forceDict["CoM_Remover"] = remover

    def addBond(self,
                i, j,  # particles connected by bond
                bondWiggleDistance=0.2,
                # Flexibility of the bond,
                # measured in distance at which energy equals kT
                distance=None,  # Equilibrium length of the bond
                bondType=None,  # Harmonic, Grosberg, ABS
                verbose=None):  # Set this to False if you're in verbose mode
                # and don't want to contaminate output by 10000 messages
        """Adds bond between two particles, allows to specify parameters

        Parameters
        ----------

        i,j : int
            Particle numbers
        bondWiggleDistance : float

            Average displacement from the equilibrium bond distance

        bondType : "Harmonic" or "Grosberg"
            Type of bond. Distance and bondWiggleDistance can be
            specified for harmonic bonds only

        verbose : bool
            Set this to False if you're in verbose mode and don't want to
            print "bond added" message

        """

        if verbose is None:
            verbose = self.verbose
        if (i >= self.N) or (j >= self.N):
            raise ValueError("\nCannot add bond with monomers %d,%d that"\
            "are beyound the polymer length %d" % (i, j, self.N))
        bondSize = float(bondWiggleDistance)
        if distance is None:
            distance = self.conlen / nm
        else:
            distance = self.conlen * distance / nm
        distance = float(distance)

        if bondType is None:
            bondType = self.bondType

        if bondType.lower() == "harmonic":
            self._initHarmonicBondForce()
            kbond = (2 * self.kT / (bondSize * self.conlen)
                ** 2) / (units.kilojoule_per_mole / nm ** 2)
            self.forceDict["HarmonicBondForce"].addBond(
                int(i), int(j), float(distance), float(kbond))

        elif bondType.lower() == "grosberg":
            self._initGrosbergBondForce()
            self.forceDict["GrosbergBondForce"].addBond(int(i), int(j), [])
        elif bondType.lower() == "abs":
            self._initAbsBondForce()
            self.forceDict["AbsBondForce"].addBond(int(i), int(
                j), [float(bondWiggleDistance), float(distance)])
        elif bondType.lower() == "abslim":
            self._initAbsDistanceLimitation()
            self.forceDict["AbsLimitation"].addBond(int(i), int(
                j), [float(bondWiggleDistance), float(distance)])

        else:
            self.exitProgram("Bond type not known")
        if verbose == True:
            print "%s bond added between %d,%d, wiggle %lf dist %lf" % (
                bondType, i, j, float(bondWiggleDistance), float(distance))

    def addHarmonicPolymerBonds(self, wiggleDist=0.05):
        """Adds harmonic bonds connecting polymer chains
        wiggleDist controls the distance at which
        energy of the bond equals kT
        """

        for i in self.chains:
            for j in xrange(i[0], i[1] - 1):
                self.addBond(j, j + 1, wiggleDist, distance=1,
                    bondType="Harmonic", verbose=False)
                self.bondsForException.append((j, j + 1))

            if self.mode == "ring":
                self.addBond(i[0], i[1] - 1, wiggleDist,
                    distance=1, bondType="Harmonic")
                self.bondsForException.append((i[0], i[1] - 1))
                if self.verbose == True:
                    print "ring bond added", i[0], i[1] - 1

        self.metadata["HarmonicPolymerBonds"] = repr({"wiggleDist": wiggleDist})

    def addGrosbergPolymerBonds(self, k=30):
        """Adds FENE bonds according to Halverson-Grosberg paper.
        (Halverson, Jonathan D., et al. "Molecular dynamics simulation study of
         nonconcatenated ring polymers in a melt. I. Statics."
         The Journal of chemical physics 134 (2011): 204904.)

        This method has a repulsive potential build-in,
        so that Grosberg bonds could be used with truncated potentials.
        Is of no use unless you really need to simulate Grosberg-type system.

        Parameters
        ----------
        k : float, optional
            Arbitrary parameter; default value as in Grosberg paper.

         """

        for i in self.chains:
            for j in xrange(i[0], i[1] - 1):
                self.addBond(j, j + 1, bondType="Grosberg")
                self.bondsForException.append((j, j + 1))
            if self.mode == "ring":
                self.addBond(i[0], i[1] - 1, bondType="Grosberg")
                self.bondsForException.append((i[0], i[1] - 1))
                if self.verbose == True:
                    print "ring bond added", i[0], i[1] - 1
        self.metadata["GorsbergPolymerForce"] = repr({"k": k})

    def addStiffness(self, k=1.5):
        """Adds harmonic angle bonds. k specifies energy in kT at one radian
        If k is an array, it has to be of the length N.
        Xth value then specifies stiffness of the angle centered at
        monomer number X.
        Values for ends of the chain will be simply ignored.

        Parameters
        ----------

        k : float or list of length N
            Stiffness of the bond.
            If list, then determines stiffness of the bond at monomer i.
            Potential is k * alpha^2 * 0.5 * kT
        """
        try:
            k[0]
        except:
            k = numpy.zeros(self.N, float) + k
        stiffForce = self.mm.CustomAngleForce(
            "kT*angK * (theta - 3.141592) * (theta - 3.141592) * (0.5)")
        self.forceDict["AngleForce"] = stiffForce
        for i in self.chains:
            for j in xrange(i[0] + 1, i[1] - 1):
                stiffForce.addAngle(j - 1, j, j + 1, [float(k[j])])
            if self.mode == "ring":
                stiffForce.addAngle(i[1] - 2, i[1] - 1, i[0], [k[i[1] - 1]])
                stiffForce.addAngle(i[1] - 1, i[0], i[0] + 1, [k[i[0]]])

        stiffForce.addGlobalParameter("kT", self.kT)
        stiffForce.addPerAngleParameter("angK")
        self.metadata["AngleForce"] = repr({"stiffness": k})

    def addGrosbergStiffness(self, k=1.5):
        """Adds stiffness according to the Grosberg paper.
        (Halverson, Jonathan D., et al. "Molecular dynamics simulation study of
         nonconcatenated ring polymers in a melt. I. Statics."
         The Journal of chemical physics 134 (2011): 204904.)

        Parameters are synchronized with normal stiffness

        If k is an array, it has to be of the length N.
        Xth value then specifies stiffness of the angle centered at
        monomer number X.
        Values for ends of the chain will be simply ignored.

        Parameters
        ----------

        k : float or N-long list of floats
            Synchronized with regular stiffness.
            Default value is very flexible, as in Grosberg paper.
            Default value maximizes entanglement length.

        """
        try:
            k[0]
        except:
            k = numpy.zeros(self.N, float) + k
        stiffForce = self.mm.CustomAngleForce(
            "GRk * kT * (1 - cos(theta - 3.141592))")
        self.forceDict["AngleForce"] = stiffForce

        stiffForce.addGlobalParameter("kT", self.kT)
        stiffForce.addPerAngleParameter("GRk")
        for i in self.chains:
            for j in xrange(i[0] + 1, i[1] - 1):
                stiffForce.addAngle(j - 1, j, j + 1, [k[j]])
            if self.mode == "ring":
                stiffForce.addAngle(i[1] - 2, i[1] - 1, i[0], [k[i[1] - 1]])
                stiffForce.addAngle(i[1] - 1, i[0], i[0] + 1, [k[i[0]]])

        self.metadata["GrosbergAngleForce"] = repr({"stiffness": k})


    def addGrosbergRepulsiveForce(self, trunc=None, radiusMult=1.):
        """This is the fastest non-transparent repulsive force.
        Done according to the paper:
        (Halverson, Jonathan D., et al. "Molecular dynamics simulation study of
         nonconcatenated ring polymers in a melt. I. Statics."
         The Journal of chemical physics 134 (2011): 204904.)
        Parameters
        ----------

        trunc : None or float
             truncation energy in kT, used for chain crossing.
             Value of 1.5 yields frequent passing,
             3 - average passing, 5 - rare passing.

        """
        radius = self.conlen * radiusMult
        self.metadata["GrosbergRepulsiveForce"] = repr({"trunc": trunc})
        nbCutOffDist = radius * 2. ** (1. / 6.)
        if trunc is None:
            repul_energy = "4 * REPe * ((REPsigma/r)^12 - (REPsigma/r)^6) + REPe"
        else:
            repul_energy = (
                "step(REPcut2 - REPU) * REPU"
                " + step(REPU - REPcut2) * REPcut2 * (1 + tanh(REPU/REPcut2 - 1));"
                "REPU = 4 * REPe * ((REPsigma/r2)^12 - (REPsigma/r2)^6) + REPe;"
                "r2 = (r^10. + (REPsigma03)^10.)^0.1")
        self.forceDict["Nonbonded"] = self.mm.CustomNonbondedForce(
            repul_energy)
        repulforceGr = self.forceDict["Nonbonded"]
        repulforceGr.addGlobalParameter('REPe', self.kT)

        repulforceGr.addGlobalParameter('REPsigma', radius)
        repulforceGr.addGlobalParameter('REPsigma03', 0.3 * radius)
        if trunc is not None:
            repulforceGr.addGlobalParameter('REPcut', self.kT * trunc)
            repulforceGr.addGlobalParameter('REPcut2', 0.5 * trunc * self.kT)
        for _ in range(self.N):
            repulforceGr.addParticle(())

        repulforceGr.setCutoffDistance(nbCutOffDist)

    def addPolynomialRepulsiveForce(self, trunc=3.0, radiusMult=1.):
        """This is a simple polynomial repulsive potential. It has the value
        of `trunc` at zero, stays flat until 0.6-0.7 and then drops to zero
        together with its first derivative at r=1.0.

        Parameters
        ----------

        trunc : float
            the energy value around r=0

        """
        radius = self.conlen * radiusMult
        self.metadata["PolynomialRepulsiveForce"] = repr({"trunc": trunc})
        nbCutOffDist = radius
        repul_energy = (
            "rsc12 * (rsc2 - 1.0) * REPe / REPemin + REPe;"
            "rsc12 = rsc4 * rsc4 * rsc4;"
            "rsc4 = rsc2 * rsc2;"
            "rsc2 = rsc * rsc;"
            "rsc = r / REPsigma * REPrmin;")
        self.forceDict["Nonbonded"] = self.mm.CustomNonbondedForce(
            repul_energy)
        repulforceGr = self.forceDict["Nonbonded"]

        repulforceGr.addGlobalParameter('REPe', trunc * self.kT)
        repulforceGr.addGlobalParameter('REPsigma', radius)
        # Coefficients for x^8*(x*x-1)
        #repulforceGr.addGlobalParameter('REPemin', 256.0 / 3125.0)
        #repulforceGr.addGlobalParameter('REPrmin', 2.0 / np.sqrt(5.0))
        # Coefficients for x^12*(x*x-1)
        repulforceGr.addGlobalParameter('REPemin', 46656.0 / 823543.0)
        repulforceGr.addGlobalParameter('REPrmin', np.sqrt(6.0 / 7.0))
        for _ in range(self.N):
            repulforceGr.addParticle(())

        repulforceGr.setCutoffDistance(nbCutOffDist)

    def addSmoothSquareWellForce(self,
        repulsionEnergy=3.0, repulsionRadius=1.,
        attractionEnergy=0.5, attractionRadius=2.0,
        ):
        """
        This is a simple and fast polynomial force that looks like a smoothed
        version of the square-well potential. The energy equals `repulsionEnergy`
        around r=0, stays flat until 0.6-0.7, then drops to zero together
        with its first derivative at r=1.0. After that it drop down to
        `attractionEnergy` and gets back to zero at r=`attractionRadius`.

        The energy function is based on polynomials of 12th power. Both the
        function and its first derivative is continuous everywhere within its
        domain and they both get to zero at the boundary.

        Parameters
        ----------

        repulsionEnergy: float
            the heigth of the repulsive part of the potential.
            E(0) = `repulsionEnergy`
        repulsionRadius: float
            the radius of the repulsive part of the potential.
            E(`repulsionRadius`) = 0,
            E'(`repulsionRadius`) = 0
        attractionEnergy: float
            the depth of the attractive part of the potential.
            E(`repulsionRadius`/2 + `attractionRadius`/2) = `attractionEnergy`
        attractionEnergy: float
            the maximal range of the attractive part of the potential.

        """
        nbCutOffDist = self.conlen * attractionRadius
        self.metadata["PolynomialAttractiveForce"] = repr({"trunc": repulsionEnergy})
        energy = (
            "step(REPsigma - r) * Erep + step(r - REPsigma) * Eattr;"
            ""
            "Erep = rsc12 * (rsc2 - 1.0) * REPe / emin12 + REPe;"
            "rsc12 = rsc4 * rsc4 * rsc4;"
            "rsc4 = rsc2 * rsc2;"
            "rsc2 = rsc * rsc;"
            "rsc = r / REPsigma * rmin12;"
            ""
            "Eattr = - rshft12 * (rshft2 - 1.0) * ATTRe / emin12 - ATTRe;"
            "rshft12 = rshft4 * rshft2 * rshft4;"
            "rshft4 = rshft2 * rshft2;"
            "rshft2 = rshft * rshft;"
            "rshft = (r - REPsigma - ATTRdelta) / ATTRdelta * rmin12"

            )
        self.forceDict["Nonbonded"] = self.mm.CustomNonbondedForce(
            energy)
        repulforceGr = self.forceDict["Nonbonded"]

        repulforceGr.addGlobalParameter('REPe', repulsionEnergy * self.kT)
        repulforceGr.addGlobalParameter('REPsigma', repulsionRadius * self.conlen)

        repulforceGr.addGlobalParameter('ATTRe', attractionEnergy * self.kT)
        repulforceGr.addGlobalParameter('ATTRdelta',
            self.conlen * (attractionRadius - repulsionRadius) / 2.0)
        # Coefficients for the minimum of x^12*(x*x-1)
        repulforceGr.addGlobalParameter('emin12', 46656.0 / 823543.0)
        repulforceGr.addGlobalParameter('rmin12', np.sqrt(6.0 / 7.0))

        for _ in range(self.N):
            repulforceGr.addParticle(())

        repulforceGr.setCutoffDistance(nbCutOffDist)

    def addSmoothSquareWellTailedForce(self,
        repulsionEnergy=3.0, repulsionRadius=1.,
        attractionEnergy=0.5, attractionRadius=2.0,
        tailEnergy=0.1, tailRadius=3.0,
        ):
        """
        This is almost the same potential as in `addSmoothSquareWellTailedForce`.
        The only difference is that the attractive part of the potential
        flattens out to the value of `tailEnergy` at r=`attractionRadius` and
        then goes quadratically to zero at `tailRadius`.
        Please, refer to the documentation for `addSmoothSquareWellForce`
        for the details of the repulsive and attractive parts of the potential.

        Parameters
        ----------

        kwargs:
            same as in `addSmoothSquareWellForce`.
        tailEnergy:
            the depth of the tail part of the potential.
        tailRadius:
            the maximal range of the tail part of the potential.
        """

        self.metadata["PolynomialAttractiveForce"] = repr({"trunc": repulsionEnergy})
        energy = (
            "step(REPsigma - r) * Erep "
            "+ step(r - REPsigma) * step(REPsigma + ATTRdelta - r) * Eattr_inner"
            "+ step(r - REPsigma - ATTRdelta) * step(REPsigma + 2.0 * ATTRdelta - r) * Eattr_outer"
            "+ step(r - REPsigma - ATTRdelta) * Etail;"
            ""
            "Erep = rsc12 * (rsc2 - 1.0) * REPe / emin12 + REPe;"
            "rsc12 = rsc4 * rsc4 * rsc4;"
            "rsc4 = rsc2 * rsc2;"
            "rsc2 = rsc * rsc;"
            "rsc = r / REPsigma * rmin12;"
            ""
            "Eattr_inner = - poly * ATTRe;"
            "Eattr_outer = - poly * (ATTRe - TAILe) - TAILe;"
            "poly = rshft12 * (rshft2 - 1.0) / emin12 + 1.0;"
            "rshft12 = rshft4 * rshft4 * rshft4;"
            "rshft4 = rshft2 * rshft2;"
            "rshft2 = rshft * rshft;"
            "rshft = (r - REPsigma - ATTRdelta) / ATTRdelta * rmin12;"
            ""
            "Etail = - TAILe * rtail * rtail * (rtail - 1.0) * (rtail - 1.0) * 16.0;"
            "rtail = (r - REPsigma - 2 * ATTRdelta) / TAILr / 2.0 + 0.5;"
            )
        self.forceDict["Nonbonded"] = self.mm.CustomNonbondedForce(
            energy)
        repulforceGr = self.forceDict["Nonbonded"]

        repulforceGr.addGlobalParameter('REPe', repulsionEnergy * self.kT)
        repulforceGr.addGlobalParameter('REPsigma', repulsionRadius * self.conlen)

        repulforceGr.addGlobalParameter('ATTRe', attractionEnergy * self.kT)
        repulforceGr.addGlobalParameter('ATTRdelta',
            self.conlen * (attractionRadius - repulsionRadius) / 2.0)

        repulforceGr.addGlobalParameter('TAILe', tailEnergy * self.kT)
        repulforceGr.addGlobalParameter('TAILr', (tailRadius - attractionRadius) * self.kT)

        # Coefficients for the minimum of x^12*(x*x-1)
        repulforceGr.addGlobalParameter('emin12', 46656.0 / 823543.0)
        repulforceGr.addGlobalParameter('rmin12', np.sqrt(6.0 / 7.0))

        for _ in range(self.N):
            repulforceGr.addParticle(())

        repulforceGr.setCutoffDistance(self.conlen * tailRadius)

    def addLennardJonesForce(
        self, cutoff=2.5, domains=False, epsilonRep=0.24, epsilonAttr=0.27,
        blindFraction=(-1), sigmaRep=None, sigmaAttr=None):

        """
        Adds a lennard-jones force, that allows for mutual attraction.
        This is the slowest force out of all repulsive.

        .. note ::
            This is the only force that allows for so-called "exceptions'.
            Exceptions allow you to change parameters of the force
            for a specific pair of particles.
            This can be used to create short-range attraction between
            pairs of particles.
            See manual for Openmm.NonbondedForce.addException.

        Parameters
        ----------

        cutoff : float, optional
            Cutoff value. Default is good.
        domains : bool, optional
            Use domains, defined by
            :py:func:'setDomains <Simulation.setDomains>'
        epsilonRep : float, optional
            Epsilon (attraction strength) for LJ-force for all particles
            (except for domain) in kT
        epsilonAttr : float, optional
            Epsilon for attractive domain (if domains are used) in kT
        blindFraction : float, 0<x<1
            Fraction of particles that are "transparent" -
            used here instead of truncation
        sigmaRep, sigmaAttr: float, optional
            Radius of particles in the LJ force. For advanced fine-tuning.

         """
        self.metadata["LennardJonesForce"] = repr({"cutoff": cutoff,
                  "domains": domains, "epsilonRep": epsilonRep,
                  "epsilonAttr": epsilonAttr, "blindFraction": blindFraction})

        if blindFraction > 0.99:
            self.exitProgram("why do you need this force without particles???"\
                             " set blindFraction between 0 and 1")
        if (sigmaRep is None) and (sigmaAttr is None):
            sigmaAttr = sigmaRep = self.conlen
        else:
            sigmaAttr = sigmaAttr * self.conlen
            sigmaRep = sigmaRep * self.conlen

        epsilonRep = epsilonRep * self.kT
        epsilonAttr = epsilonAttr * self.kT

        nbCutOffDist = self.conlen * cutoff
        self.epsilonRep = epsilonRep
        repulforce = self.mm.NonbondedForce()

        self.forceDict["Nonbonded"] = repulforce
        for i in xrange(self.N):
            particleParameters = [0., 0., 0.]

            if numpy.random.random() > blindFraction:
                particleParameters[1] = (sigmaRep)
                particleParameters[2] = (epsilonRep)

                if domains == True:
                    if self.domains[i] != 0:
                        particleParameters[1] = (sigmaAttr)
                        particleParameters[2] = (epsilonAttr)

            repulforce.addParticle(*particleParameters)

        repulforce.setCutoffDistance(nbCutOffDist)

    def addSoftLennardJonesForce(self, epsilon=0.42, trunc=2, cutoff=2.5):
        """A softened version of lennard-Jones force.
        Now we're moving to polynomial forces, so go there instead.
        """

        nbCutOffDist = self.conlen * cutoff

        repul_energy = (
            'step(REPcut2 - REPU) * REPU +'
            ' step(REPU - REPcut2) * REPcut2 * (1 + tanh(REPU/REPcut2 - 1));'
            'REPU = 4 * REPe * ((REPsigma/r2)^12 - (REPsigma/r2)^6);'
            'r2 = (r^10. + (REPsigma03)^10.)^0.1')
        self.forceDict["Nonbonded"] = self.mm.CustomNonbondedForce(
            repul_energy)
        repulforceGr = self.forceDict["Nonbonded"]
        repulforceGr.addGlobalParameter('REPe', self.kT * epsilon)

        repulforceGr.addGlobalParameter('REPsigma', self.conlen)
        repulforceGr.addGlobalParameter('REPsigma03', 0.3 * self.conlen)
        repulforceGr.addGlobalParameter('REPcut', self.kT * trunc)
        repulforceGr.addGlobalParameter('REPcut2', 0.5 * trunc * self.kT)

        for _ in range(self.N):
            repulforceGr.addParticle(())

        repulforceGr.setCutoffDistance(nbCutOffDist)

    def addMutualException(self, particles):
        """used to exclude a bunch of particles
        from calculation of nonbonded force

        Parameters
        ----------
        particles : list
            List of particles for whom to exclude nonbonded force.
        """
        for i in particles:  # xrange(len(particles)):
            for j in particles:  # xrange(len(particles)):
                if j > i:
                    self.bondsForException.append((i, j))

    def addInteraction(self, i, j, epsilon, sigma=None, length=3):
        """Adds attractive short-range interaction of strength epsilon
        between particles i,j and a few neighboring particles
        requires :py:func:'LennardJones Force<Simulation.addLennardJonesForce>'

        Parameters
        ----------
        i,j : int
            Interacting particles
        epsilon : float
            LJ strength
        sigma : float, optional
            LJ length. If you increase it past 1.5, note the cutoff!
        length : int, optional, default = 3
            Number of particles around i,j that also attract each other

        """
        if "interactions" not in self.metadata:
            self.metadata["interactions"] = []
        self.metadata["interactions"].append((i, j))

        if type(self.forceDict["Nonbonded"]) != self.mm.NonbondedForce:
            self.exit("Cannot add interactions"\
                      " without Lennard-Jones nonbonded force")

        if sigma is None:
            sigma = 1.1 * self.conlen
        epsilon = epsilon * units.kilocalorie_per_mole
        if (min(i, j) < length) or (max(i, j) > self.N - length):
            print "!!!!!!!!!bond with %d and %d is out of range!!!!!" % (i, j)
            return
        repulforce = self.forceDict["Nonbonded"]
        for t1 in xrange(i - length / 2, i + (length - length / 2)):
            for t2 in xrange(j - length / 2, j + (length - length / 2)):
                repulforce.addException(t1, t2, 0, sigma, epsilon, True)
                if self.verbose == True:
                    print "Exception added between"\
                    " particles %d and %d" % (t1, t2)

        for tt in xrange(i - length, i + length):
            repulforce.setParticleParameters(
                tt, 0, self.conlen, self.epsilonRep)
        for tt in xrange(j - length, j + length):
            repulforce.setParticleParameters(
                tt, 0, self.conlen, self.epsilonRep)

    def addCylindricalConfinement(self, r, bottom=None, k=0.1, top=9999):
        "As it says."

        if bottom == True:
            warnings.warn(DeprecationWarning(
                "Use bottom=0 instead of bottom = True! "))
            bottom = 0

        self.metadata["CylindricalConfinement"] = repr({"r": r,
            "bottom": bottom, "k": k, "top" : top})

        if bottom is not None:
            extforce2 = self.mm.CustomExternalForce(
                "step(r-CYLaa) * CYLkb * (sqrt((r-CYLaa)*(r-CYLaa) + CYLt*CYLt) - CYLt)"
                "+ step(-z + CYLbot) * CYLkb * (sqrt((z - CYLbot)^2 + CYLt^2) - CYLt) "
                "+ step(z - CYLtop) * CYLkb * (sqrt((z - CYLtop)^2 + CYLt^2) - CYLt);"
                "r = sqrt(x^2 + y^2 + CYLtt^2)")
        else:
            extforce2 = self.mm.CustomExternalForce(
                "step(r-CYLaa) * CYLkb * (sqrt((r-CYLaa)*(r-CYLaa) + CYLt*CYLt) - CYLt);"
                "r = sqrt(x^2 + y^2 + CYLtt^2)")

        self.forceDict["CylindricalConfinement"] = extforce2
        for i in xrange(self.N):
            extforce2.addParticle(i, [])
        extforce2.addGlobalParameter("CYLkb", k * self.kT / nm)
        extforce2.addGlobalParameter("CYLtop", top * self.conlen)
        if bottom is not None:
            extforce2.addGlobalParameter("CYLbot", bottom * self.conlen)
        extforce2.addGlobalParameter("CYLkt", self.kT)
        extforce2.addGlobalParameter("CYLweired", nm)
        extforce2.addGlobalParameter("CYLaa", (r - 1. / k) * nm)
        extforce2.addGlobalParameter("CYLt", (1. / (10 * k)) * nm)
        extforce2.addGlobalParameter("CYLtt", 0.01 * nm)

    def addSphericalConfinement(self,
                r="density",  # radius... by default uses certain density
                k=5.,  # How steep the walls are
                density=.3):  # target density, measured in particles
                                #per cubic nanometer (bond size is 1 nm)
        """Constrain particles to be within a sphere.
        With no parameters creates sphere with density .3

        Parameters
        ----------
        r : float or "density", optional
            Radius of confining sphere. If "density" requires density,
            or assumes density = .3
        k : float, optional
            Steepness of the confining potential, in kT/nm
        density : float, optional, <1
            Density for autodetection of confining radius.
            Density is calculated in particles per nm^3,
            i.e. at density 1 each sphere has a 1x1x1 cube.
        """
        self.metadata["SphericalConfinement"] = repr({"r": r, "k": k,
            "density": density})

        spherForce = self.mm.CustomExternalForce(
            "step(r-SPHaa) * SPHkb * (sqrt((r-SPHaa)*(r-SPHaa) + SPHt*SPHt) - SPHt) "
            ";r = sqrt(x^2 + y^2 + z^2 + SPHtt^2)")
        self.forceDict["SphericalConfinement"] = spherForce

        for i in xrange(self.N):
            spherForce.addParticle(i, [])
        if r == "density":
            r = (3 * self.N / (4 * 3.141592 * density)) ** (1 / 3.)

        self.sphericalConfinementRadius = r
        if self.verbose == True:
            print "Spherical confinement with radius = %lf" % r
        #assigning parameters of the force
        spherForce.addGlobalParameter("SPHkb", k * self.kT / nm)
        spherForce.addGlobalParameter("SPHaa", (r - 1. / k) * nm)
        spherForce.addGlobalParameter("SPHt", (1. / k) * nm / 10.)
        spherForce.addGlobalParameter("SPHtt", 0.01 * nm)
        return r

    def excludeSphere(self, r=5, position=(0, 0, 0)):
        """Excludes particles from a sphere of radius r at certain position.
        """

        spherForce = self.mm.CustomExternalForce(
            "step(EXaa-r) * EXkb * (sqrt((r-EXaa)*(r-EXaa) + EXt*EXt) - EXt) ;"
            "r = sqrt((x-EXx)^2 + (y-EXy)^2 + (z-EXz)^2 + EXtt^2)")
        self.forceDict["ExcludeSphere"] = spherForce

        for i in xrange(self.N):
            spherForce.addParticle(i, [])

        self.sphericalConfinementRadius = r
        if self.verbose == True:
            print "Spherical confinement with radius = %lf" % r
        #assigning parameters of the force
        spherForce.addGlobalParameter("EXkb", 2 * self.kT / nm)
        spherForce.addGlobalParameter("EXaa", (r - 1. / 3) * nm)
        spherForce.addGlobalParameter("EXt", (1. / 3) * nm / 10.)
        spherForce.addGlobalParameter("EXtt", 0.01 * nm)
        spherForce.addGlobalParameter("EXx", position[0] * self.conlen)
        spherForce.addGlobalParameter("EXy", position[1] * self.conlen)
        spherForce.addGlobalParameter("EXz", position[2] * self.conlen)

    def addLaminaAttraction(self, width=1, depth=1, r=None):
        """Attracts one domain to the lamina. Infers radius
        from spherical confinement, that has to be initialized already.

        Parameters
        ----------

        width : float, optional
            Width of attractive layer next to the lamina, nm.
        depth : float, optional
            Depth of attractive potential in kT
        r : float, optional
            Radius of an attractive cage. If not specified, inferred
            from previously defined spherical potential.
        """

        self.metadata["laminaAttraction"] = repr({"width": width,
            "depth": depth, "r": r})
        laminaForce = self.mm.CustomExternalForce(
            "step(LAMr-LAMaa + LAMwidth) * step(LAMaa + LAMwidth - LAMr) "
            "* LAMdepth * (LAMr-LAMaa + LAMwidth) * (LAMaa + LAMwidth - LAMr) "
            "/ (LAMwidth * LAMwidth);"
            "LAMr = sqrt(x^2 + y^2 + z^2 + LAMtt^2)")
        self.forceDict["Lamina attraction"] = laminaForce

        #adding all the particles on which force acts
        for i in xrange(self.N):
            if self.domains[i] > 0.5:
                laminaForce.addParticle(i, [])
        if r is None:
            try:
                r = self.sphericalConfinementRadius
            except:
                raise ValueError("No spherical confinement radius defined"\
                                 " yet. Apply spherical confinement first!")
        if self.verbose == True:
            print "Lamina attraction added with r = %d" % r

        laminaForce.addGlobalParameter("LAMaa", r * nm)
        laminaForce.addGlobalParameter("LAMwidth", width * nm)
        laminaForce.addGlobalParameter("LAMdepth", depth * self.kT)
        laminaForce.addGlobalParameter("LAMtt", 0.01 * nm)

    def tetherParticles(self, particles, k=30):
        """tethers particles in the 'particles' array.
        Increase k to tether them stronger, but watch the system!

        Parameters
        ----------

        particles : list of ints
            List of particles to be tethered (fixed in space)
        k : int, optional
            Steepness of the tethering potential.
            Values >30 will require decreasing potential,
            but will make tethering rock solid.
        """
        self.metadata["TetheredParticles"] = repr({"particles": particles, "k": k})
        if "Tethering Force" not in self.forceDict:
            tetherForce = self.mm.CustomExternalForce(
              " TETHkb * ((x - TETHx0)^2 + (y - TETHy0)^2 + (z - TETHz0)^2)")
            self.forceDict["Tethering Force"] = tetherForce
        else:
            tetherForce = self.forceDict["Tethering Force"]

        #assigning parameters of the force
        tetherForce.addGlobalParameter("TETHkb", k * self.kT / nm)
        tetherForce.addPerParticleParameter("TETHx0")
        tetherForce.addPerParticleParameter("TETHy0")
        tetherForce.addPerParticleParameter("TETHz0")
        for i in particles:  # adding all the particles on which force acts
            i = int(i)
            coordinates = self.data[i]
            tetherForce.addParticle(i, list(coordinates))
            if self.verbose == True:
                print "particle %d tethered! " % i

    def addGravity(self, k=0.1, cutoff=None):
        """adds force pulling downwards in z direction
        When using cutoff, acts only when z>cutoff"""
        self.metadata["gravity"] = repr({"k": k, "cutoff": cutoff})
        if cutoff is None:
            gravity = self.mm.CustomExternalForce("kG * z")
        else:
            gravity = self.mm.CustomExternalForce(
                "kG * (z - cutoffG) * step(z - cutoffG)")
            gravity.addGlobalParameter("cutoffG", cutoff * nm)
        gravity.addGlobalParameter("kG", k * self.kT / (nm))

        for i in xrange(self.N):
            gravity.addParticle(i, [])
        self.forceDict["Gravity"] = gravity

    def addPullForce(self, particles, forces):
        """adds force pulling on each particle
        """

        pullForce = self.mm.CustomExternalForce(
            "PULLx * x + PULLy * y + PULLz * z")
        pullForce.addPerParticleParameter("PULLx")
        pullForce.addPerParticleParameter("PULLy")
        pullForce.addPerParticleParameter("PULLz")
        for num, force in map(None, particles, forces):
            force = [float(i) * (self.kT / self.conlen) for i in force]
            pullForce.addParticle(num, force)
        self.forceDict["PullForce"] = pullForce

    def _applyForces(self):

        if self.forcesApplied == True:
            return
        """Applies all the forces in the forcedict.
        Forces should not be modified after that, unless you do it carefully
        (see openmm reference)."""

        exc = self.bondsForException
        print "Number of exceptions:", len(exc)

        if len(exc) > 0:
            exc = numpy.array(exc)
            exc = numpy.sort(exc, axis=1)
            exc = [tuple(i) for i in exc]
            exc = list(set(exc))  # only unique pairs are left

        for i in self.forceDict.keys():  # Adding exceptions
            force = self.forceDict[i]
            if hasattr(force, "addException"):
                print 'Add exceptions for {0} force'.format(i)
                for pair in exc:
                    force.addException(int(pair[0]),
                        int(pair[1]), 0, 0, 0, True)
            elif hasattr(force, "addExclusion"):
                print 'Add exclusions for {0} force'.format(i)
                for pair in exc:
                    #force.addExclusion(*pair)
                    force.addExclusion(int(pair[0]), int(pair[1]))

            if hasattr(force, "CutoffNonPeriodic") and hasattr(
                                                    force, "CutoffPeriodic"):
                if self.PBC:
                    force.setNonbondedMethod(force.CutoffPeriodic)
                    print "Using periodic boundary conditions!!!!"
                else:
                    force.setNonbondedMethod(force.CutoffNonPeriodic)
            print "adding force ", i, self.system.addForce(self.forceDict[i])

        self.context = self.mm.Context(
            self.system, self.integrator, self.platform)
        self.initPositions()
        self.initVelocities()
        self.forcesApplied = True
        if hasattr(self, "storage") and hasattr(self, "metadata"):
            self.storage["metadata"] = self.metadata

    def initVelocities(self, mult=1):
        """Initializes particles velocities

        Parameters
        ----------
        mult : float, optional
            Multiply velosities by this. Is good for a cold/hot start.
        """
        try:
            self.context
        except:
            raise ValueError("No context, cannot set velocs."\
                             "Initialize context before that")

        sigma = units.sqrt(self.kT / self.system.getParticleMass(
            0))  # calculating mean velocity
        velocs = units.Quantity(mult * numpy.random.normal(
            size=(self.N, 3)), units.meter) * (sigma / units.meter)
        #Guide to simtk.unit: 1. Always use units.quantity.
        #2. Avoid dimensionless shit.
        #3. If you have to, create fake units, as done here with meters
        self.context.setVelocities(velocs)

    def initPositions(self):
        """Sends particle coordinates to OpenMM system.
        If system has exploded, this is
         used in the code to reset coordinates. """

        print "Positions... ",
        try:
            self.context
        except:
            raise ValueError("No context, cannot set velocs."\
                             " Initialize context before that")

        self.context.setPositions(self.data)
        print " loaded!",
        state = self.context.getState(getPositions=True, getEnergy=True)
            #get state of a system: positions, energies
        eP = state.getPotentialEnergy() / self.N / self.kT
        print "potential energy is %lf" % eP

    def reinitialize(self, mult=1):
        """Reinitializes the OpenMM context object.
        This should be called if low-level parameters,
        such as forces, have changed.

        Parameters
        ----------
        mult : float, optional
            mult to be passed to
             :py:func:'initVelocities <Simulation.initVelocities>'
        """
        self.context.reinitialize()
        self.initPositions()
        self.initVelocities(mult)

    def localEnergyMinimization(self, tolerance=1, maxIterations=0):
        "A wrapper to the build-in OpenMM Local Energy Minimization"
        print "Performing local energy minimization"
        self._applyForces()
        oldName = self.name
        self.name = "minim"

        self.state = self.context.getState(getPositions=False,
                                           getEnergy=True)
        eK = (self.state.getKineticEnergy() / self.N / self.kT)
        eP = self.state.getPotentialEnergy() / self.N / self.kT
        locTime = self.state.getTime()
        print "before minimization eK={0}, eP={1}, time={2}".format(eK, eP, locTime)

        self.mm.LocalEnergyMinimizer.minimize(
            self.context, tolerance, maxIterations)

        self.state = self.context.getState(getPositions=False,
                                           getEnergy=True)
        eK = (self.state.getKineticEnergy() / self.N / self.kT)
        eP = self.state.getPotentialEnergy() / self.N / self.kT
        locTime = self.state.getTime()
        print "after minimization eK={0}, eP={1}, time={2}".format(eK, eP, locTime)

        self.name = oldName

    def energyMinimization(self, stepsPerIteration=100,
                           maxIterations=1000,
                           failNotConverged=True):
        """Runs system at smaller timestep and higher collision
        rate to resolve possible conflicts.

        Now we're moving towards local energy minimization,
        this is here for backwards compatibility.
        """
        warnings.warn(DeprecationWarning("Maybe use local energy minimization instead - it is better!"))
        print "Performing energy minimization"
        self._applyForces()
        oldName = self.name
        self.name = "minim"
        if (maxIterations is True) or (maxIterations is False):
            raise ValueError(
                "Please stop using the old notation and read the new energy minimization code")
        if (failNotConverged is not True) and (failNotConverged is not False):
            raise ValueError(
                "Please stop using the old notation and read the new energy minimization code")

        def_step = self.integrator.getStepSize()
        def_fric = self.integrator.getFriction()

        def minimizeDrop():
            drop = 10.
            for dummy in xrange(maxIterations):
                if drop < 1:
                    drop = 1.
                if drop > 10000:
                    raise RuntimeError("Timestep too low. Perhaps, "\
                                       "something is wrong!")

                self.integrator.setStepSize(def_step / float(drop))
                self.integrator.setFriction(def_fric * drop)
                #self.reinitialize()
                numAttempts = 5
                for attempt in xrange(numAttempts):
                    a = self.doBlock(stepsPerIteration, increment=False,
                        reinitialize=False)
                    #self.initVelocities()
                    if a == False:
                        drop *= 2
                        print "Drop increased to {0}".format(drop)
                        self.initVelocities()
                        break
                    if attempt == numAttempts - 1:
                        if drop == 1.:
                            return 0
                        drop /= 2
                        print "Drop decreased to {0}".format(drop)
                        self.initVelocities()
            return -1

        if failNotConverged and (minimizeDrop() == -1):
            raise RuntimeError(
                "Reached maximum number of iterations and still not converged\n"\
                "increase maxIterations or set failNotConverged=False")
        self.name = oldName
        self.integrator.setFriction(def_fric)
        self.integrator.setStepSize(def_step)
        #self.reinitialize()
        print "Finished energy minimization"

    def doBlock(self, steps=None, increment=True, num=None, reinitialize=True):
        """performs one block of simulations, doing steps timesteps,
        or steps_per_block if not specified.

        Parameters
        ----------

        steps : int or None
            Number of timesteps to perform. If not specified, tries to
            infer it from self.steps_per_block
        increment : bool, optional
            If true, will not increment self.steps counter
        num : int or None, optional
            If specified, will split the block in num subblocks.
            Default value is 10.
        """

        if self.forcesApplied == False:
            self._applyForces()
            self.forcesApplied = True
        if increment == True:
            self.step += 1
        if steps is None:
            steps = self.steps_per_block
        if (increment == True) and ((self.step % 50) == 25):
            self.printStats()

        for attempt in xrange(6):
            print "%s  bl=%d" % (self.name, self.step),
            if num is None:
                num = steps / 5 + 1
            a = time.time()
            for _ in xrange(steps / num):
                self.integrator.step(num)  # integrate!
                sys.stdout.flush()
            if (steps % num) > 0:
                self.integrator.step(steps % num)

            # get state of a system: positions, energies
            self.state = self.context.getState(getPositions=True,
                                               getEnergy=True)

            b = time.time()
            coords = self.state.getPositions(asNumpy=True)
            newcoords = coords / nm
            # calculate energies in KT/particle
            eK = (self.state.getKineticEnergy() / self.N / self.kT)
            eP = self.state.getPotentialEnergy() / self.N / self.kT

            if self.velocityReinitialize:
                if eK > 2.4:
                    print "(i)",
                    self.initVelocities()
            print "pos[1]=[%.1lf %.1lf %.1lf]" % tuple(newcoords[0]),

            if ((numpy.isnan(newcoords).any()) or (eK > 200) or
                (numpy.isnan(eK)) or (numpy.isnan(eP))):

                self.context.setPositions(self.data)
                self.initVelocities()
                if reinitialize == False:
                    print "eK={0}, eP={1}".format(eK, eP)
                    return False
                print "eK={0}, eP={1}, trying one more time at step {2} ".format(eK, eP, self.step)
            else:
                dif = numpy.sqrt(numpy.mean(numpy.sum((newcoords -
                    self.getData()) ** 2, axis=1)))
                print "dr=%.2lf" % (dif,),
                self.data = coords
                print "kin=%.2lf pot=%.2lf" % (eK,
                    eP), "Rg=%.3lf" % self.RG(),
                print "SPS=%.0lf" % (steps / (float(b - a)))
                break

            if attempt in [3, 4]:
                self.localEnergyMinimization()
            if attempt == 5:
                self.exitProgram("exceeded number of attempts")
        return True

    def printStats(self):
        """Prints detailed statistics of a system.
        Will be run every 50 steps
        """
        state = self.context.getState(getPositions=True,
            getVelocities=True, getEnergy=True)

        eP = state.getPotentialEnergy()
        pos = numpy.array(state.getPositions() / nm)
        bonds = numpy.sqrt(numpy.sum(numpy.diff(pos, axis=0) ** 2, axis=1))
        sbonds = numpy.sort(bonds)
        vel = state.getVelocities()
        mass = self.system.getParticleMass(0)
        vkT = numpy.array(vel / units.sqrt(self.kT / mass), dtype=float)
        self.velocs = vkT
        EkPerParticle = 0.5 * numpy.sum(vkT ** 2, axis=1)

        cm = numpy.mean(pos, axis=0)
        centredPos = pos - cm[None, :]
        dists = numpy.sqrt(numpy.sum(centredPos ** 2, axis=1))
        per95 = numpy.percentile(dists, 95)
        den = (0.95 * self.N) / ((4. * numpy.pi * per95 ** 3) / 3)
        per5 = numpy.percentile(dists, 5)
        den5 = (0.05 * self.N) / ((4. * numpy.pi * per5 ** 3) / 3)
        x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
        minmedmax = lambda x: (x.min(), numpy.median(x), x.mean(), x.max())

        print
        print "Statistics for the simulation %s, number of particles: %d, "\
        " number of chains: %d,  mode:  %s" % (
            self.name, self.N, len(self.chains), self.mode)
        print
        print "Statistics for particle position"
        print "     mean position is: ", numpy.mean(
            pos, axis=0), "  Rg = ", self.RG()
        print "     median bond size is ", numpy.median(bonds)
        print "     three shortest/longest (<10)/ bonds are ", sbonds[
            :3], "  ", sbonds[sbonds < 10][-3:]
        if (sbonds > 10).sum() > 0:
            print "longest 10 bonds are", sbonds[-10:]

        print "     95 percentile of distance to center is:   ", per95
        print "     density of closest 95% monomers is:   ", den
        print "     density of the core monomers is:   ", den5
        print "     min/median/mean/max coordinates are: "
        print "     x: %.2lf, %.2lf, %.2lf, %.2lf" % minmedmax(x)
        print "     y: %.2lf, %.2lf, %.2lf, %.2lf" % minmedmax(y)
        print "     z: %.2lf, %.2lf, %.2lf, %.2lf" % minmedmax(z)
        print
        print "Statistics for velocities:"
        print "     mean kinetic energy is: ", numpy.mean(
            EkPerParticle), "should be:", 1.5
        print "     fastest particles are (in kT): ", numpy.sort(
            EkPerParticle)[-5:]

        print
        print "Statistics for the system:"
        print "     Forces are: ", self.forceDict.keys()
        print "     Number of exceptions:  ", len(self.bondsForException)
        print
        print "Potential Energy Ep = ", eP / self.N / self.kT

    def show(self, shifts=[0., 0.2, 0.4, 0.6, 0.8]):
        """shows system in rasmol by drawing spheres
        draws 4 spheres in between any two points (5 * N spheres total)
        """

        #if you want to change positions of the spheres along each segment,
        #change these numbers: e.g. [0,.1, .2 ...  .9] will draw 10 spheres,
        # and this will look better

        data = self.getData()
        if len(data[0]) != 3:
            data = numpy.transpose(data)
        if len(data[0]) != 3:
            print "wrong data!"
            return
        #determining the 95 percentile distance between particles,
        meandist = numpy.percentile(numpy.sqrt(
            numpy.sum(numpy.diff(data, axis=0) ** 2, axis=1)), 95)
        #rescaling the data, so that bonds are of the order of 1.
        #This is because rasmol spheres are of the fixed diameter.
        data /= meandist

        if self.N > 1000:  # system is sufficiently large
            count = 0
            for _ in xrange(100):
                a, b = numpy.random.randint(0, self.N, 2)
                dist = numpy.sqrt(numpy.sum((data[a] - data[b]) ** 2))
                if dist < 1.3:
                    count += 1
            if count > 100:
                raise RuntimeError(
                    "Too many particles are close together. "\
                    "This will cause rasmol to choke")

        rascript = tempfile.NamedTemporaryFile()
        # writing the rasmol script. Spacefill controls radius of the sphere.
        rascript.write("""wireframe off
        color temperature
        spacefill 100
        background white
        """)
        rascript.flush()

        # creating the array, linearly chanhing from -225 to 225
        # to serve as an array of colors
        colors = numpy.array([int((j * 450.) / (len(data))) -
            225 for j in xrange(len(data))])

        #creating spheres along the trajectory
        newData = numpy.zeros(
            (len(data) * len(shifts) - (len(shifts) - 1), 4))
        for i in xrange(len(shifts)):
            newData[i:-1:len(shifts), :3] = data[:-1] * shifts[
                i] + data[1:] * (1 - shifts[i])
            newData[i:-1:len(shifts), 3] = colors[:-1]
        newData[-1, :3] = data[-1]
        newData[-1, 3] = colors[-1]

        towrite = tempfile.NamedTemporaryFile()
        towrite.write("%d\n\n" % (len(newData)))
        # number of atoms and a blank line after is a requirement of rasmol

        for i in newData:
            towrite.write("CA\t%lf\t%lf\t%lf\t%d\n" % tuple(i))
        towrite.flush()
        "TODO: rewrite using subprocess.popen"

        if os.name == "posix":  # if linux
            os.system("rasmol -xyz %s -script %s" % (
                towrite.name, rascript.name))
        else:  # if windows
            os.system("C:/RasWin/raswin.exe -xyz %s -script %s" % (
                                        towrite.name, rascript.name))

        rascript.close()
        towrite.close()


class SimulationWithCrosslinks(Simulation):
    """
    A subclass used to do simulations for the Metaphase project
    """

    def addConsecutiveRandomBonds(self, loopSize, bondWiggle, bondLength=0.,
                                  smeerLoopSize=0.2, distanceBetweenBonds=2,
                                  verbose=False):
        shift = int(loopSize * smeerLoopSize)
        if shift == 0:
            shift = 1
        begin = numpy.random.randint(distanceBetweenBonds)
        while True:
            b1 = begin
            b2 = begin + loopSize + numpy.random.randint(shift)
            if b2 > self.N - 3:
                if (self.N - b1) > (5 * distanceBetweenBonds + 5):
                    b2 = self.N - 1 - numpy.random.randint(distanceBetweenBonds)
                else:
                    break

            self.addBond(b1, b2, bondWiggle, bondLength,
                         verbose=verbose)
            begin = b2 + numpy.random.randint(distanceBetweenBonds)
            if self.verbose == True:
                print "bond added between %d and %d" % (b1, b2)

    def addDoubleRandomLengthBonds(self, bondlength, bondRange, distance):
        begin = 4
        started = True
        past = 0
        while True:
            past
            b1 = begin
            b2 = begin + numpy.random.randint(
                0.5 * bondlength, 1.7 * bondlength)
            if b2 > self.N - 4:
                break
            self.addBond(b1, b2, bondRange, distance)
            if self.verbose == True:
                print "bond added between %d and %d" % (b1, b2)
            if started == False:
                self.addBond(past, b2, bondRange, distance)
                if self.verbose == True:
                    print "bond added between %d and %d" % (past, b2)
                past = b1
            started = False
            begin = b2

    def addAttractionToTheCore(self, k, r0, coreParticles=[]):

        """Attracts a subset of particles to the core,
         repells the rest from the core"""

        attractForce = self.mm.CustomExternalForce(
            " COREk * ((COREr - CORErn) ^ 2)  ; "\
            "COREr = sqrt(x^2 + y^2 + COREtt^2)")
        attractForce.addGlobalParameter(
            "COREk", k * self.kT / (self.conlen * self.conlen))
        attractForce.addGlobalParameter("CORErn", r0 * self.conlen)
        attractForce.addGlobalParameter("COREtt", 0.001 * self.conlen)
        self.forceDict["CoreAttraction"] = attractForce
        for i in coreParticles:
            attractForce.addParticle(int(i), [])

        if r0 > 0.1:

            excludeForce = self.mm.CustomExternalForce(
                " CORE2k * ((CORE2r - CORE2rn) ^ 2) * step(CORE2rn - CORE2r) ;"
                "CORE2r = sqrt(x^2 + y^2 + CORE2tt^2)")
            excludeForce.addGlobalParameter("CORE2k", k *
                self.kT / (self.conlen * self.conlen))
            excludeForce.addGlobalParameter("CORE2rn", r0 * self.conlen)
            excludeForce.addGlobalParameter("CORE2tt", 0.001 * self.conlen)
            self.forceDict["CoreExclusion"] = excludeForce
            for i in xrange(self.N):
                excludeForce.addParticle(i, [])

    def fixParticlesZCoordinate(self, particles, zCoordinates, k=0.3,
                                useOtherAxis="z", mode="abs", gap=None):
        """Limits position of a set of particles in z coordinate

        Parameters
        ----------
        particles : list
            List of particles to be fixed.
        zCoordinates : list, or tuple of length 2
            If has length of particles, then should contain all Z coordinates
            If has length 2, then contains z coordinates of first and
            Nth particles, and the rest is approximated linearly.
        k : float, optional
            Strength of attraction, measured in kT/(bondlength)
        useOtherAxis : "x","y" or "z", optional
            Apply the same in the other dimension
        """

        if not len(particles) == len(zCoordinates):
            assert len(zCoordinates) == 2
            start, stop = tuple(zCoordinates)
            zCoordinates = []
            for par in particles:
                zCoordinates.append(start + float(
                    stop - start) * (par / float(self.N)))

        if (mode == "abs") and (gap is None):
            zFixForce = self.mm.CustomExternalForce(
            "ZFIXk * (sqrt((%s - ZFIXr0)^2 + ZFIXa^2) - ZFIXa)" % (
                                                           useOtherAxis,))
            zFixForce.addGlobalParameter("ZFIXk", k * self.kT / (self.conlen))
        elif (mode == "abs") and (gap is not None):
            zFixForce = self.mm.CustomExternalForce(
            "ZFIXk * step(%s - ZFIXr0 - ZFIXgap * 0.5) *"\
            " (sqrt((%s - ZFIXr0 - ZFIXgap * 0.5)^2 + ZFIXa^2) - ZFIXa) + "\
            "ZFIXk * step(-%s + ZFIXr0 - ZFIXgap * 0.5) * "\
            "(sqrt((-%s + ZFIXr0 - ZFIXgap * 0.5)^2 + ZFIXa^2) - ZFIXa)"\
            % (useOtherAxis, useOtherAxis, useOtherAxis, useOtherAxis))

            zFixForce.addGlobalParameter("ZFIXk", k * self.kT / (self.conlen))
            zFixForce.addGlobalParameter("ZFIXgap", self.conlen * gap)

        elif (mode == "quadratic") and (gap is None):
            zFixForce = self.mm.CustomExternalForce(
                "ZFIXk * ((%s - ZFIXr0)^2)" % (useOtherAxis,))
            zFixForce.addGlobalParameter("ZFIXk", k * self.kT /
                (self.conlen * self.conlen))
        elif (mode == "quadratic") and (gap is not None):
            zFixForce = self.mm.CustomExternalForce(
            "ZFIXk * (step(%s - ZFIXr0 - ZFIXgap * 0.5) * "\
            "(%s - ZFIXr0 - ZFIXgap * 0.5)^2 +  "\
            "step(-%s + ZFIXr0 - ZFIXgap * 0.5) * "\
            "(-%s + ZFIXr0 - ZFIXgap * 0.5)^2)" \
            % (useOtherAxis, useOtherAxis, useOtherAxis, useOtherAxis))

            zFixForce.addGlobalParameter("ZFIXk", k * self.kT /
                (self.conlen * self.conlen))
            zFixForce.addGlobalParameter("ZFIXgap", self.conlen * gap)
        else:
            raise RuntimeError("Not implemented")

        zFixForce.addPerParticleParameter("ZFIXr0")

        zFixForce.addGlobalParameter("ZFIXa", 0.05 * self.conlen)
        for par, zcoor in map(None, particles, zCoordinates):
            zFixForce.addParticle(int(par), [float(zcoor)])
        self.forceDict["fixZCoordinates"] = zFixForce


class ExperimentalSimulation(Simulation):
    "contains some experimental features"

    def quickLoad(self, data, mode="chain", Nchains=1,
                  trunc=None, confinementDensity="NoConfinement"):
        """quickly loads a set of repulsive chains,
        possibly adds spherical confinement"""
        self.setup()
        self.load(data)
        self.setLayout(mode, Nchains)
        self.addHarmonicPolymerBonds()
        self.addSimpleRepulsiveForce(trunc=trunc)
        if type(confinementDensity) != str:
            self.addSphericalConfinement(density=confinementDensity)

    def createWalls(self, left=None, right=None, k=0.5):
        "creates walls at x = left, x = right, x direction only"
        if left is None:
            left = self.data[0][0] + 1. * nm
        else:
            left = 1. * nm * left
        if right is None:
            right = self.data[-1][0] - 1. * nm
        else:
            right = 1. * nm * right

        if self.verbose == True:
            print "left wall created at ", left / (1. * nm)
            print "right wall created at ", right / (1. * nm)

        extforce2 = self.mm.CustomExternalForce(
            " WALLk * (sqrt((x - WALLright) * (x-WALLright) + WALLa * WALLa ) - WALLa) * step(x-WALLright) "
            "+ WALLk * (sqrt((x - WALLleft) * (x-WALLleft) + WALLa * WALLa ) - WALLa) * step(WALLleft - x) ")
        extforce2.addGlobalParameter("WALLk", k * self.kT / nm)
        extforce2.addGlobalParameter("WALLleft", left)
        extforce2.addGlobalParameter("WALLright", right)
        extforce2.addGlobalParameter("WALLa", 1 * nm)
        for i in xrange(self.N):
            extforce2.addParticle(i, [])
        self.forceDict["WALL Force"] = extforce2

    def addSphericalWell(self, r=10, depth=1):
        """pushes particles towards a boundary
        of a cylindrical well to create uniform well coverage"""

        extforce4 = self.mm.CustomExternalForce(
            "WELLdepth * (((sin((WELLr * 3.141592 * 0.5) / WELLwidth)) ^ 10)  -1) * step(-WELLr + WELLwidth);"
            "WELLr = sqrt(x^2 + y^2 + z^2 + WELLtt^2)")
        self.forceDict["Well attraction"] = extforce4

        #adding all the particles on which force acts
        for i in xrange(self.N):
            if self.domains[i] > 0.5:
                extforce4.addParticle(i, [])
        if r is None:
            try:
                r = self.sphericalConfinementRadius * 0.5
            except:
                exit("No spherical confinement radius defined yet."\
                     " Apply spherical confinement first!")
        if self.verbose == True:
            print "Well attraction added with r = %d" % r

        #assigning parameters of the force
        extforce4.addGlobalParameter("WELLwidth", r * nm)
        extforce4.addGlobalParameter("WELLdepth", depth * self.kT)
        extforce4.addGlobalParameter("WELLtt", 0.01 * nm)


class YeastSimulation(Simulation):
    """
    This class is maintained by Geoff to do simulations for the Yeast project
    """

    def addNucleolus(self, k=1, r=None):
        "method"
        if r is None:
            r = self.sphericalConfinementRadius

        extforce3 = self.mm.CustomExternalForce(
            "step(r-NUCaa) * NUCkb * (sqrt((r-NUCaa)*(r-NUCaa) + NUCt*NUCt) - NUCt);"
            "r = sqrt(x^2 + y^2 + (z + NUCoffset )^2 + NUCtt^2)")

        self.forceDict["NucleolusConfinement"] = extforce3
        #adding all the particles on which force acts
        if self.verbose == True:
            print "NUCleolus confinement from radius = %lf" % r
        #assigning parameters of the force
        extforce3.addGlobalParameter("NUCkb", k * self.kT / nm)
        extforce3.addGlobalParameter("NUCaa", (r - 1. / k) * nm * 1.75)
        extforce3.addGlobalParameter("NUCoffset", (r - 1. / k) * nm * 1.1)
        extforce3.addGlobalParameter("NUCt", (1. / k) * nm / 10.)
        extforce3.addGlobalParameter("NUCtt", 0.01 * nm)
        for i in xrange(self.N):
            extforce3.addParticle(i, [])

    def addLaminaAttraction(self, width=1, depth=1, r=None, particles=None):
        extforce3 = self.mm.CustomExternalForce(
            "-1 * step(LAMr-LAMaa + LAMwidth) * step(LAMaa + LAMwidth - LAMr) * LAMdepth"
            "* abs( (LAMr-LAMaa + LAMwidth) * (LAMaa + LAMwidth - LAMr)) / (LAMwidth * LAMwidth);"
            "LAMr = sqrt(x^2 + y^2 + z^2 + LAMtt^2)")
        self.forceDict["Lamina attraction"] = extforce3

        #re-defines lamina attraction based on particle index instead of domains.

        #adding all the particles on which force acts
        if particles is None:
            for i in xrange(self.N):
                extforce3.addParticle(i, [])
                if self.verbose == True:
                    print "particle %d laminated! " % i

        else:
            for i in particles:
                extforce3.addParticle(i, [])
                if self.verbose == True:
                    print "particle %d laminated! " % i

        if r is None:
            try:
                r = self.sphericalConfinementRadius
            except:
                exit("No spherical confinement radius defined yet."\
                     "Apply spherical confinement first!")

        if self.verbose == True:
            print "Lamina attraction added with r = %d" % r

        #assigning parameters of the force
        extforce3.addGlobalParameter("LAMaa", r * nm)
        extforce3.addGlobalParameter("LAMwidth", width * nm)
        extforce3.addGlobalParameter("LAMdepth", depth * self.kT)
        extforce3.addGlobalParameter("LAMtt", 0.01 * nm)


class GrandeSimulation(Simulation,
                       SimulationWithCrosslinks,
                       ExperimentalSimulation):
    pass
