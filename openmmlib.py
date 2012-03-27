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

Implemented forces 
------------------

All forces of the system are contained in the self.ForceDict dictionary. 
After the force is added, different methods are free to modify parameters of the force. 
Once the system is started, all forces are automatically applied and cannot be modified.   

Two types of bond forces are harmonic bond force and FENE-type bond as described in Grosberg papers.  
Individual bonds can be added using :py:func:`addBond <Simulation.addBond>`, while polymer bonds can 
be added using :py:func:`addHarmonicPolymerBonds <Simulation.addHarmonicPolymerBonds>`, etc. 

Repulsive force can be of 3 types: simple repulsife force U = 1/r^12; 
Grosberg repulsive force - a faster and better implementation of repulsive force; and 
LennardJones Force, that can be attractive and allows to specify extra
attraction/repulsion between any pairs of particles.

Stiffness force can be harmonic, or "special" Grosberg force, kept only for compatibility 
with the systems used in Grosberg forces

External forces include spherical confinement, cylindrical confinement, 
attraction to "lamina"- surface of a cylinder, gravity, etc. 


Functions
---------
Functions depend on each other, and have to be applied in certain groups  

1. 


:py:func:`load <Simulation.load>`  ---    Mandatory

:py:func:`setup <Simulation.setup>`  ---   Mandatory

:py:func:`setLayout <Simulation.setLayout>` --- Mandatory, after self.load()

:py:func:`saveFolder <Simulation.saveFolder>`  ---  Optional (default is folder with the code) 

2. 
 
self.add___Force()  --- Use any set of forces

Then go all the addBond, and other modifications and tweaks of forces. 

3. 

Before running actual simulation, it is advised to resolve all possible conflict by 
doing :py:func:`energyMinimization <Simulation.energyMinimization>`

4. 
 

:py:func:`doBlock <Simulation.doBlock>`  --- the actual simulation
 
:py:func:`save <Simulation.save>` --- saves conformation
  
  
Frequently-used settings - where to specify them? 
-------------------------------------------------

Select GPU ("0" or "1") - :py:func:`setup <Simulation.setup>`
 
Select chain/ring - :py:func:`setup <Simulation.setup>`
 
Select timestep or collision rate - :py:class:`Simulation`


-------------------------------------------------------------------------------

"""
# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php


import numpy
import scipy.stats
import cPickle
import sys,os
import time
import joblib
import tempfile
os.environ["LD_LIBRARY_PATH"] = "/usr/local/cuda/lib64:/usr/local/openmm/lib"

import simtk.openmm as openmm
import simtk.unit as units

nm = units.meter * 1e-9
fs = units.second * 1e-15
ps = units.second * 1e-12 



class Simulation():
    """Base class for openmm simulations

    """
    def __init__(self,timestep=80,thermostat=0.001,temperature = 300 * units.kelvin, 
                 verbose = False,
                 velocityReinitialize = True,  #reinitialize velocities at every block if E_kin is more than 2.4 
                  name = "sim"):   #name to print out
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
            Use it if you run simulations one after another and want to see what's going on. 
    
                
        """         
        
        self.name = name        
        self.timestep = timestep * fs
        self.collisionRate = thermostat * ps        
        self.temperature = temperature
        self.verbose = verbose
        self.velocityReinialize = velocityReinitialize
        self.loaded = False  #check if the data is loaded
        self.forcesApplied = False
        self.folder = "."         
        self.metadata = {}
        
    def setup(self,platform="OpenCL", PBC = False,PBCbox = None,GPU = "0",verbose = True):           
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
             
        verbose : bool, optional 
            Shout out loud about every change. 
          
        """
        self.step = 0
        if PBC == True: 
            self.metadata["PBC"] = True             

        self.kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA # Boltzmann constant
        self.kT = self.kB * self.temperature # thermal energy        
        self.mass = 100.0 * units.amu   #All masses are the same, changing them would be difficult in this formalism 
        self.bondsForException = []
        self.mm = openmm
        self.conlen = 1. * nm
        self.system = self.mm.System()
        self.PBC = PBC
        
        if self.PBC == True:   #if periodic boundary conditions
            if PBCbox == None:            
                data = self.getData()
                data -= numpy.min(data,axis = 0)
                self.setData(data)
                datasize = 1.1 * (2+(numpy.max(self.getData(),axis = 0) - numpy.min(self.getData(),axis = 0))) #size of the system plus some overhead                          
                self.SolventGridSize = (datasize / 1.1)  - 2                                                                
                print "density is ", self.N / (datasize[0] * datasize[1] * datasize[2])
            else:
                PBCbox = numpy.array(PBCbox)                                 
                datasize = PBCbox 
            
            self.metadata["PBCbox"] = PBCbox
            self.system.setDefaultPeriodicBoxVectors([datasize[0],0.,0.],[0.,datasize[1],0.],[0.,0.,datasize[2]])
            self.BoxSizeReal = datasize 
            print "system size is %lfx%lfx%lf nm, solvent grid size is   %lfx%lfx%lf in dimensionless units" % tuple(list(datasize) + list(self.SolventGridSize))
            time.sleep(5)
        
        
        self.GPU = GPU  #setting default GPU  

        if platform.lower() == "opencl":
            platformObject = openmm.Platform.getPlatformByName('OpenCL')
            platformObject.setPropertyDefaultValue('OpenCLDeviceIndex', self.GPU)
        elif platform.lower() == "reference":
            platformObject =  openmm.Platform.getPlatformByName('Reference')
        elif platform.lower() == "cuda":
            platformObject = openmm.Platform.getPlatformByName('Cuda')
            platformObject.setPropertyDefaultValue('CudaDevice', self.GPU)
        else:
            self.exit("\n!!!!!!!!!!unknown platform!!!!!!!!!!!!!!!\n")            
        self.platform = platformObject
                    
        self.forceDict = {}  #Dictionary to store forces
        try:
            for _ in xrange(self.N):
                self.system.addParticle(self.mass)
            print self.N, "particles loaded"
        except: pass                                 
        self.integrator = openmm.LangevinIntegrator(self.temperature, self.collisionRate, self.timestep)

    def saveFolder(self,folder):        
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
        

    def exitProgram(self,line):
        """Prints error line and exits program
        
        Parameters
        ----------
        
        line : str
            Line to print 
            
        """
        print line
        print "--------------> Bye <---------------"
        exit()
        
    def setLayout(self,mode = "chain",chains = None,Nchains = 1):
        """sets layout of chains for chains or rings. 
        By default makes one chain. You can change it to one ring (mode = ring).  
        You can either have chains/rings of equal length (mode =, Nchains =).   
        Or you can have chains or rings of different lengthes (mode =, chains =).  
        You can't have a mix of rings and chains
        
        Parameters
        ----------
        
        mode : "chain" or "ring" 
            Does the system consist of rings or chains? 
            
        chains : None or ((0,L1),(L1,L2),(L2,L3)...)
            Specifies exact chain/ring start/end positions if chains are of different lengths.
            
        Nchains : int 
            Number of chains, if they all are of the same lengths. 
            Ignored if chains is specified exactly.
            
        """
        N = self.N
        if mode in ["chain","ring"]:            
            if chains != None: 
                self.chains = chains           
            else: 
                self.chains = []
                for i in xrange(Nchains):
                    self.chains.append(((N* i)/Nchains,(N * (i+1))/Nchains))
        self.mode = mode
        #print self.N, chains
        layout = {"chains":chains,"mode":mode,"Nchains":Nchains}
        self.metadata["layout"] = layout 
    
    
    def getLayout(self):
        "returns configuration of chains"
        return self.chains            
    
    def load(self,filename,   #Input filename, or input data array   
             center = False    #Shift center of mass to zero? 
             ):        
        """loads data from file. 
        Accepts text files, joblib files or pure data as Nx3 or 3xN array
        
        Parameters
        ----------
        
        filename : joblib file, or text file name, or Nx3 or 3xN numpy array
            Input filename or array with data
        center : bool, optional 
            Move center of mass to zero before starting the simulation 
        """
        if type(filename) == str:             
            try:
                line0 = open(filename).readline() 
                try: N = int(line0)
                except ValueError: raise TypeError("Cannot read text file... reading pickle file")                         
                
                lines = open(filename).readlines()
                data = [[float(i) for i in j.split()] for j in lines if len(j) > 3]
                
                if len(data) != N:
                    raise ValueError("N does not correspond to the number of lines!")
            
            except TypeError:                
                mydict = dict(joblib.load(filename))
                data = mydict.pop("data") 
                self.oldMetadata = data
        else:
            data = filename
                    
        data = numpy.asarray(data,float)
        if len(data) == 3: 
            data = numpy.transpose(data)
        if len(data[0]) != 3: 
            self.exitProgram("strange data file")        
        if numpy.isnan(data).any(): self.exitProgram("\n!!!!!!!!!!file contains NANS!!!!!!!!!\n")
        if center == True:
            av = numpy.mean(data,0)
            data -= av
            
        if center == "zero":
            minvalue = numpy.min(data,0)
            data -= minvalue
                                        
        self.setData(data)                                     
        
        if self.verbose == True:
            print "center of mass is", numpy.mean(self.data,0)
            print "Radius of gyration is," , self.RG()
        
        try: 
            for i in xrange(self.N):
                self.system.addParticle(self.mass)
            if self.verbose == True:print "%d particles loaded" % self.N
            self.loaded = True
        except: pass



    def save(self,filename = None):        
        "Saves conformation plus some metadata. Metadata is not interpreted by this library, and is for your reference"
        self.metadata["data"] = self.getData()
        self.metadata["timestep"] = self.timestep / fs
        self.metadata["Collision rate"] = self.collisionRate / ps                 
        if filename == None: 
            f = os.path.join(self.folder , "block%d.dat" % self.step)
        else:
            f = os.path.join(self.folder , filename)
        joblib.dump(self.metadata,filename = f,compress = 3)
        
    def getData(self):
        "Returns an Nx3 array of positions"
        return numpy.asarray(self.data / nm, dtype = "float32") 
    
    def getScaledData(self):        
        """Returns data, scaled back to PBC box """
        alldata = self.getData()
        boxsize = numpy.array(self.BoxSizeReal)
        mults = numpy.array(alldata / boxsize[None,:],int)
        return alldata - mults * boxsize[None,:]
    
    def setData(self,data):
        """Sets particle positions
        
        Parameters
        ----------
        
        data : Nx3 array-line 
            Array of positions with distance ~1 between connected atoms. 
        """                 
        data = numpy.asarray(data,dtype = "float")
        self.data = units.Quantity(data,nm)
        self.N = len(self.data)    
            
    def RG(self):
        """
        Returns
        -------
        
        Gyration ratius in units of length (bondlength). 
        """
        return numpy.sqrt(numpy.sum(numpy.var(numpy.array(self.data/nm) ,0)))
    
    def RMAX(self):
        """
        Returns
        -------
        Distance to the furthest from the origin particle. 
        
        """
        return numpy.max(numpy.sqrt(numpy.sum((numpy.array(self.data/nm) )**2,1)))
        
    def dist(self,i,j):
        """
        Calculates distance between particles i and j
        """                
        dif = (self.data[i] - self.data[j]) / (nm)
        return numpy.sqrt(sum(dif**2))

    def useDomains(self,domains = None, filename = None):
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
        

        if domains != None:
            self.domains = domains            
            
        elif filename != None:
            self.domains = cPickle.load(open(domains))
        else: self.exit("You have to specify at least some domains!")
            
        if len(self.domains) != self.N: self.exitProgram("Wrong domain lengths")                        
        cPickle.dump(self.domains,open(os.path.join(self.folder ,"/domains.dat"),'wb'))
                    
                
    def _initHarmonicBondForce(self):
        "Internal, inits harmonic forse for polymer and non-polymer bonds"
        if "HarmonicBondForce" not in self.forceDict.keys(): 
            self.forceDict["HarmonicBondForce"] = self.mm.HarmonicBondForce()
        self.bondType = "Harmonic"        


    def _initGrosbergBondForce(self):        
        "inits Grosberg FENE bond force"
        if "GrosbergBondForce" not in self.forceDict.keys():
            force = "- 0.5 * GROSk * GROSr0 * GROSr0 * log(1-(r/GROSr0)* (r / GROSr0)) + (4 * GROSe * ((GROSs/r)^12 - (GROSs/r)^6) + GROSe) * step(GROScut - r)"
            bondforce3  = openmm.CustomBondForce(force)            
            bondforce3.addGlobalParameter("GROSk",30 * self.kT / (self.conlen * self.conlen))
            bondforce3.addGlobalParameter("GROSr0",self.conlen * 1.5)
            bondforce3.addGlobalParameter('GROSe',self.kT)
            bondforce3.addGlobalParameter('GROSs',self.conlen)
            bondforce3.addGlobalParameter("GROScut",self.conlen * 2.**(1./6.))
            self.forceDict["GrosbergBondForce"] = bondforce3
                   
    def addBond(self,
                i,j,   #particles connected by bond
                bondWiggleDistance,  # -----> Harmonic only!!! <--------  
                                     #Flexibility of the bond,  measured in distance at which energy equals kT
                distance = None,     #Equilibrium length of the bond -----> Harmonic only!!! <--------
                bondType = None,     #Harmonic or Grosberg  
                verbose = None):     #Set this to False if you're in verbose mode and don't want to contaminate output by 10000 messages
        """Adds bond between two particles, allows to specify parameters for Harmonic bonds only!
        
        Parameters
        ----------
        
        i,j : int 
            Particle numbers
        bondWiggleDistance : float
            Average displacement from the equilibrium bond distance
        bondType : "Harmonic" or "Grosberg" 
            Type of bond. Distance and bondWiggleDistance can be specified for harmonic bonds only
        verbose : bool 
            Set this to False if you're in verbose mode and don't want to print "bond added" message
            
        """
        
        
        if verbose == None:
            verbose = self.verbose 
        
        bondSize = float(bondWiggleDistance)
        if distance == None: distance = self.conlen / nm
        else:  distance = self.conlen * distance / nm
        
        distance  = float(distance)        
        if bondType == None: 
            bondType = self.bondType             
        if bondType.lower() == "harmonic":
            self._initHarmonicBondForce()
            kbond = ( 2 * self.kT / (bondSize * self.conlen) ** 2 ) / (units.kilojoule_per_mole / nm ** 2 )
            self.forceDict["HarmonicBondForce"].addBond(int(i),int(j),float(distance),float(kbond))
            if verbose == True: print "Harmonic bond added between %d,%d, params %lf %lf" % (i,j,float(distance),float(kbond))
        elif bondType.lower() == "grosberg":
            self.initGrosbergForce()
            self.forceDict["GrosbergForce"].addBond(int(i),int(j),[])
        else: self.exitProgram("Bond type not known")

          
    def addHarmonicPolymerBonds(self,wiggleDist = 0.05):
        """Adds harmonic bonds connecting polymer chains
        wiggleDist controls the distance at which energy of the bond equals kT""" 
              
        for i in self.chains:
            for j in xrange(i[0],i[1] - 1):
                self.addBond(j,j+1,wiggleDist,distance = 1, bondType = "Harmonic",verbose = False)
                self.bondsForException.append((j,j+1))

            if self.mode == "ring":
                self.addBond(i[0],i[1] -1,wiggleDist,distance = 1, bondType = "Harmonic")
                self.bondsForException.append((i[0],i[1] -1))            
                if self.verbose == True: print "ring bond added", i[0],i[1]-1
                
        self.metadata["HarmonicPolymerBonds"] = {"wiggleDist":wiggleDist}
    
    def addGrosbergPolymerBonds(self,k = 30):
        """Adds FENE bonds according to Grosberg paper. 
        This method has a repulsive potential build-in, so that Grosberg bonds could be used with truncated potentials. 
        Is of no use unless you really need to simulate Grosberg-type system.
        
        Parameters
        ----------
        k : float, optional
            Arbitrary parameter; default value as in Grosberg paper.   
        
         """
        self._initGrosbergBondForce()
        force = self.forceDict["GrosbergBondForce"]
        for i in self.chains:
            for j in xrange(i[0],i[1] - 1):
                force.addBond(j,j+1,[])
                self.bondsForException.append((j,j+1))               
            if self.mode == "ring":
                force.addBond(i[0],i[1] - 1,[])
                self.bondsForException.append((i[0],i[1] - 1))
                if self.verbose == True: print "ring bond added", i[0],i[1]-1
        self.metadata["GorsbergPolymerForce"] = {"k":k}
                        

    def addStiffness(self,k = 40):
        """Adds harmonic angle bonds. k specifies energy in kT at one radian
        
        Parameters
        ----------
        
        k : float
            Potential is k * alpha^2 * 0.5 * kT
        """
        myforce= self.mm.CustomAngleForce("k * (theta - 3.141592) * (theta - 3.141592) * (0.5)")
        self.forceDict["AngleForce"] = myforce 
        for i in self.chains:
            for j in xrange(i[0]+1,i[1] - 1):
                myforce.addAngle(j-1,j,j+1,[])
        myforce.addGlobalParameter("k",k*self.kT)
        self.metadata["AngleForce"] = {"stiffness":k}

    def addGrosbergStiffness(self, k = 1.5):
        """Adds stiffness according to the Grosberg paper. Parameters are synchronized with normal stiffness
        
        Parameters
        ----------
        
        k : float
            Synchronized with regular stiffness. Default value is very flexible, as in Grosberg paper. Default value maximizes entanglement length.  
            
        """

        myforce= self.mm.CustomAngleForce("k * (1 - cos(theta - 3.141592))")
        self.forceDict["AngleForce"] = myforce 
        for i in self.chains:
            for j in xrange(i[0]+1,i[1] - 1):
                myforce.addAngle(j-1,j,j+1,[])
        myforce.addGlobalParameter("k",k * self.kT)
        self.metadata["GrosbergAngleForce"] = {"stiffness":k}

        
    def addSimpleRepulsiveForce(self,cutoff = 1.7,trunc = None,rep = 0.26):
        """Creates a repulsive force between all particles. 
        
        .. warning:: 
            This force is about to be deprecated. GrosbergRepulsiveForce is more efficient and equally powerful. 
    
        
        Parameters
        ----------
        
        cutoff : float, optional, default value is good. 
            Cutoff distance. Small values are not adviced. 
            
        trunc : float or None, optional
            Cutoff energy, used to allow chain passing. Measured in kT. Value of 2.5 yields frequent passing, 4 - average passing, 6 - rare passing.  
        rep : float, optional
            Strength of repulsive potential : U = rep * 1/r^12. Default value is good. 
         
        """
        self.metadata["SimgleRepulsiveForce"] = {"cutoff":cutoff,"trunc":trunc,"rep":rep}
        nbCutOffDist = self.conlen * cutoff    #repulsive part saturates quickly 
        if trunc == None: repul_energy = 'REPepsilon*(REPsigma/r)^12'
        else: repul_energy = '1/((1/REPcutoff) + (1/REPU + 0.0001 * REPcutoff));REPU=REPepsilon*(REPsigma/r)^12;r2 = (r^10. + (0.3 * REPs)^10.)^0.1 '
        #last equation is to avoid NANs when r is close to zero  
        
        epsilonRep = rep * units.kilocalorie_per_mole
        sigmaRep = 1.06 * self.conlen
        self.forceDict["Nonbonded"] = openmm.CustomNonbondedForce(repul_energy)
        repulforce = self.forceDict["Nonbonded"]
        if trunc!=None: repulforce.addGlobalParameter('REPcutoff',trunc * self.kT)
        repulforce.addGlobalParameter('REPepsilon',epsilonRep)
        repulforce.addGlobalParameter('REPsigma',sigmaRep) 
        for _ in range(self.N):
            repulforce.addParticle(())                
        if self.PBC: 
            repulforce.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
            print "Using periodic boundary conditions!!!!!!!!!!"            
        else: repulforce.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic)
        repulforce.setCutoffDistance(nbCutOffDist)
    
    def addGrosbergRepulsiveForce(self,trunc=None):
        """This is the fastest repulsive force.
        
        Parameters
        ----------
        
        trunc : None or float
             truncation energy in kT, used for chain crossing.
             Value of 2.5 yields frequent passing, 4 - average passing, 6 - rare passing. 
        
        """
        self.metadata["GrosbergRepulsiveForce"] = {"trunc":trunc}            
        nbCutOffDist = self.conlen * 2.**(1./6.)
        if trunc == None:     
            repul_energy = "4 * REPe * ((REPs/r)^12 - (REPs/r)^6) + REPe"
        else: 
            repul_energy = "1 / (1 / REPcut + 1 / (REPU0 + REPa * REPcut) ) - REPcut * (REPa/(REPa+1)) ;REPU0 = 4 * REPe * ((REPs/r2)^12 - (REPs/r2)^6) + REPe;r2 = (r^10. + (0.3 * REPs)^10.)^0.1"            
        self.forceDict["Nonbonded"] = openmm.CustomNonbondedForce(repul_energy)
        repulforce = self.forceDict["Nonbonded"]
        repulforce.addGlobalParameter('REPe',self.kT)
        repulforce.addGlobalParameter('REPs',self.conlen)
        if trunc != None: 
            repulforce.addGlobalParameter('REPcut',self.kT * trunc)        
            repulforce.addGlobalParameter('REPa',0.001)
        for _ in range(self.N):
            repulforce.addParticle(())                
        if self.PBC: 
            repulforce.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
            print "Using periodic boundary conditions!!!!!!!!11"            
        else: repulforce.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic)
        repulforce.setCutoffDistance(nbCutOffDist)
        
    def addLennardJonesForce(self,cutoff = 2.5, 
                            domains = False,
                            epsilonRep = 0.24,epsilonAttr = 0.27,   #parameters for LJ force 
                            blindFraction = -1,
                            sigmaRep = None, sigmaAttr = None):
        """
        Adds a lennard-jones force, that allows for mutual attraction. This is the slowest force out of all repulsive. 
        
        .. note:: 
            This is the only force that allows for so-called "exceptions'. 
            Exceptions allow you to change parameters of the force for a specific pair of particles. 
            This can be used to create short-range attraction between pairs of particles.
            See manual for Openmm.NonbondedForce.addException.  
        
        Parameters
        ----------
        
        cutoff : float, optional 
            Cutoff value. Default is good. 
        domains : bool, optional 
            Use domains, defined by :py:func:'setDomains <Simulation.setDomains>'
        epsilonRep : float, optional  
            Epsilon (attraction strength) for LJ-force for all particles (except for domain) 
        epsilonAttr : float, optional
            Epsilon for attractive domain (if domains are used)
        blindFraction : float, 0<x<1
            Fraction of particles that are "transparent" - used here instead of truncation 
        sigmaRep, sigmaAttr: float, optional
            Radius of particles in the LJ force. For advanced fine-tuning. 
        
         """
        self.metadata["LennardJonesForce"] = {"cutoff":cutoff,"domains":domains,"epsilonRep":epsilonRep, "epsilonAttr":epsilonAttr,"blindFraction":blindFraction,
                                              "sigmaRep":sigmaRep, "sigmaAttr":sigmaAttr}
        if blindFraction > 0.99: self.exitProgram ("why do you need this force without particles??? set blindFraction between 0 and 1") 
        if (sigmaRep == None) and (sigmaAttr == None):
            sigmaAttr = sigmaRep = self.conlen / nm
        else:
            sigmaAttr  = sigmaAttr * self.conlen  / nm
            sigmaRep = sigmaRep * self.conlen  / nm                    
        epsilonRep = epsilonRep *  units.kilocalorie_per_mole / units.kilojoule_per_mole
        epsilonAttr = epsilonAttr * units.kilocalorie_per_mole / units.kilojoule_per_mole
        nbCutOffDist = self.conlen * cutoff
        self.epsilonRep = epsilonRep                 
        repulforce = openmm.NonbondedForce()
        
        self.forceDict["Nonbonded"] = repulforce
        if domains == False:
            for i in xrange(self.N):
                if numpy.random.random() > blindFraction: 
                    repulforce.addParticle(0,sigmaRep,epsilonRep)
                else:
                    repulforce.addParticle(0,0,0)
        else:
            for i in xrange(self.N):
                if numpy.random.random() > blindFraction:
                    if self.domains[i] == 0:
                        repulforce.addParticle(0,sigmaRep,epsilonRep)
                    else:
                        repulforce.addParticle(0,sigmaAttr,epsilonAttr)
                else:
                    repulforce.addParticle(0,0,0)
        if self.PBC: 
            repulforce.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
            print "Using periodic boundary conditions!!!!"
        else: repulforce.setNonbondedMethod(openmm.NonbondedForce.CutoffNonPeriodic)
        repulforce.setCutoffDistance(nbCutOffDist)
        
    def addMutualException(self, particles):        
        """used to exclude a bunch of particles from calculation of nonbonded force
        
        Parameters
        ----------
        particles : list 
            List of particles for whom to exclude nonbonded force. 
        """
        for i in particles: #xrange(len(particles)):
            for j in particles:#xrange(len(particles)):
                if j> i:        
                    self.bondsForException.append((i,j))

        
    def addInteraction(self,i,j,epsilon,sigma = None, length = 3):
        """Adds attractive short-range interaction of strength epsilon between particles i,j and a few neighboring particles
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
        if not self.metadata.has_key("interactions"): 
            self.metadata["interactions"] = []
        self.metadata["interactions"].append((i,j))
        
        if type(self.forceDict["Nonbonded"]) != self.mm.NonbondedForce:
            self.exit("Cannot add interactions without Lennard-Jones nonbonded force")
            
        if sigma == None: sigma = 1.1 * self.conlen
        epsilon = epsilon * units.kilocalorie_per_mole
        if (min(i,j) < length) or (max(i,j) > self.N - length): 
            print "!!!!!!!!!bond with %d and %d is out of range!!!!!" % (i,j)
            return        
        repulforce = self.forceDict["Nonbonded"]
        for t1 in xrange(i-length/2, i+(length - length/2)):
            for t2 in xrange(j-length/2, j+(length - length/2 )):
                repulforce.addException(t1,t2,0,sigma,epsilon,True)
                if self.verbose == True: 
                    print "Exception added between particles %d and %d" % (t1,t2)
                    
        for tt in xrange(i-length, i+length):
            repulforce.setParticleParameters(tt,0,self.conlen, self.epsilonRep)
        for tt in xrange(j-length, j+length):
            repulforce.setParticleParameters(tt,0,self.conlen, self.epsilonRep)

    def addCylindricalConfinement(self,r,bottom = True,k=0.1,weired = False):
        "As it says. Weird was used for Natasha simulations... and is weird."
        
        self.metadata["CylindricalConfinement"] = {"r":r,"bottom":bottom,"k":k,"weird":weired}
        if bottom == True: extforce2 = openmm.CustomExternalForce("step(r-CYLaa) * CYLkb * (sqrt((r-CYLaa)*(r-CYLaa) + CYLt*CYLt) - CYLt) + step(-z) * CYLkb * (sqrt(z^2 + CYLt^2) - CYLt) ;r = sqrt(x^2 + y^2 + CYLtt^2)")
        else: extforce2 = openmm.CustomExternalForce("step(r-CYLaa) * CYLkb * (sqrt((r-CYLaa)*(r-CYLaa) + CYLt*CYLt) - CYLt) ;r = sqrt(x^2 + y^2 + CYLtt^2)") 
        if weired == True: extforce2 = openmm.CustomExternalForce(" 0.6 * CYLkt * CYLaa*CYLaa / (CYLaa * CYLaa + r * r) + step(r-CYLaa) * CYLkb * (sqrt((r-CYLaa)*(r-CYLaa)  + CYLt*CYLt) - CYLt) + step(-z) * CYLkb * (sqrt(z^2 + CYLt^2) - CYLt) ;r = sqrt(x^2 + y^2 + CYLtt^2)")

        self.forceDict["CylindricalConfinement"] = extforce2
        for i in xrange(self.N): extforce2.addParticle(i,[])
        extforce2.addGlobalParameter("CYLkb",k*self.kT/nm);
        extforce2.addGlobalParameter("CYLkt",self.kT);
        extforce2.addGlobalParameter("CYLweired",nm)
        extforce2.addGlobalParameter("CYLaa",(r - 1./k)*nm);
        extforce2.addGlobalParameter("CYLt", (1./(10*k))*nm);
        extforce2.addGlobalParameter("CYLtt",0.01*nm);
    def addSphericalConfinement(self,r="density",   #radius... by default uses certain density  
                                k = 5.,             #How steep the walls are
                                density = .3):      #target density, measured in particles per cubic nanometer (bond size is 1 nm) 
        """Constrain particles to be within a sphere. With no parameters creates sphere with density .3 
        
        Parameters
        ----------
        r : float or "density", optional
            Radius of confining sphere. If "density" requires density, or assumes density = .3 
        k : float, optional 
            Steepness of the confining potential, in kT/nm 
        density : float, optional, <1 
            Density for autodetection of confining radius. 
            Density is calculated in particles per nm^3, i.e. at density 1 each sphere has a 1x1x1 cube.             
        """
        self.metadata["SphericalConfinement"] = {"r":r,"k":k,"density":density}        

        extforce2 = openmm.CustomExternalForce("step(r-SPHaa) * SPHkb * (sqrt((r-SPHaa)*(r-SPHaa) + SPHt*SPHt) - SPHt) ;r = sqrt(x^2 + y^2 + z^2 + SPHtt^2)")
        self.forceDict["SphericalConfinement"] = extforce2
        
        for i in xrange(self.N): extforce2.addParticle(i,[])
        if density != None:
            r = (3 * self.N / (4 * 3.141592 * density)) ** (1/3.)
             
        self.sphericalConfinementRadius = r
        if self.verbose == True: print "Spherical confinement with radius = %lf" % r 
        #assigning parameters of the force 
        extforce2.addGlobalParameter("SPHkb",k*self.kT/nm);
        extforce2.addGlobalParameter("SPHaa",(r - 1./k)*nm);
        extforce2.addGlobalParameter("SPHt",(1./k)*nm/10.);
        extforce2.addGlobalParameter("SPHtt",0.01*nm);

    def addLaminaAttraction(self,width = 1,depth = 1, r = None):
        """Attracts one domain to the lamina. Infers radius from spherical confinement, that has to be initialized already.
        
        Parameters
        ----------
        
        width : float, optional 
            Width of attractive layer next to the lamina, nm.  
        depth : float, optional
            Depth of attractive potential in kT 
        r : float, optional 
            Radius of an attractive cage. If not specified, inferred from previously defined spherical potential. 
        """
        
        self.metadata["laminaAttraction"] = {"width":width,"depth":depth,"r":r}
        extforce3 = openmm.CustomExternalForce("step(LAMr-LAMaa + LAMwidth) * step(LAMaa + LAMwidth - LAMr) * LAMdepth * (LAMr-LAMaa + LAMwidth) * (LAMaa + LAMwidth - LAMr) / (LAMwidth * LAMwidth)  ;LAMr = sqrt(x^2 + y^2 + z^2 + LAMtt^2)")
        self.forceDict["Lamina attraction"] = extforce3
        
        #adding all the particles on which force acts
        for i in xrange(self.N): 
            if self.domains[i] > 0.5: extforce3.addParticle(i,[])
        if r == None:
            try: r = self.sphericalConfinementRadius
            except: raise ValueError("No spherical confinement radius defined yet. Apply spherical confinement first!")
        if self.verbose == True: print "Lamina attraction added with r = %d"% r 
         
        extforce3.addGlobalParameter("LAMaa", r * nm);
        extforce3.addGlobalParameter("LAMwidth",width*nm);
        extforce3.addGlobalParameter("LAMdepth",depth * self.kT);        
        extforce3.addGlobalParameter("LAMtt",0.01*nm);

    
    def tetherParticles(self,particles, k = 30):        
        """tethers particles in the 'particles' array. Increase k to tether them stronger, but watch the system!
        
        Parameters
        ----------
        
        particles : list of ints 
            List of particles to be tethered (fixed in space) 
        k : int, optional 
            Steepness of the tethering potential. Values >30 will require decreasing potential, 
            but will make tethering rock solid. 
        """
        self.metadata["TetheredParticles"] = {"particles":particles,"k":k}
         
        extforce2 = openmm.CustomExternalForce(" TETHkb * ((x - TETHx0)^2 + (y - TETHy0)^2 + (z - TETHz0)^2)")
        self.forceDict["Tethering Force"] = extforce2

        #assigning parameters of the force 
        extforce2.addGlobalParameter("TETHkb",k*self.kT/nm);
        extforce2.addPerParticleParameter("TETHx0")
        extforce2.addPerParticleParameter("TETHy0")
        extforce2.addPerParticleParameter("TETHz0")
        for i in particles:#adding all the particles on which force acts
            coordinates = self.data[i]
            extforce2.addParticle(i,list(coordinates))
            if self.verbose == True: print "particle %d tethered! " % i 
        

    def addGravity(self,k = 0.1,cutoff = None):        
        """adds force pulling downwards in z direction
        When using cutoff, acts only when z>cutoff"""
        self.metadata["gravity"] = {"k":k,"cutoff":cutoff}        
        if cutoff == None:
            extforce3 = openmm.CustomExternalForce("kG * z")
        else: 
            extforce3 = openmm.CustomExternalForce("kG * (z - cutoffG) * step(z - cutoffG)")
            extforce3.addGlobalParameter("cutoffG", cutoff * nm)
        extforce3.addGlobalParameter("kG",k * self.kT / (nm))
            
        for i in xrange(self.N): extforce3.addParticle(i,[])
        self.forceDict["Gravity"] = extforce3
                                  
    def _applyForces(self):
        """Applies all the forces in the forcedict. 
        Forces should not be modified after that, unless you do it carefully (see openmm reference)."""
        exc = self.bondsForException        
        print "Number of exceptions:", len(exc)        
        
        if len(exc) > 0:
            exc = numpy.array(exc) 
            exc = numpy.sort(exc,axis = 1) 
            exc = [tuple(i) for i in exc]
            exc = list(set(exc)) #only unique pairs are left
             
        for i in self.forceDict.keys():   #Adding exceptions 
            force = self.forceDict[i]
            if type(force) == openmm.NonbondedForce:
                for pair in exc: 
                    force.addException(int(pair[0]),int(pair[1]),0,0,0,True)
            elif type(force) == openmm.CustomNonbondedForce:
                for pair in exc:
                    #force.addExclusion(*pair) 
                    force.addExclusion(int(pair[0]),int(pair[1]))                                     
            print "adding force ",i,self.system.addForce(self.forceDict[i])        
        self.context = openmm.Context(self.system, self.integrator, self.platform)
        self.forcesApplied = True
        
    def initVelocities(self,mult=1):
        """Initializes particles velocities
        
        Parameters 
        ----------
        mult : float, optional
            Multiply velosities by this. Is good for a cold/hot start. 
        """
        try: self.context
        except: raise ValueError("No context, cannot set velocs. Initialize context before that") 
        
        sigma = units.sqrt(self.kT/ self.system.getParticleMass(0) )   #calculating mean velocity        
        velocs = units.Quantity(mult * numpy.random.normal(size=(self.N,3)), units.meter) * (sigma / units.meter)
        #Guide to simtk.unit: 1. Always use units.quantity. 2. Avoid dimensionless shit. 3. If you have to, create fake units, as done here with meters                          
        self.context.setVelocities(velocs)
        
    def initPositions(self):
        """Sends particle coordinates to OpenMM system. 
        If system has exploded, this is used in the code to reset coordinates. """
        self.context.setPositions(self.data)
        state = self.context.getState(getPositions=True,getEnergy = True)  #get state of a system: positions, energies
        eP = state.getPotentialEnergy()/self.N/self.kT        
        print "Positions loaded, potential energy is %lf" % eP
        
    def reinitialize(self,mult = 1):
        """Reinitializes the OpenMM context object. 
        This should be called if low-level parameters, such as timestep, thermostat or forces, have changed.          
        
        Parameters
        ----------
        mult : float, optional 
            mult to be passed to :py:func:'initVelocities <Simulation.initVelocities>'
        """
        self.context.reinitialize()
        self.initPositions()
        self.initVelocities(mult)
        
    def energyMinimization(self,steps = 1000,twoStage = False,collisionRate = None):
        """Runs system at smaller timestep and higher collision rate to resolve possible conflicts.
        Does 10 or 15 (two-stage) blocks. 
        
        Parameters
        ----------
        steps : int, optional 
            Number of steps to make
        twoStage : bool, optional 
            Do additional minimization with low thermostat to better reolve conflicts. 
        collisionRage : bool, optional 
            Override default collision rage. Can be used for heavy minimization. 
            
        """
        print "Performing energy minimization"
        if self.forcesApplied == False:
            self._applyForces()
            self.forcesApplied = True
        def_step = self.integrator.getStepSize()
        self.integrator.setStepSize(def_step/15)
        def_fric = self.integrator.getFriction()
        if collisionRate == None: self.integrator.setFriction(1.3)
        else: self.integrator.setFriction(collisionRate)
        self.reinitialize()
        
        for _ in xrange(10): 
            self.doBlock(steps = steps, increment = False)
            self.initVelocities()
            
        if twoStage == True:
            self.integrator.setFriction(0.1)
            self.reinitialize()
            for _ in xrange(5):            
                self.doBlock(steps = steps, increment = False)
                self.initVelocities()

        self.integrator.setFriction(def_fric)
        self.integrator.setStepSize(def_step)
        self.reinitialize()
        print "Finished energy minimization"
                
    def doBlock(self,steps = None,increment = True,num = None):
        """performs one block of simulations, doing steps timesteps, or steps_per_block if not specified. 
        
        Parameters
        ----------
        
        steps : int or None 
            Number of timesteps to perform. If not specified, tries to infer it from self.steps_per_block
        increment : bool, optional
            If true, will not increment self.steps counter
        num : int or None, optional 
            If specified, will split the block in num subblocks. 
            Default value is 10. 
            
        
        """
        
        if self.forcesApplied == False: 
            self._applyForces()
            self.forcesApplied = True
        if increment == True: self.step += 1
        if steps == None: steps = self.steps_per_block
        if (increment == True) and ((self.step % 50) == 25): self.printStats()
        
        for attempt in xrange(6):
            print "%s  block=%d" % (self.name,self.step),
            if num == None:
                num = steps/10 + 1
            a = time.time()
            for _ in xrange(steps/num):
                self.integrator.step(num)  #integrate!
                print ".",
                sys.stdout.flush()
            self.integrator.step(steps%num)
            
            self.state = self.context.getState(getPositions=True,getEnergy = True)  #get state of a system: positions, energies
            b = time.time()
            coords = self.state.getPositions(asNumpy=True)
            newcoords = coords / nm
            eK = self.state.getKineticEnergy()/self.N/self.kT                          #calculate energies in KT/particle
            eP = self.state.getPotentialEnergy()/self.N/self.kT
            if self.velocityReinialize == True:
                if eK > 2.4:
                    self.initVelocities()                       
            print " Coord[1]=[%.1lf %.1lf %.1lf] " % tuple(newcoords[0]),            
            
            if (numpy.isnan(newcoords).any()) or (eK > 20) or (numpy.isnan(eK) ) or (numpy.isnan(eP)):
                self.context.setPositions(self.data)
                self.initVelocities()
                print "trying one more time at step # %i" % self.step
            else:
                dif = numpy.sqrt(numpy.mean(numpy.sum((newcoords - self.getData())**2,axis = 1)))
                print "sq(MSD)=%.2lf" % (dif,), 
                self.data = coords
                print " %.2lf kin, %.2lf pot, %.2lf tot," % (eK,eP,eK+eP),  " Rg=%.3lf" % self.RG(),
                print "SPS=%.0lf:"%(steps/(float(b-a)))                
                break
            if attempt in [3,4]:
                self.energy_minimization(10)
            if attempt == 5:                                
                self.exitProgram("exceeded number of attmpts")
                
    def printStats(self):
        """Prints detailed statistics of a system. 
        Will be run every 50 steps        
        """
        ss = scipy.stats #shortcut
        state = self.context.getState(getPositions = True, getVelocities = True, getEnergy = True)
        
        eP = state.getPotentialEnergy()
        pos = numpy.array(state.getPositions()/nm)
        bonds = numpy.sum(numpy.diff(pos,axis = 0) **2, axis = 1)
        sbonds = numpy.sort(bonds)
        vel = state.getVelocities()         
        mass = self.system.getParticleMass(0)
        vkT = numpy.array(vel / units.sqrt(self.kT / mass),dtype = float)
        EkPerParticle = 0.5 * numpy.sum(vkT**2,axis = 1)
        EkSimulated = 0.5 * numpy.sum((numpy.random.randn(*vkT.shape)) **2, axis = 1)
        
        cm = numpy.mean(pos,axis = 0)
        centredPos = pos - cm[None,:]
        dists = numpy.sqrt(numpy.sum(centredPos**2, axis = 1))
        per95 = numpy.percentile(dists, 95)
        den = (0.95 * self.N) / ((4. * numpy.pi * per95**3) / 3 )
        
        print
        print "Statistics for the simulation %s, number of particles: %d,  number of chains: %d,  mode:  %s" % (self.name, self.N, len(self.chains),self.mode)
        print
        print "Mean position is: ", numpy.mean(pos,axis = 0), "  Rg = ",self.RG()
        print "     mean bond size is ", numpy.mean(bonds)
        print "     three shortest/longest bonds are ", sbonds[:3],"  ",sbonds[-3:]
        print "     95 percentile of distance to center is:   ",per95
        print "     density of closest 95% monomers is:   ", den        
        print
        print "Statistics for velocities:"      
        print "     mean kinetic energy is: ", numpy.mean(EkPerParticle), "should be:", 1.5 
        print "     fastest particles are: ", numpy.sort(EkPerParticle)[-5:]
        print "     fastest particles in simulated distribution:   ", numpy.sort(EkSimulated)[-5:]
        print "     moments of velosities distribution are:", numpy.mean(EkPerParticle), numpy.var(EkPerParticle), ss.skew(EkPerParticle), scipy.stats.kurtosis(EkPerParticle)
        print "     moments of simulated  distribution are:", numpy.mean(EkSimulated), numpy.var(EkSimulated), ss.skew(EkSimulated), scipy.stats.kurtosis(EkSimulated)
        print
        print "Statistics for the system:"
        print "     Forces are: ", self.forceDict.keys()
        print "     Number of exceptions:  ", len(self.bondsForException)
        print 
        print "Potential Energy Ep = ", eP/self.N/self.kT        
        
    def show(self):
        """shows system in rasmol by drawing spheres
        draws 4 spheres in between any two points (5 * N spheres total)
        """
              
        #if you want to change positions of the spheres along each segment, change these numbers
        #e.g. [0,.1, .2 ...  .9] will draw 10 spheres, and this will look better
        shifts = [0.,0.2,0.4,0.6,0.8]
        data = self.getData()        
        if len(data[0]) != 3: 
            data = numpy.transpose(data)
        if len(data[0]) != 3:
            print "wrong data!"
            return        
        #determining the 95 percentile distance between particles,  
        meandist = numpy.percentile(numpy.sqrt(numpy.sum(numpy.diff(data,axis = 0)**2,axis = 1)),95)
        #rescaling the data, so that bonds are of the order of 1. This is because rasmol spheres are of the fixed diameter. 
        data /= meandist                
        rascript = tempfile.NamedTemporaryFile() #writing the rasmol script. Spacefill controls radius of the sphere.
        rascript.write("""wireframe off 
        color temperature
        spacefill 100 
        background white
        """)
        rascript.flush()        
        
        #creating the array, linearly chanhing from -225 to 225, to serve as an array of colors  
        colors = numpy.array([int((j*450.)/(len(data)))-225 for j in xrange(len(data))])    
        
        #creating spheres along the trajectory            
        newData = numpy.zeros((len(data) * len(shifts) - (len(shifts) - 1) ,4))  
        for i in xrange(len(shifts)):            
            newData[i:-1:len(shifts),:3] = data[:-1] * shifts[i] + data[1:] * ( 1 - shifts[i])            
            newData[i:-1:len(shifts),3] = colors[:-1]
        newData[-1,:3] = data[-1]
        newData[-1,3] = colors[-1]
                    
        towrite = tempfile.NamedTemporaryFile()
        towrite.write("%d\n\n"%(len(newData)))  #number of atoms and a blank line after is a requirement of rasmol
            
        for i in newData:                     
            towrite.write("CA\t%lf\t%lf\t%lf\t%d\n" % tuple(i)) 
        towrite.flush()
        "TODO: rewrite using subprocess.popen" 
        if os.name == "posix":  #if linux 
            os.system("rasmol -xyz %s -script %s" % (towrite.name, rascript.name))
        else:     #if windows 
            os.system("C:/RasWin/raswin.exe -xyz %s -script %s" % (towrite.name, rascript.name))         #For windows you might need to change the place where your rasmol file is
        rascript.close()
        towrite.close() 
        

class SimulationWithCrosslinks(Simulation):
    
    def addRandomCrosslinks(self,num=10):
        data = numpy.array(self.getData())
        for _ in xrange(num):
            a = numpy.random.randint(0,self.N)            
            cur = data[a]
            dist = numpy.sqrt(numpy.sum((data - cur[None,:])**2,axis = 1))
            attempt = dist < 1.9
             
            attempt[a] = False
            attempt[a-1] = False
            attempt[a-2] = False
            attempt[a+1] = False
            attempt[a+2] = False
            if attempt.sum() <= 0: 
                print "no atoms"
                continue
            nums = numpy.nonzero(attempt)[0]
            
            b = nums[numpy.random.randint(0,len(nums))]

            if self.verbose == True: print "random bond added",  self.dist(a,b)
            print a,b
            self.addBond(a,b,0.5)
            
            
    def addPerChainRandomCrosslinks(self,num=10,chain = 0):
        data = numpy.array(self.getData())
        for _ in xrange(num):
            while True: 
                mychain = self.chains[chain]
                a = numpy.random.randint(mychain[0],mychain[1])            
                cur = data[a]
                dist = numpy.sqrt(numpy.sum((data - cur[None,:])**2,axis = 1))
                attempt = dist < 1.9
                 
                attempt[a] = False
                attempt[a-1] = False
                attempt[a-2] = False
                attempt[a+1] = False
                attempt[a+2] = False
                attempt[:mychain[0]] = False
                attempt[mychain[1]:] = False

                if attempt.sum() <= 0: 
                    print "no atoms"
                    continue
                nums = numpy.nonzero(attempt)[0]
                
                b = nums[numpy.random.randint(0,len(nums))]
    
                if self.verbose == True: print "random bond added",  self.dist(a,b)
                print a,b
                self.addBond(a,b,0.5)
                break
    
    def addConsecutiveRandomBonds(self,bondlength,bondRange,distance,smeer = 0.2):
        shift = int(bondlength * smeer)
        begin = numpy.random.randint(shift)
        while True:
            b1 = begin
            b2 = begin + bondlength
            if b2 > self.N - 3: break
            self.addBond(b1,b2,bondRange,distance)
            begin = begin + bondlength + numpy.random.randint(shift) + shift/2
            if self.verbose == True: print "bond added between %d and %d" % (b1,b2)
            
    def addDoubleRandomLengthBonds(self,bondlength,bondRange,distance):
        begin = 4
        started = True
        past = 0
        while True:
            past
            b1 = begin
            b2 = begin + numpy.random.randint(0.5*bondlength,1.7*bondlength)
            if b2 > self.N -4: break
            self.addBond(b1,b2,bondRange,distance)
            if self.verbose == True: print "bond added between %d and %d" % (b1,b2)
            if started == False:
                self.addBond(past,b2,bondRange,distance)
                if self.verbose == True: print "bond added between %d and %d" % (past,b2)
                past = b1
            started = False
            begin = b2
    
    

class ExperimentalSimulation(Simulation):
    "contains some experimental features"

    def quickLoad(self,data,mode = "chain",Nchains = 1,trunc = None,confinementDensity = "NoConfinement"):        
        "quickly loads a set of repulsive chains, possibly adds spherical confinement"
        self.setup()
        self.load(data)
        self.setLayout(mode,Nchains)
        self.addHarmonicPolymerBonds()
        self.addSimpleRepulsiveForce(trunc = trunc)
        if type(confinementDensity) != str: self.addSphericalConfinement( density = confinementDensity)

    def createWalls(self,left = None, right = None,k=0.5):
        "creates walls at x = left, x = right, x direction only"
        if left == None: 
            left = self.data[0][0] + 1. * nm
        else:
            left = 1. * nm * left
        if right == None: 
            right = self.data[-1][0] - 1. * nm
        else:
            right = 1. * nm * right 
        
        if self.verbose == True:
            print "left wall created at ", left / (1. * nm) 
            print "right wall created at ", right / ( 1. * nm)
        
        extforce2 = openmm.CustomExternalForce(" WALLk * (sqrt((x - WALLright) * (x-WALLright) + WALLa * WALLa ) - WALLa) * step(x-WALLright) + WALLk * (sqrt((x - WALLleft) * (x-WALLleft) + WALLa * WALLa ) - WALLa) * step(WALLleft - x) ")
        extforce2.addGlobalParameter("WALLk",k*self.kT/nm)
        extforce2.addGlobalParameter("WALLleft",left)
        extforce2.addGlobalParameter("WALLright",right)
        extforce2.addGlobalParameter("WALLa",1 * nm)
        for i in xrange(self.N):
            extforce2.addParticle(i,[])
        self.forceDict["WALL Force"] = extforce2

    
    def addSphericalWell(self,r = 10,depth = 1):
        "pushes particles towards a boundary of a cylindrical well to create uniform well coverage"       

        extforce4 = openmm.CustomExternalForce("WELLdepth * (((sin((WELLr * 3.141592 * 0.5) / WELLwidth)) ^ 10)  -1) * step(-WELLr + WELLwidth)  ;WELLr = sqrt(x^2 + y^2 + z^2 + WELLtt^2)")
        self.forceDict["Well attraction"] = extforce4
        
        
        #adding all the particles on which force acts
        for i in xrange(self.N): 
            if self.domains[i] > 0.5: extforce4.addParticle(i,[])
        if r == None:
            try: r = self.sphericalConfinementRadius * 0.5
            except: exit("No spherical confinement radius defined yet. Apply spherical confinement first!")
        if self.verbose == True: print "Well attraction added with r = %d"% r 
         
        #assigning parameters of the force 
        extforce4.addGlobalParameter("WELLwidth",r*nm);
        extforce4.addGlobalParameter("WELLdepth",depth * self.kT);        
        extforce4.addGlobalParameter("WELLtt",0.01*nm);


    def add_cubic_grid(self,step,percent = 0.5,mass = 0.1):
        """In PBC this function adds a cubic grid of atoms (solvent),
            - separated by step units of length
            - with the probability percent
            -with mass = mass
            -in a box defined by PBC
            """

        print "Be sure that changing to MULT=1 didn't affect it... "
        raw_input("type any key to continue...") 
        size = self.SolventGridSize  #using box that was written in setup()
        newdata = []                 #array of solvent molecules
        count = [int(i / step) for i in size]  #grid size
        myarray = numpy.zeros(tuple(count),int) + 1  #grid itself
        for particle in self.getData():              #exclude monomers 
            p = [int(i) for i in particle/step]
            p = numpy.array(p)
            lb = p - 1
            lb[lb<0] = 0 
            
            myarray[lb[0] : p[0] + 2,lb[1]  : p[1] + 2,lb[2]: p[2] + 2] = 0 
        myarray[numpy.random.random(tuple(count)) > percent] = 0
        print "Using solvent with %d particles," % numpy.sum(myarray) 
        for i in xrange(count[0]):
            for j in xrange(count[1]):
                for k in xrange(count[2]):
                    if myarray[i,j,k] == 1:
                        newdata.append((step * i,step * j,step * k))
        newdata = numpy.array(newdata)
        for i in newdata: 
            self.system.addParticle(self.mass * mass)
        pastN = len(self.data)
        fulldata = numpy.concatenate([self.getData(),newdata])
        self.setData(fulldata)
        self.domains = numpy.zeros(self.N,int)
        self.domains[:pastN] = 1
        self.domains[pastN:] = 0
                    
        
class YeastSimulation(Simulation):
    """
    This is a multi-line definition. 
    I use it to test some bugs. 
    
    This class is some crazy Geoff's work.
    Don't even look at it.  
    """

    
    def addNucleolus(self, k = 1, r =  None):
        "method"
        if r==None: r =  self.sphericalConfinementRadius
        
        
        extforce3 = openmm.CustomExternalForce("step(r-NUCaa) * NUCkb * (sqrt((r-NUCaa)*(r-NUCaa) + NUCt*NUCt) - NUCt) ;r = sqrt(x^2 + y^2 + (z + NUCoffset )^2 + NUCtt^2)")
        
        self.forceDict["NucleolusConfinement"] = extforce3
        #adding all the particles on which force acts
        if self.verbose == True: print "NUCleolus confinement from radius = %lf" % r 
        #assigning parameters of the force 
        extforce3.addGlobalParameter("NUCkb",k*self.kT/nm);
        extforce3.addGlobalParameter("NUCaa", (r - 1./k)*nm * 1.75 );
        extforce3.addGlobalParameter("NUCoffset", (r - 1./k)*nm * 1.1 );
        extforce3.addGlobalParameter("NUCt",(1./k)*nm/10.);
        extforce3.addGlobalParameter("NUCtt",0.01*nm);
        for i in xrange(self.N): extforce3.addParticle(i,[])


    def addLaminaAttraction(self,width = 1,depth = 1, r = None, particles = None):
        "method"
        extforce3 = openmm.CustomExternalForce("-1 * step(LAMr-LAMaa + LAMwidth) * step(LAMaa + LAMwidth - LAMr) * LAMdepth * abs( (LAMr-LAMaa + LAMwidth) * (LAMaa + LAMwidth - LAMr)) / (LAMwidth * LAMwidth)  ;LAMr = sqrt(x^2 + y^2 + z^2 + LAMtt^2)")
        self.forceDict["Lamina attraction"] = extforce3
        
        #re-defines lamina attraction based on particle index instead of domains.
        
        #adding all the particles on which force acts
        if particles == None:
            for i in xrange(self.N): 
                extforce3.addParticle(i,[])
                if self.verbose == True: print "particle %d laminated! " % i 

        else:
            for i in particles:
                extforce3.addParticle(i,[])
                if self.verbose == True: print "particle %d laminated! " % i 
        
        if r == None:
            try: r = self.sphericalConfinementRadius
            except: exit("No spherical confinement radius defined yet. Apply spherical confinement first!")
        
        if self.verbose == True: print "Lamina attraction added with r = %d"% r 
         
        #assigning parameters of the force 
        extforce3.addGlobalParameter("LAMaa", r * nm);
        extforce3.addGlobalParameter("LAMwidth",width*nm);
        extforce3.addGlobalParameter("LAMdepth",depth * self.kT);        
        extforce3.addGlobalParameter("LAMtt",0.01*nm);

                      
                
#spiral = create_spiral(3,5.4,6000)        
#rw = createRW(60000)


