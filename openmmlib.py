import numpy
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
    def __init__(self,start_mode = "cold",timestep=25,thermostat=0.01,verbose = False,velocityReinitialize = True, name = "sim"):
        self.name = name        
        self.timestep = timestep * fs
        self.collisionRate = thermostat * ps        
        self.verbose = verbose
        self.velocityReinialize = velocityReinitialize
        self.loaded = False  #check if the data is loaded
        self.forcesApplied = False
        self.folder = "."
        self.nm = nm 
        self.metadata = {}
    
        
    def add_cubic_grid(self,step,percent = 0.5,mass = 0.1):
        print "Be sure that changing to MULT=1 didn't affect it... "
        raw_input("type any key to continue...") 
        """In PBC this function adds a cubic grid of atoms (solvent),
            - separated by step units of length
            - with the probability percent
            -with mass = mass
            -in a box defined by PBC
            """
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
        
                        
        
    def quickLoad(self,data,mode = "chain",Nchains = 1,trunk = None):
        "quickly loads one chain"
        self.setup()
        self.load(data)
        self.setLayout(mode,Nchains)
        self.addHarmonicPolymerBonds()
        self.addSimpleRepulsiveForce(trunk = trunk)
        
    def setup(self,platform="OpenCL", PBC = False,density = False,PBCbox = None):           
        "sets up the system, autodetects the size of PBC box unless passed explicitely"
        self.step = 0
        if PBC == True: 
            self.metadata["PBC"] = True             
        try: myplatform = self.platformName
        except: 
            myplatform = platform
            self.platformName = platform 
        
        try: self.temperature
        except: self.temperature = 300 * units.kelvin  #we can redefine the temperature sometimes... if we need it. 
        try: self.collisionRate
        except: self.collisionRate = 0.01 * ps
        self.kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA # Boltzmann constant
        self.kT = self.kB * self.temperature # thermal energy
        self.beta = 1.0 / self.kT # inverse temperature
        self.mass = 100.0 * units.amu
        self.bondsForException = []
        self.mm = openmm
        self.conlen = 1. * self.nm
        self.system = self.mm.System()
        self.PBC = PBC
        if self.PBC == True:
            if PBCbox == None:            
                data = self.getData()
                data -= numpy.min(data,axis = 0)
                self.setData(data)
                datasize = 1.1 * (2+(numpy.max(self.getData(),axis = 0) - numpy.min(self.getData(),axis = 0))) #size of the system plus some overhead                          
                self.SolventGridSize = (datasize / 1.1)  - 2
                realBoxSize = datasize
                realVolume = realBoxSize[0] * realBoxSize[1] * realBoxSize[2]
                mydensity = self.N / realVolume
                print "density is ", mydensity
            else:
                PBCbox = numpy.array(PBCbox) 
                                
                datasize = PBCbox 
            self.metadata["PBCbox"] = PBCbox
            self.system.setDefaultPeriodicBoxVectors([datasize[0],0.,0.],[0.,datasize[1],0.],[0.,0.,datasize[2]])
            self.BoxSizeReal = datasize 
            print "system size is %lfx%lfx%lf nm, solvent grid size is   %lfx%lfx%lf in dimensionless units" % tuple(list(datasize) + list(self.SolventGridSize))
            time.sleep(5)
             
            
            
        
        try: self.GPU 
        except: self.GPU = "0"  #setting default GPU  

        if myplatform == "OpenCL":
            platform = openmm.Platform.getPlatformByName('OpenCL')
            platform.setPropertyDefaultValue('OpenCLDeviceIndex', self.GPU)
        elif myplatform == "Reference":
            platform =  openmm.Platform.getPlatformByName('Reference')
        elif platform == "CUDA":
            platform = openmm.Platform.getPlatformByName('Cuda')
            platform.setPropertyDefaultValue('CudaDevice', self.GPU)
        else:
            raw_input("\n!!!!!!!!!!unknown platform!!!!!!!!!!!!!!!\n")
            exit()
        self.platform = platform
            
        "clearing ForceDict to recreate it"
        self.forceDict = {}
        try:
            for _ in xrange(self.N):
                self.system.addParticle(self.mass)
            print self.N, "particles loaded"
        except: pass
                                 
        self.integrator = openmm.LangevinIntegrator(self.temperature, self.collisionRate, self.timestep)
        
    def use_platform(self,platform):
        "sets the platform"
        self.platformName = platform
        self.setup()

    def exit(self,line):
        print line
        print "--------------> Bye <---------------"
        exit()
        
            
    
    def load(self,filename,center = False):        
        "loads file from  an Nx3 or 3xN or filename"
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
                print "Working with joblib file"
                mydict = dict(joblib.load(filename))
                data = mydict.pop("data")
                self.oldMetadata = data 
        else:
            data = filename
                    
        data = numpy.asarray(data,float)
        if len(data) == 3: 
            data = numpy.transpose(data)
        if len(data[0]) != 3: 
            self.exit("strange data file")        
        if numpy.isnan(data).any(): self.exit("\n!!!!!!!!!!file contains NANS!!!!!!!!!\n")
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

    
    def RG(self):
        "returns radius of gyration"
        return numpy.sqrt(numpy.sum(numpy.var(numpy.array(self.data/self.nm) ,0)))
    
    def RMAX(self):
        "returns distance to the furthest particle"
        return numpy.max(numpy.sqrt(numpy.sum((numpy.array(self.data/self.nm) )**2,1)))
        
    def setLayout(self,mode = "chain",chains = None,Nchains = 1):
        "sets layout of chains for chains or rings. "
        "chains: [0 100],[100 200],... etc"
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

    def useDomains(self,mode,chromosomes = None, domains = "domains.dat"):
        "in combination with LennardJonseForce this allows to create domains from Hi-C eigenvectors"
        self.domains = numpy.zeros(self.N,int)

        if mode == "file":
            domains = cPickle.load(open(domains))
            if len(domains) != self.N: exit("Wrong domain length")
            self.domains = domains
        if mode == "domains":
            if type(domains) == type("bla"): domains = cPickle.load(open(domains))
            if len(domains) != self.N: exit("Wrong domain length")
            self.domains = domains
            cPickle.dump(domains,open(self.folder + "/domains.dat",'wb'))            
                
    def initHarmonicBondForce(self):
        "Inits harmonic forse for polymer and non-polymer bonds"
        if "HarmonicBondForce" not in self.forceDict.keys(): 
            self.forceDict["HarmonicBondForce"] = self.mm.HarmonicBondForce()
        self.bondType = "Harmonic"        

          
    def addHarmonicPolymerBonds(self,wiggleDist = 0.05):
        "adds harmonic polymer bonds"      
        for i in self.chains:
            for j in xrange(i[0],i[1] - 1):
                self.addBond(j,j+1,wiggleDist,distance = 1, bondType = "Harmonic",verbose = False)

            if self.mode == "ring":
                self.addBond(i[0],i[1] -1,wiggleDist,distance = 1, bondType = "Harmonic")            
                if self.verbose == True: print "ring bond added", i[0],i[1]-1
                
        self.metadata["HarmonicPolymerBonds"] = {"wiggleDist":wiggleDist}

    
    def addGrosbergPolymerBonds(self,k = 30):
        self.initGrosbergBondForce()
        force = self.forceDict["GrosbergBondForce"]
        for i in self.chains:
            for j in xrange(i[0],i[1] - 1):
                force.addBond(j,j+1,[])               
            if self.mode == "ring":
                force.addBond(i[0],i[1] - 1,[])
                print "chain bond added", i[0],i[1]-1
        self.metadata["GorsbergPolymerForce"] = {"k":k}
        
                

    def addStifness(self,k = 40):
        myforce= self.mm.CustomAngleForce("k * (theta - 3.141592) * (theta - 3.141592)")
        self.forceDict["AngleForce"] = myforce 
        for i in self.chains:
            for j in xrange(i[0]+1,i[1] - 1):
                myforce.addAngle(j-1,j,j+1,[])
        myforce.addGlobalParameter("k",k)
        self.metadata["AngleForce"] = {"stifness":k}


    def addGrosbergStifness(self, k = 1.5):
        myforce= self.mm.CustomAngleForce("k * (1 - cos(theta - 3.141592))")
        self.forceDict["AngleForce"] = myforce 
        for i in self.chains:
            for j in xrange(i[0]+1,i[1] - 1):
                myforce.addAngle(j-1,j,j+1,[])
        myforce.addGlobalParameter("k",k * self.kT)
        self.metadata["GrosbergAngleForce"] = {"stifness":k}

        
    def addSimpleRepulsiveForce(self,cutoff = 1.7,trunk = None,rep = 0.26):
        "creates a force where each particle repells eath other particle"
        nbCutOffDist = self.conlen * cutoff    #repulsive part saturates quickly 
        if trunk == None: repul_energy = 'REPepsilon*(REPsigma/r)^12'
        else: repul_energy = '1/((1/REPcutoff) + (1/REPU + 0.0001 * REPcutoff));REPU=REPepsilon*(REPsigma/r)^12'
        
        epsilonRep = rep * units.kilocalorie_per_mole
        sigmaRep = 1.06 * self.conlen
        
        #creating  a custom repulsive force
        self.forceDict["Nonbonded"] = openmm.CustomNonbondedForce(repul_energy)
        repulforce = self.forceDict["Nonbonded"]
        if trunk!=None: repulforce.addGlobalParameter('REPcutoff',trunk * self.kT)
        repulforce.addGlobalParameter('REPepsilon',epsilonRep)
        repulforce.addGlobalParameter('REPsigma',sigmaRep) 
        for n in range(self.N):
            repulforce.addParticle(())        
        # Set non-periodic boundary conditions with cutoff.
        if self.PBC: 
            repulforce.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
            print "Using periodic boundary conditions!!!!!!!!11"
            
        else: repulforce.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic)
        repulforce.setCutoffDistance(nbCutOffDist)
        self.metadata["SimgleRepulsiveForce"] = {"cutoff":cutoff,"trunk":trunk,"rep":rep}
        
    
    def addGrosbergRepulsiveForce(self,trunk=None):
        nbCutOffDist = self.conlen * 2.**(1./6.)
        if trunk == None:     
            repul_energy = "4 * REPe * ((REPs/r)^12 - (REPs/r)^6) + REPe"
        else: 
            repul_energy = "1 / (1 / REPcut + 1 / (REPU0 + REPa * REPcut) ) - REPcut * (REPa/(REPa+1)) ;REPU0 = 4 * REPe * ((REPs/r2)^12 - (REPs/r2)^6) + REPe;r2 = (r^10. + (0.3 * REPs)^10.)^0.1"
            
        self.forceDict["Nonbonded"] = openmm.CustomNonbondedForce(repul_energy)
        repulforce = self.forceDict["Nonbonded"]
        repulforce.addGlobalParameter('REPe',self.kT)
        repulforce.addGlobalParameter('REPs',self.conlen)
        if trunk != None: 
            repulforce.addGlobalParameter('REPcut',self.kT * trunk)        
            repulforce.addGlobalParameter('REPa',0.001)
        for n in range(self.N):
            repulforce.addParticle(())        
        # Set non-periodic boundary conditions with cutoff.
        if self.PBC: 
            repulforce.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
            print "Using periodic boundary conditions!!!!!!!!11"            
        else: repulforce.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic)
        repulforce.setCutoffDistance(nbCutOffDist)
        
        
        
    def addMutualException(self, particles):
        for i in particles: #xrange(len(particles)):
            for j in particles:#xrange(len(particles)):
                if j> i:        
                    self.bondsForException.append((i,j))



        
    def dist(self,i,j):
        "returns distance between particles i and j"                
        dif = self.data[i] - self.data[j]
        dif = dif / (self.nm)
        dif = numpy.sqrt(sum(dif**2))        
        return dif
        
    def addLennardJonesForce(self,cutoff = 2.5, 
                            domains = False,
                            epsilonRep = 0.24,epsilonAttr = 0.27,   #parameters for LJ force 
                            blindFraction = -1,
                            sigmaRep = None, sigmaAttr = None):
        """adds repulsive-attractive force. epsilonRep < 0.24 is mainly repulsive, >0.24 becomes attractive.   !!!check in now!!! 
         You can specify parameters differently per domains. 
         If you want to skip random set of particles, modify "blind percent" - it'll skip particles with probability blindPercent 
         """
        if blindFraction > 0.99: self.exit ("why do you need this force without particles??? set blindFraction between 0 and 1") 
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
        
        
    def addInteraction(self,i,j,epsilon,sigma = None, length = 3):
        """adds attractive short-range interaction of strength epsilon between particles i,j
        requires LennardJones Force
        """
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

            
        
        
        
            
            
    def addSphericalConfinement(self,r="density",  k = 1., density = None):
        "constrain particles to be within a sphere. With no parameters creates sphere with density 0.25. "
        
        extforce2 = openmm.CustomExternalForce("step(r-SPHaa) * SPHkb * (sqrt((r-SPHaa)*(r-SPHaa) + SPHt*SPHt) - SPHt) ;r = sqrt(x^2 + y^2 + z^2 + SPHtt^2)")
        self.forceDict["SphericalConfinement"] = extforce2
        
        #adding all the particles on which force acts
        for i in xrange(self.N): extforce2.addParticle(i,[])
        if r == "density":
            r = 13.  * ((self.N/4000.) ** 0.333333333333333)
        if density != None:
            r = (3 * self.N / (4 * 3.141592 * density)) ** (1/3.)
            k = 5
             
        self.sphericalConfinementRadius = r
        if self.verbose == True: print "SPherical confinement with radius = %lf" % r 
        #assigning parameters of the force 
        extforce2.addGlobalParameter("SPHkb",k*self.kT/self.nm);
        extforce2.addGlobalParameter("SPHaa",(r - 1./k)*self.nm);
        extforce2.addGlobalParameter("SPHt",(1./k)*self.nm/10.);
        extforce2.addGlobalParameter("SPHtt",0.01*self.nm);


    ### MARKED FOR MIGRATION TO YEAST CLASS ####
    def addNucleolus(self, k = 1, r =  None):
        if r==None: r =  self.sphericalConfinementRadius
        
        
        extforce3 = openmm.CustomExternalForce("step(r-NUCaa) * NUCkb * (sqrt((r-NUCaa)*(r-NUCaa) + NUCt*NUCt) - NUCt) ;r = sqrt(x^2 + y^2 + (z + NUCoffset )^2 + NUCtt^2)")
        
        self.forceDict["NucleolusConfinement"] = extforce3
        #adding all the particles on which force acts
        if self.verbose == True: print "NUCleolus confinement from radius = %lf" % r 
        #assigning parameters of the force 
        extforce3.addGlobalParameter("NUCkb",k*self.kT/self.nm);
        extforce3.addGlobalParameter("NUCaa", (r - 1./k)*self.nm * 1.75 );
        extforce3.addGlobalParameter("NUCoffset", (r - 1./k)*self.nm * 1.1 );
        extforce3.addGlobalParameter("NUCt",(1./k)*self.nm/10.);
        extforce3.addGlobalParameter("NUCtt",0.01*self.nm);
        for i in xrange(self.N): extforce3.addParticle(i,[])


    def addLaminaAttraction(self,width = 1,depth = 1, r = None):
        extforce3 = openmm.CustomExternalForce("step(LAMr-LAMaa + LAMwidth) * step(LAMaa + LAMwidth - LAMr) * LAMdepth * (LAMr-LAMaa + LAMwidth) * (LAMaa + LAMwidth - LAMr) / (LAMwidth * LAMwidth)  ;LAMr = sqrt(x^2 + y^2 + z^2 + LAMtt^2)")
        self.forceDict["Lamina attraction"] = extforce3

        #adding all the particles on which force acts
        for i in xrange(self.N): 
            if self.domains[i] > 0.5: extforce3.addParticle(i,[])
        if r == None:
            try: r = self.sphericalConfinementRadius
            except: exit("No spherical confinement radius defined yet. Apply spherical confinement first!")
        if self.verbose == True: print "Lamina attraction added with r = %d"% r 
         
        #assigning parameters of the force 
        extforce3.addGlobalParameter("LAMaa", r * self.nm);
        extforce3.addGlobalParameter("LAMwidth",width*self.nm);
        extforce3.addGlobalParameter("LAMdepth",depth * self.kT);        
        extforce3.addGlobalParameter("LAMtt",0.01*self.nm);


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
        extforce4.addGlobalParameter("WELLwidth",r*self.nm);
        extforce4.addGlobalParameter("WELLdepth",depth * self.kT);        
        extforce4.addGlobalParameter("WELLtt",0.01*self.nm);
        
    
    def tetherParticles(self,particles, k = 30):
        "tethers particles in the 'particles' array. Increase k to tether them stronger, but watch the system!" 
        extforce2 = openmm.CustomExternalForce(" TETHkb * ((x - TETHx0)^2 + (y - TETHy0)^2 + (z - TETHz0)^2)")
        self.forceDict["Tethering Force"] = extforce2

        #assigning parameters of the force 
        extforce2.addGlobalParameter("TETHkb",k*self.kT/self.nm);
        extforce2.addPerParticleParameter("TETHx0")
        extforce2.addPerParticleParameter("TETHy0")
        extforce2.addPerParticleParameter("TETHz0")
        for i in particles:#adding all the particles on which force acts
            coordinates = self.data[i]
            extforce2.addParticle(i,list(coordinates))
            if self.verbose == True: print "particle %d tethered! " % i 
        
        
        
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
        extforce2.addGlobalParameter("WALLk",k*self.kT/self.nm)
        extforce2.addGlobalParameter("WALLleft",left)
        extforce2.addGlobalParameter("WALLright",right)
        extforce2.addGlobalParameter("WALLa",1 * nm)
        for i in xrange(self.N):
            extforce2.addParticle(i,[])
        self.forceDict["WALL Force"] = extforce2
        
            
    def addCylindricalConfinement(self,r="density",bottom = True,k=0.1,weired = False):
        "as it says" 
        
        if bottom == True: extforce2 = openmm.CustomExternalForce("step(r-CYLaa) * CYLkb * (sqrt((r-CYLaa)*(r-CYLaa) + CYLt*CYLt) - CYLt) + step(-z) * CYLkb * (sqrt(z^2 + CYLt^2) - CYLt) ;r = sqrt(x^2 + y^2 + CYLtt^2)")
        else: extforce2 = openmm.CustomExternalForce("step(r-CYLaa) * CYLkb * (sqrt((r-CYLaa)*(r-CYLaa) + CYLt*CYLt) - CYLt) ;r = sqrt(x^2 + y^2 + CYLtt^2)") 
        if weired == True: extforce2 = openmm.CustomExternalForce(" 0.6 * CYLkt * CYLaa*CYLaa / (CYLaa * CYLaa + r * r) + step(r-CYLaa) * CYLkb * (sqrt((r-CYLaa)*(r-CYLaa)  + CYLt*CYLt) - CYLt) + step(-z) * CYLkb * (sqrt(z^2 + CYLt^2) - CYLt) ;r = sqrt(x^2 + y^2 + CYLtt^2)")
        self.forceDict["CylindricalConfinement"] = extforce2

        #adding all the particles on which force acts
        for i in xrange(self.N): extforce2.addParticle(i,[])
        if r == "density":
            r = 13. * 3 * ((self.N/4000.) ** 0.333333333333333)
         
        #assigning parameters of the force 
        extforce2.addGlobalParameter("CYLkb",k*self.kT/self.nm);
        extforce2.addGlobalParameter("CYLkt",self.kT);
        extforce2.addGlobalParameter("CYLweired",self.nm)
        extforce2.addGlobalParameter("CYLaa",(r - 1./k)*self.nm);
        extforce2.addGlobalParameter("CYLt", (1./(10*k))*self.nm);
        extforce2.addGlobalParameter("CYLtt",0.01*self.nm);
        
    def getData(self):
        "use this to return data, not something else"
        return self.data / self.nm
    def getScaledData(self):        
        "returns back data rescaled to the simulations box"
        alldata = self.getData()
        boxsize = numpy.array(self.BoxSizeReal)
        mults = numpy.array(alldata / boxsize[None,:],int)
        return alldata - mults * boxsize[None,:]
    def setData(self,data):
        """use this, not self.data, to set data.
         internal realization may change later
        """
        data = numpy.asarray(data,dtype = "float")
        self.data = units.Quantity(data,self.nm)
        self.N = len(self.data)    
    def addGravity(self,k = 0.1,cutoff = None):
        """adds force pulling downwards in z direction
        When using cutoff, acts only when z>cutoff
        """
        
        if cutoff == None:
            extforce3 = openmm.CustomExternalForce("kG * z")
            extforce3.addGlobalParameter("kG",k * self.kT / (self.nm))
            for i in xrange(self.N): extforce3.addParticle(i,[])
            self.forceDict["Gravity"] = extforce3
        else:
            extforce3 = openmm.CustomExternalForce("kG * (z - cutoffG) * step(z - cutoffG)")
            extforce3.addGlobalParameter("kG",k * self.kT / (self.nm))
            extforce3.addGlobalParameter("cutoffG", cutoff * self.nm)
            for i in xrange(self.N): extforce3.addParticle(i,[])
            self.forceDict["Gravity"] = extforce3
     


    def initGrosbergBondForce(self):
        "inits Grosberg FENE bond force"
        if "GrosbergBondForce" not in self.forceDict.keys():
            force = "- 0.5 * GROSk * GROSr0 * GROSr0 * log(1-(r/GROSr0)* (r / GROSr0))"
            bondforce3  = openmm.CustomBondForce(force)
            
            bondforce3.addGlobalParameter("GROSk",30 * self.kT / (self.conlen * self.conlen))
            bondforce3.addGlobalParameter("GROSr0",self.conlen * 1.5)
            self.forceDict["GrosbergBondForce"] = bondforce3

        
        
    def addBond(self,i,j,bondWiggleDistance,distance = None, bondType = None,verbose = None):
        
        #print (self.data[i,:]- self.data[j,:]), ((self.data[i,:]- self.data[j,:])**2).sum()
        #if  ((self.data[i,:]- self.data[j,:])**2).sum()**.5 > 1.4: self.exit ('asdfa') # sanity check on distance between molecules
        if verbose == None:
            verbose = self.verbose 
        
        bondSize = float(bondWiggleDistance)
        if distance == None: distance = self.conlen / nm
        else:  distance = self.conlen * distance / nm
        
        distance  = float(distance)        
        if bondType == None: 
            bondType = self.bondType             
        if bondType == "Harmonic":
            self.initHarmonicBondForce()
            kbond = ( 2 * self.kT / (bondSize * self.conlen) ** 2 ) / (units.kilojoule_per_mole / nm ** 2 )
            self.forceDict["HarmonicBondForce"].addBond(int(i),int(j),float(distance),float(kbond))
            if verbose == True: print "Harmonic bond added between %d,%d, params %lf %lf" % (i,j,float(distance),float(kbond))
        elif bondType == "Grosberg":
            self.initGrosbergForce()
            self.forceDict["GrosbergForce"].addBond(i,j,[])
        else: self.exit("Bond type not known")
        
    def addRandomCrosslinks(self,num=10):
        data = numpy.array(self.getData())
        for i in xrange(num):
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
        for i in xrange(num):
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
            
             
            
         
            
            
                                  
    def applyForces(self):
        ##making comprehensive list of exceptions
        exc = self.bondsForException
        
        print "exeptionlength", len(exc)
        
        if len(exc) > 0:
            exc = numpy.array(exc) 
            exc = numpy.sort(exc,axis = 1) 
            exc = [tuple(i) for i in exc]
            exc = list(set(exc)) 

        for i in self.forceDict.keys():

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
    def initVelocities(self,mult = 1,massmult = 'auto'):
        "sets velocities."
        if massmult == "auto":
            try:
                massmult = self.massmult
            except:
                massmult = 1 
        try: self.context
        except: print "No context, cannot set velocs. Initialize context before that"
        
        sigma = units.Quantity(numpy.zeros(1,float), units.sqrt(self.kT.unit/units.amu))
        mass = self.system.getParticleMass(0)
        sigma[0] = units.sqrt(self.kT/(mass * massmult))   #calculating mean velocity
        velocs = mult * units.Quantity(numpy.random.normal(size=(self.N,3)) * sigma,sigma.unit)
        
        #set initial configuration in the OpenMM
        self.context.setVelocities(velocs)
    def initPositions(self):
        "sets positions"
        self.context.setPositions(self.data)
        state = self.context.getState(getPositions=True,getEnergy = True)  #get state of a system: positions, energies
        eP = state.getPotentialEnergy()/self.N/self.kT
        
        print "Positions loaded, potential energy is %lf" % eP
    def reinitialize(self,mult = 1):
        "call this every time you change some parameters of the system between blocks"
        self.context.reinitialize()
        self.initPositions()
        self.initVelocities(mult)
    def saveFolder(self,folder):
        "sets the folder where to save data"
        if folder[-1] == "/": folder = folder[:-1]
        if os.path.exists(folder) == False:
            os.mkdir(folder)
        self.folder = folder
    def energy_minimization(self,steps = 1000,twoStage = False,thermostat = None):
        "run this first, to minimize energy."
        print "Performing energy minimization"
        if self.forcesApplied == False:
            self.applyForces()
            self.forcesApplied = True
        def_step = self.integrator.getStepSize()
        self.integrator.setStepSize(def_step/15)
        def_fric = self.integrator.getFriction()
        if thermostat == None: self.integrator.setFriction(1.3)
        else: self.integrator.setFriction(thermostat)
        self.reinitialize()
        
        for i in xrange(10): 
            self.doBlock(steps = steps, increment = False)
            self.initVelocities()
            
        if twoStage == True:
            self.integrator.setFriction(0.1)
            self.reinitialize()
            for i in xrange(5):            
                self.doBlock(steps = steps, increment = False)
                self.initVelocities()

        self.integrator.setFriction(def_fric)
        self.integrator.setStepSize(def_step)
        self.reinitialize()
        print "Finished energy minimization"        
    def doBlock(self,steps = None,increment = True,num = None):
        """performs one block of simulations, doing steps timesteps, or steps_per_block if not specified. 
        "increment" controls where we increase the self.step or not
        """
        
        if self.forcesApplied == False: 
            self.applyForces()
            self.forcesApplied = True
        if increment == True: self.step += 1
        if steps == None: steps = self.steps_per_block
        
        N = self.N
        kT = self.kT
        
        for attempt in xrange(6):
            print "%s  block=%d" % (self.name,self.step),
            if num == None:
                num = steps/10 + 1
            a = time.time()
            for i in xrange(steps/num):
                self.integrator.step(num)  #integrate!
                #state = self.context.getState(getEnergy = True)
                #eK = state.getKineticEnergy()/self.N/self.kT
                #if numpy.isnan(eK): break
                #if eK > 50: break
                print ".",
                sys.stdout.flush()
            self.integrator.step(steps%num)
            state = self.context.getState(getPositions=True,getEnergy = True)  #get state of a system: positions, energies
            b = time.time()
            coords = state.getPositions(asNumpy=True)
            newcoords = coords / self.nm
            eK = state.getKineticEnergy()/self.N/self.kT                          #calculate energies in KT/particle
            eP = state.getPotentialEnergy()/self.N/self.kT
            if self.velocityReinialize == True:
                if eK > 2.4:
                    self.initVelocities()   
            print " Coord[1]= [%.1lf %.1lf %.1lf] " % tuple(newcoords[0]),
            
            
            
            if (numpy.isnan(newcoords).any()) or (eK > 20) or (numpy.isnan(eK) ) or (numpy.isnan(eP)):
                self.context.setPositions(self.data)
                self.initVelocities()
                print "trying one more time at step # %i" % self.step
            else:
                self.data = coords
                print " %.2lf kin, %.2lf pot, %.2lf tot," % (eK,eP,eK+eP),  " Rg= %.3lf" % self.RG(),
                print "SPS = %.0lf:"%(steps/(float(b-a)))                
                break
            if attempt in [3,4]:
                self.energy_minimization(10)
            if attempt == 5:                                
                self.exit("exceeded number of attmpts")
                
    def addConsecutiveRandomBonds(self,bondlength,range,distance,smeer = 0.2):
        shift = int(bondlength * smeer)
        begin = numpy.random.randint(shift)
        while True:
            b1 = begin
            b2 = begin + bondlength
            if b2 > self.N - 3: break
            self.addBond(b1,b2,range,distance)
            begin = begin + bondlength + numpy.random.randint(shift) + shift/2
            if self.verbose == True: print "bond added between %d and %d" % (b1,b2)
            
    def addDoubleRandomLengthBonds(self,bondlength,range,distance):
        begin = 4
        started = True
        past = 0
        while True:
            past
            b1 = begin
            b2 = begin + numpy.random.randint(0.5*bondlength,1.7*bondlength)
            if b2 > self.N -4: break
            self.addBond(b1,b2,range,distance)
            if self.verbose == True: print "bond added between %d and %d" % (b1,b2)
            if started == False:
                self.addBond(past,b2,range,distance)
                if self.verbose == True: print "bond added between %d and %d" % (past,b2)
                past = b1
            started = False
            begin = b2
    def show(self,coloring = "chain",chain = "all" ):
        """shows system in rasmol by drawing spheres
        draws 4 spheres in between any two points (5 * N spheres total)"""
              
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
  
        if os.name == "posix":  #if linux 
            os.system("rasmol -xyz %s -script %s" % (towrite.name, rascript.name))
        else:     #if windows 
            os.system("C:/RasWin/raswin.exe -xyz %s -script %s" % (towrite.name, rascript.name))         #For windows you might need to change the place where your rasmol file is
        rascript.close()
        towrite.close() 
        
    
    
    def save(self,filename = None):
        "saves conformation plus some metadata. Metadata is not interpreted by this library, and is for your reference"
        self.metadata["data"] = self.getData()
        self.metadata["timestep"] = self.timestep / fs
        self.metadata["Collision rate"] = self.collisionRate / ps         
        
        if filename == None: 
            f = os.path.join(self.folder , "block%d.dat" % self.step)
        else:
            f = os.path.join(self.folder , filename)
        joblib.dump(self.metadata,filename = f,compress = 3)
                
        
class YeastSimulation(Simulation):
    
    def addNucleolus(self, k = 1, r =  None):
        if r==None: r =  self.sphericalConfinementRadius
        
        
        extforce3 = openmm.CustomExternalForce("step(r-NUCaa) * NUCkb * (sqrt((r-NUCaa)*(r-NUCaa) + NUCt*NUCt) - NUCt) ;r = sqrt(x^2 + y^2 + (z + NUCoffset )^2 + NUCtt^2)")
        
        self.forceDict["NucleolusConfinement"] = extforce3
        #adding all the particles on which force acts
        if self.verbose == True: print "NUCleolus confinement from radius = %lf" % r 
        #assigning parameters of the force 
        extforce3.addGlobalParameter("NUCkb",k*self.kT/self.nm);
        extforce3.addGlobalParameter("NUCaa", (r - 1./k)*self.nm * 1.75 );
        extforce3.addGlobalParameter("NUCoffset", (r - 1./k)*self.nm * 1.1 );
        extforce3.addGlobalParameter("NUCt",(1./k)*self.nm/10.);
        extforce3.addGlobalParameter("NUCtt",0.01*self.nm);
        for i in xrange(self.N): extforce3.addParticle(i,[])


    def addLaminaAttraction(self,width = 1,depth = 1, r = None, particles = None):
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
        extforce3.addGlobalParameter("LAMaa", r * self.nm);
        extforce3.addGlobalParameter("LAMwidth",width*self.nm);
        extforce3.addGlobalParameter("LAMdepth",depth * self.kT);        
        extforce3.addGlobalParameter("LAMtt",0.01*self.nm);

                      
                
#spiral = create_spiral(3,5.4,6000)        
#rw = createRW(60000)


