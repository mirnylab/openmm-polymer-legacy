import sys
#sys.path.append("/net/evolution/home/magus/OpenMM")


import os, cPickle
from math import sin,cos,sqrt
import numpy
from array import array
from numutils import createRW, coolAverage
from inout import removeAxes

def load(file,center=True):
    "loads xyz file "
    f = open(file, 'r')
    lines = f.readlines()
    N = int(lines[0])
    #print "loading half glubule!!!"
        
    datax = array('f', [0. for i in xrange(N)])
    datay = array('f', [0. for i in xrange(N)])
    dataz = array('f', [0. for i in xrange(N)])
    for i in xrange(1, N + 1):
        line = lines[i].split()
        datax[i - 1] = float(line[0])
        datay[i - 1] = float(line[1])
        dataz[i - 1] = float(line[2])
    if center==False:
        return numpy.array([numpy.array(datax), numpy.array(datay), numpy.array(dataz)])
    datax,datay,dataz = numpy.array(datax),numpy.array(datay),numpy.array(dataz)
    diffs = (datax[0:N-1]-datax[1:N])**2+(datay[0:N-1]-datay[1:N])**2+(dataz[0:N-1]-dataz[1:N])**2
    diffs = numpy.abs(1-numpy.sqrt(diffs))
    if numpy.max(diffs) > 0.6:
        print "error, %lf at %s" % (numpy.max(diffs),file)
        
    
    
    datax,datay,dataz = numpy.array(datax),numpy.array(datay),numpy.array(dataz)
    datax -= numpy.sum(datax)/len(datax)
    datay -= numpy.sum(datay)/len(datay)
    dataz -= numpy.sum(dataz)/len(dataz)
    #print numpy.sum(datax),numpy.sum(datay),numpy.sum(dataz)
    
    
    return numpy.transpose(numpy.array([datax,datay,dataz]))


def save(data,file,doint=False):
    if len(data) != 3: 
        data = numpy.transpose(data)
    if len(data) !=3: 
        print "wrong data"
        exit()
    f = open(file, 'w')
    f.write(str(len(data[0])) + "\n")
    for i in xrange(len(data[0])):
        if (doint == False):
            for j in xrange(len(data)):
                f.write(str(data[j][i]) + " ")
            f.write("\n")
        else:
            for j in xrange(len(data)):
                f.write(str(int(data[j][i])) + " ")
            f.write("\n")
            
    f.close()




def center(data):
    if len(data) == 3: data = numpy.transpose(data)
    data -= numpy.min(data,0)[None,:]
    data -= 0.5*numpy.max(data,0)[None,:]
    return data
    
def run_in_separate_process(func, *args, **kwds):

    pread, pwrite = os.pipe()
    pid = os.fork()
    if pid > 0:
        os.close(pwrite)
        with os.fdopen(pread, 'rb') as f:
            status, result = cPickle.load(f)
        os.waitpid(pid, 0)
        if status == 0:
            return result
        else:
            raise result
    else: 
        os.close(pread)
        try:
            result = func(*args, **kwds)
            status = 0
        except Exception, exc:
            result = exc
            status = 1
        with os.fdopen(pwrite, 'wb') as f:
            try:
                cPickle.dump((status,result), f, cPickle.HIGHEST_PROTOCOL)
            except cPickle.PicklingError, exc:
                cPickle.dump((2,exc), f, cPickle.HIGHEST_PROTOCOL)
        os._exit(0)

def create_spiral(r1,r2,N):
    Pi = 3.141592
    points = []
    finished = [False]
    def rad(phi):
        return phi/(2*Pi)
    def ang(rad):
        return 2*Pi*rad
    def coord(phi):
        r = rad(phi)
        return (r*sin(phi),r*cos(phi))
    def fullcoord(phi,z):
        c = coord(phi)
        return [c[0],c[1],z]
    def dist(phi1,phi2):
        c1 = coord(phi1)
        c2 = coord(phi2)
        d = sqrt((c1[1]-c2[1])**2+(c1[0]-c2[0])**2) 
        return d
    def nextphi(phi):
        
        phi1 = phi
        phi2 = phi + 0.7*Pi
        mid = phi2
        
        while abs(dist(phi,mid)-1 ) > 0.00001:            
            mid = (phi1 + phi2)/2.
            if dist(phi,mid) > 1 : phi2 = mid
            else: phi1 = mid
        return mid
    def prevphi(phi):
        
        phi1 = phi
        phi2 = phi - 0.7*Pi
        mid = phi2
        
        while abs(dist(phi,mid)-1 ) > 0.00001:            
            mid = (phi1 + phi2)/2.
            if dist(phi,mid) > 1 : phi2 = mid
            else: phi1 = mid
        return mid
    
    def add_point(point,points=points,finished=finished):        
        if (len(points) == N) or (finished[0] == True):
            points = numpy.array(points)            
            finished[0] = True
            print "finished!!!"
        else: 
            points.append(point)
            
    z = 0
    forward = True
    curphi = ang(r1)
    add_point(fullcoord(curphi, z))
    while True:
        if finished[0] == True:
            return numpy.transpose(points) 
        if forward == True:
            curphi = nextphi(curphi)
            add_point(fullcoord(curphi, z))
            if(rad(curphi) > r2):
                forward = False
                z+=1
                add_point(fullcoord(curphi, z))
        else: 
            curphi = prevphi(curphi)
            add_point(fullcoord(curphi, z))
            if(rad(curphi) < r1):
                forward = True
                z+=1
                add_point(fullcoord(curphi, z))

        
def resolve(data,chains = None,mode = "chain",steps = 3000):
    
    
    def all(mode = mode, data = data, chains = chains):
        
        from myopenmm import Simulation
        a = Simulation(timestep = 5, thermostat = 0.5,name = "resolve")
        a.setup()
        a.load(data)        
        if chains != None: a.setLayout(mode=mode, chains = chains)
        else: a.setLayout(mode = mode)
        a.addHarmonicPolymerBonds()
        a.addSimpleRepulsiveForce(trunk = 20,rep = 0.6)
        a.energy_minimization(steps, twoStage = True)
        return a.getData(), a.getLayout()

    return run_in_separate_process(all)


    
    
    

def create_sausage(N,ratio = 4,enlarge = 1):
        
    "creates sausage by dropping a spiral into the cylinder "
    D =( 5* N / (3. * ratio)) ** 0.33333 + 1
    D = D * enlarge
    print D
    print ratio * D
    spiral = create_spiral(3.1,3.3,N)
    def all():        
        import myopenmm        
        a = myopenmm.Simulation(timestep = 13, thermostat =0.1,name = "sausage")
        a.setup()
        a.load(spiral)
        a.setLayout(mode = "chain")
        a.addCylindricalConfinement(r=D/2.,bottom = True,k=2.5)
        a.addHarmonicPolymerBonds(0.06)
        a.addStifness(10)


        a.addSimpleRepulsiveForce()
        a.addGravity(0.3,cutoff = D * ratio )
        a.applyForces()
        a.steps_per_block = 4000
        a.energy_minimization(steps = 100)
        
        while a.getData()[:,2].max() >  D * ratio + 15:
            #print a.getData()[:,2].max(),  a.getData()[:,1].min(),             (ratio * D )
            a.doBlock()
            #a.initVelocities()
            #a.show()
        
        for i in xrange(5):
            a.doBlock()
        return a.getData(),a.getLayout()
    
    
        
            
    return   run_in_separate_process(all)


def createGrosberg():
    from myopenmm import Simulation
    a = Simulation(timestep = 25, thermostat = 0.004)
    a.load("/home/magus/evo/trajectories/globule_creation/6000_SAW/crumpled1.dat")
    a.saveFolder(".")
    a.setup(platform = "OpenCL")
    a.verbose = True
    a.setLayout()
    a.addSphericalConfinement(density = 0.8)
    a.addGrosbergPolymerBonds()
    #a.addSimpleRepulsiveForce()
    a.addGrosbergRepulsiveForce()
    a.addGrosbergStifness()
    
    #    a.addSphericalConfinement()
    #    a.useDomains(mode = "chromosomes",chromosomes = [2])
    #    print (a.domains == 1).sum()
    
    #    a.addLaminaAttraction(4,-.7,r=10)
    
    #a.addInteraction(1000,1200,0.3)
    #a.addInteraction(1800,2000,0.6)
    #a.addInteraction(2400,2600,0.9)
    #a.addInteraction(3000,3200,1.2)    
    #a.addInteraction(3600,3800,1.5)
    a.energy_minimization(steps = 1000,twoStage= True)
    for i in xrange(50000):
        a.doBlock(5000)
        a.save()
    
    
    
    import myopenmm
    res = resolve("/home/magus/evo/trajectories/globule_creation/4000_equilibrium/equilibrium1.dat", steps = 1200)
    a = myopenmm.Simulation(timestep = 23, thermostat =0.005,name = "sausage")
    a.setup()
    a.verbose = True
    a.load(res[0])
    a.setLayout(mode = "chain")
    a.addHarmonicPolymerBonds(0.07)
    a.addStifness(15)
    a.addSphericalConfinement(1, -0.05)
    a.saveFolder("/home/magus/workspace/testnucl2/pasta/run3")
    a.addSimpleRepulsiveForce()
    a.energy_minimization(600)
    for i in xrange(100):
        a.doBlock(3000)
        a.show() 
    exit()
    
    
    






def create_sausage_2(N,ratio = 4,enlarge = 1):
        
    "creates sausage by dropping a spiral into the cylinder "
    D =( 5* N / (3. * ratio)) ** 0.33333 + 1
    D = D * enlarge
    print D
    print ratio * D
    spiral = create_spiral(2.7,4.7,N)
    def all():        
        import myopenmm        
        a = myopenmm.Simulation(timestep = 23, thermostat =0.02,name = "sausage")
        a.setup()
        a.verbose = True
        a.load(spiral)
        a.setLayout(mode = "chain")
        a.saveFolder("/home/magus/workspace/testnucl2/pasta/run3")
        a.save("base")
        
        

        a.addCylindricalConfinement(r=D/2.,bottom = True,k=2.5,weired = True)
        a.addHarmonicPolymerBonds(0.07)
        
        
        #just pizza
        #a.addDoubleRandomLengthBonds(500,1.3,1.3)
        #a.addConsecutiveRandomBonds(150, 2,2,smeer = 0.6)
        #a.addConsecutiveRandomBonds(350, 2,2,smeer = 0.6)
        

        #real pizza
        #a.addDoubleRandomLengthBonds(2000,1.4,1.4)
        #a.addDoubleRandomLengthBonds(500,1.2,1.2)
        #a.addConsecutiveRandomBonds(35,1.5,1.5,smeer = 0.06)
        #a.addConsecutiveRandomBonds(220,1.5,1.5,smeer = 0.06)
        
              
        
        a.addStifness(5)


        a.addSimpleRepulsiveForce(1.35, 1.7)
        a.addGravity(0.06,cutoff = D * ratio )
        #a.addDoubleRandomLengthBonds(300,1.5,1.5)
        a.applyForces()
        a.steps_per_block = 3000
        a.energy_minimization(steps = 200)
        
        while a.getData()[:,2].max() >  D * ratio + 15:
            #print a.getData()[:,2].max(),  a.getData()[:,1].min(),             (ratio * D )
            a.doBlock()
            
            a.save()
            #if numpy.random.random() < 0.05: a.show()
            #a.initVelocities()
            #a.show()
        
        for i in xrange(1000):
            a.doBlock()
            a.save()
        return a.getData(),a.getLayout()
    
    
        
            
    return   run_in_separate_process(all)




def shrink(data,chains,radius = "density",mode = "chain"):
    
    def all(radius = radius, mode = mode, data = data, chains = chains):
        
        from myopenmm import Simulation
        a = Simulation(timestep = 10, thermostat = 0.2,name = "shrink")
        a.setup()
        a.load(data)
        if radius == "density":
            radius = 13. * ((a.N/4000.) ** 0.333333333333333) 
        
        a.setLayout(mode=mode, chains = chains)
        a.addHarmonicPolymerBonds()
        a.addSimpleRepulsiveForce()
        a.addSphericalConfinement(radius-3,k=0.2)
        
        while a.RMAX() > radius  + 5:
            a.doBlock(1000)
        for i in xrange(3): a.doBlock(1000)        
        return a.getData(), a.getLayout()
    return run_in_separate_process(all)




    

def cubicShrink(i):
    def all():        
        from myopenmm import Simulation
        a = Simulation(timestep = 10, thermostat = 0.5,name = "shrink")
        a.setup()
        a.load("/home/magus/evo/trajectories/globule_creation/32000_SAW_expanded/crumpled%d.dat" % i )
        a.setLayout("chain")

        a.addGrosbergPolymerBonds()
        a.addGrosbergRepulsiveForce()
        a.addSphericalConfinement(density = 0.95)
        a.energy_minimization(300)
        for j in xrange(20):
            a.doBlock(100)
        a.save("/home/magus/evo/trajectories/globule_creation/32000_SAW_compressed/crumpled%d.dat" % i )
    return run_in_separate_process(all)


        
    
        
        
        


    




def cubic_grid_example(data):
    from myopenmm import Simulation
    a = Simulation(timestep = 5, thermostat = 0.001)
    #a.load("/home/magus/evo/trajectories/saws/300/saw3.dat",center = "zero")
    #a.load("/home/magus/evo/GO43_test_solvent/expanded-1.dat",center = "zero")
    a.load(myrw)
    myN = len(a.getData())
    a.setup(PBC = True)
    a.verbose = True
    a.setLayout()
    a.add_cubic_grid(step = 0.6, percent = 0.5, mass = 0.2)
    a.addHarmonicPolymerBonds(dist = 0.1)
    a.addLennardJonesForce(epsilonRep = 0.01,sigmaAttr = 1,sigmaRep = 0.6,epsilonAttr = 1.2)
    
    #a.addSphericalConfinement()
    a.saveFolder("/home/magus/evo/GO43_test_solvent")
    save(a.getScaledData()[:myN],a.folder + "/expanded%d.dat" % -1)
    a.energy_minimization(steps = 30)
    a.massmult = 0.2
    a.initVelocities()
    
    for i in xrange(200000):
        a.doBlock(300 ,num = 100)
        save(a.getScaledData()[:myN],a.folder + "/expanded%d.dat" % i)

        
     
#diffuse_with_three_contacts(None)
    

def do_brush():
    from myopenmm import Simulation
    a = Simulation(timestep = 25, thermostat = 0.005)
    a.load("/home/magus/evo/trajectories/RWs/32000/rw1.dat")
    ##a.load("/home/magus/evo/GO43_test_solvent/expanded-1.dat",center = "zero")
    #a.load(myrw)
    #myN = len(a.getData())
    a.setup()
    a.verbose = True
    a.setLayout()
  
    a.addHarmonicPolymerBonds(dist = 0.1)
    a.addSimpleRepulsiveForce(trunk = 14)
    
    
    a.addConsecutiveRandomBonds(27, 0.6,0.,smeer = 0.06)
    a.addConsecutiveRandomBonds(1300, 1.5,0.,smeer = 0.06)
    #a.addSphericalConfinement()
    a.saveFolder("/home/magus/evo/GO_test")
    #save(a.getScaledData()[:myN],a.folder + "/expanded%d.dat" % -1)
    a.energy_minimization(steps = 500)
    #a.massmult = 0.2
    for i in xrange(10):
        a.doBlock(20000)
        a.save()


#data = []
#base = -300.
#for i in xrange(100):
#    data.append((base,(-1)**i * sqrt(2.)/4.,0 ))
#    base += sqrt(2.)/2
#creates initial configuration without knot

def wall_attraction_RW():
    from myopenmm import Simulation
    a = Simulation(timestep = 30, thermostat = 0.2)
    a.load("data.xyz")
    a.saveFolder(".")
    a.setup(PBC = True)
    a.verbose = True
    a.setLayout()
    a.addHarmonicPolymerBonds(dist = 0.1)            
    
    a.useDomains(mode = "chromosomes",chromosomes = [2])
    print (a.domains == 1).sum()
  
    a.addLaminaAttraction(4,-.7,r=10)

    #a.addInteraction(1000,1200,0.3)
    #a.addInteraction(1800,2000,0.6)
    #a.addInteraction(2400,2600,0.9)
    #a.addInteraction(3000,3200,1.2)    
    #a.addInteraction(3600,3800,1.5)
    a.energy_minimization(steps = 200)
    for i in xrange(100):
        a.doBlock(3000,num = 2000)
        a.save()
     
    


def create_pure_RW():
    from myopenmm import Simulation
    a = Simulation(timestep = 20, thermostat = 0.0015)
    a.load("/home/magus/evo/trajectories/globule_creation/32000_RW/crumpled3.dat")
    a.saveFolder("/home/magus/evo/GO46_create_pure_rw")
    a.setup()
    a.verbose = True
    a.setLayout()
    a.addHarmonicPolymerBonds(dist = 0.1)
    a.energy_minimization()
    for i in xrange(20000):
        a.doBlock(50000,num = 10000)
        a.save()
        
    
#create_pure_RW()


def diffuse_with_three_contacts(data):
    from myopenmmNew import Simulation
    a = Simulation(timestep = 90, thermostat = 0.002)
    a.load("/home/magus/evo/trajectories/globule_creation/4000_new_protocol/globules/crumpled1.dat")
    a.saveFolder("/home/magus/evo/GO41_interaction_test/lamina")
    a.setup()
    a.verbose = True
    a.setLayout()
    a.addGrosbergPolymerBonds() 
    a.addGrosbergRepulsiveForce()  
    a.addSphericalConfinement(10)
    a.addStifness(8)
    #a.show(coloring = "domains")
    a.energy_minimization(500,twoStage = True)
    for i in xrange(100):
        a.doBlock(5000)
        a.save()
        a.show()
    
     
diffuse_with_three_contacts(None)     
    
#sausage = shrink(*create_sausage(4000))    
#diffuse_with_three_contacts(None)

    
def knotDiffusion(filename,savedir):
    "eats filename or 'xyz' file" 
    def run(filename = filename, savedir = savedir):
        from myopenmm import Simulation
        a = Simulation(timestep = 30, thermostat = 0.012)  #create a simulation object
        a.load(filename)   #loads particle data in the system
        
                
        a.setup()          #sets up all the variables
        
        a.verbose = True
        lengthes = [6000] * 5 + [3000] * 6 + [2000] * 6
        lengthes = numpy.array(lengthes)
        shift = 0
        chains = []
        for i in lengthes:
            chains.append((shift,shift + i))
            shift += i
        
        
            
            
        a.setLayout()      #tells the MD code that this is a chain - default system. Change it for rings.
         
        a.saveFolder(savedir)   #folder where trajectory will be saved
        a.addHarmonicPolymerBonds(dist = 0.08)   #adds bonds between neighboring monomers
        #a.addHarmonicStifnessBonds(divider = 2.5)   #increase divider to the larger value to make polymer less stiff
        #do not decrease divider less that 2.5  
        a.addLennardJonesForce( epsilonRep = 0.4)              #adds repulsion between particles
        #a.tetherParticles([0,a.N - 1])           #tethers ends
        #a.createWalls()                         #creates walls at the first and last particles to prevent knots from forming at the end
        
        a.applyForces()
                           
        a.energy_minimization(100,True,10)
        a.save("expanded0.dat")
        a.save("initial.dat")
        for i in xrange(5000):
            a.doBlock(2000)
            a.save()
    run()
    #run_in_separate_process(run)
    


#knotDiffusion("/home/magus/evo/trajectories/saws/32000/saw2.dat","/home/magus/evo/GO44_opening_out/mixed")
#xit()        
def run_opening_out():
    i = 1
    data = resolve("/home/magus/evo/GO38_32k_diffusion/eq/expanded%d.dat" % (1700 + 300 * i),steps = 2000)[0]
    openingOut(data,"/home/magus/evo/GO44_opening_out/equilibrium%d" % i)
    
    
    for i in xrange(2,11):
        openingOut("/home/magus/evo/GO20_test_knot/globules32/crumpled%d.dat" % i,"/home/magus/evo/GO44_opening_out/fractal%d" % i)
        expname = "/home/magus/evo/GO38_32k_diffusion/run%d/expanded24000.dat" % i
        if os.path.exists(expname):
            openingOut(expname,"/home/magus/evo/GO44_opening_out/mixed%d" % i )
        data = resolve("/home/magus/evo/GO38_32k_diffusion/eq/expanded%d.dat" % (1700 + 300 * i),steps = 4000)[0]
        openingOut(data,"/home/magus/evo/GO44_opening_out/equilibrium%d" % i)
from base import fmap
def check_opening_out():
    cr = []
    mix = []
    eq = []
    beg,end = (32000/5) * 2,(32000/5)*3
    for i in xrange(1,11):
        
        def rad(data):
            return numpy.sqrt(numpy.sum(numpy.var(data[beg:end] ,axis = 0)))
            
        cr.append (fmap(lambda j:rad(load("/home/magus/evo/GO44_opening_out/fractal%d/expanded%d.dat" % (i,j) )) ,range(1,51),n=8))
        eq.append (fmap(lambda j:rad(load("/home/magus/evo/GO44_opening_out/equilibrium%d/expanded%d.dat" % (i,j) )) ,range(1,51),n=8))
        try:mix.append (fmap(lambda j:rad(load("/home/magus/evo/GO44_opening_out/mixed%d/expanded%d.dat" % (i,j) )) ,range(1,51),n=8))
        except: print "no mixed globule for %d" % i
    import matplotlib.pyplot as plt
    plt.figure(figsize = (6,4))
    plt.plot(coolAverage(cr),label="Fractal",color = (0.23 , 0.32 , 0.95))
    plt.plot(coolAverage(mix),label = "Mixed" , color = (0.5 , 1 , 0.13))
    plt.plot(coolAverage(eq),label = "Equilibrium",color = (1 , 0.63 , 0.13))
    removeAxes(shift = 0)
    plt.xlabel("Time",fontsize = 17)
    plt.ylabel("Loop size" , fontsize = 17)       
    
    
    ax = plt.axes()
    for xlabel_i in ax.get_xticklabels():
        xlabel_i.set_fontsize(15)
    for xlabel_i in ax.get_yticklabels():
        xlabel_i.set_fontsize(15)    
    legend = plt.legend(prop={"size":15},loc = 0)
    legend.draw_frame(False)
    plt.gcf().subplots_adjust(left=0.07, bottom=0.12, top=0.98, right=0.98)
    
        
    plt.show()
        
        
#check_opening_out()        
        
        
    
 
    
#do_brush()    

