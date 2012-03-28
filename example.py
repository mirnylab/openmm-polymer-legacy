from math import sin,cos,sqrt
import numpy
from openmmlib import Simulation

        
def exampleOpenmm():
    """
    You need to have a OpenMM-compatible GPU and OpenMM installed to run this script. 
    Otherwise you can switch to "reference" platform 
    a.setup(platform = "reference")
    But this will be extremely slow... 
    
    Installing OpenMM may be not easy too... but you can try 
    """
    
    
    a = Simulation(timestep = 80, thermostat = 0.002)
    
    a.setup(platform = "OpenCL", verbose = True)
    a.saveFolder("trajectory")   #folder where to save trajectory   
    a.load("globule")  #filename to load
    a.setLayout(mode = "chain")                     #default = chain    
    a.addSphericalConfinement(density = 0.55)    
    a.addHarmonicPolymerBonds()        
    a.addGrosbergRepulsiveForce()   #Fastest pure repulsive force
    a.addGrosbergStiffness()
    a.energyMinimization(steps = 200,twoStage= True)
            
    for _ in xrange(10):
        a.doBlock(3000)        
        a.save()
    a.printStats()
    a.show()
    
    
exampleOpenmm()
exit() 

"""Below are scripts that can be used to create starting conformation. 
create_spiral will create a spiral. 
create_sausage will create a sausage with a certain length/width ratio of a certain N. 
""" 



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

        
    
    

def create_sausage(N,ratio = 4,enlarge = 1):        
    "creates sausage by dropping a spiral into the cylinder "
    D =( 5* N / (3. * ratio)) ** 0.33333 + 1
    D = D * enlarge
    print D
    print ratio * D
    spiral = create_spiral(3.1,3.3,N)
    import openmmlib        
    a = openmmlib.Simulation(timestep = 40, thermostat =0.1,name = "sausage")
    a.setup()
    a.load(spiral)
    a.setLayout(mode = "chain")
    a.addCylindricalConfinement(r=D/2.,bottom = True,k=2.5)
    a.addHarmonicPolymerBonds(0.06)
    a.addStiffness(5)


    a.addSimpleRepulsiveForce()
    a.addGravity(0.3,cutoff = D * ratio )    
    a.steps_per_block = 4000
    a.energyMinimization(steps = 100)
    
    while a.getData()[:,2].max() >  D * ratio + 15:
        #print a.getData()[:,2].max(),  a.getData()[:,1].min(),             (ratio * D )
        a.doBlock()
        #a.initVelocities()
        #a.show()
    
    for _ in xrange(5):
        a.doBlock()
    a.show()
    
    


    
create_sausage(1000,3)
