import numpy 
import joblib
import os 
from math import sqrt, sin, cos  
from mirnylib.plotting import showPolymerRasmol


def load(filename, h5dictKey = None):
    """Universal load function for any type of data file"""
    
    if not os.path.exists(filename):
        raise IOError("File not found :( \n %s" % filename)
    
                 
    try:
        "checking for a text file"
        line0 = open(filename).readline() 
        try: N = int(line0)
        except ValueError: raise TypeError("Cannot read text file... reading pickle file")                                        
        lines = open(filename).readlines()[1:]
        data = [[float(i) for i in j.split()] for j in lines if len(j) > 3]
        
        if len(data) != N:
            raise ValueError("N does not correspond to the number of lines!")
        return numpy.array(data)
    
    except TypeError:
        pass 
    
    try:              
        "loading from a joblib file here"                
        mydict = dict(joblib.load(filename))
        data = mydict.pop("data")
        return data
    
    except:
        pass         

         

    try:
        "checking for h5dict file "
        from mirnylib.h5dict import h5dict                     
        mydict = h5dict(path = filename,mode = 'r')
        if h5dictKey == None: 
            keys = mydict.keys() 
            if len(keys) != 1: 
                raise ValueError("H5Dict has more than one key. Please specify the key.")
            h5dictKey = keys[0]
        assert h5dictKey in mydict.keys() 
        data = mydict[str(h5dictKey)]
        return data
    except IOError:
        raise IOError("Failed to open file") 
        
    

         
def save(data, filename, mode = "txt", h5dictKey = "1" ):
    
    h5dictKey = str(h5dictKey)
    mode = mode.lower()
    
    if mode == "h5dict":
        from mirnylib.h5dict import h5dict
        mydict = h5dict(filename, mode = "w")
        mydict[h5dictKey] = data                   
        return
    
    
    if mode == "joblib":
        metadata = {}
        metadata["data"] = data
        joblib.dump(metadata,filename = filename,compress = 3)
        return 
    
    if  mode == "txt":        
        lines = [str(len(data)) + "\n"]
        for particle in data: 
            lines.append("".join([str(j) + " " for j in particle]) + "\n")            
        with open(filename,'w') as myfile:
            myfile.writelines(lines)
        return 
    
    raise ValueError("Unknown mode : %s, use h5dict, joblib or txt" % mode)


def generateRandomLooping(length = 10000, oneMoverPerBp = 1000 ,numSteps = 100):
    
    
    N = length  
    myarray = numpy.zeros(N,int)    
    movers = []
    onsetRate = length / float(oneMoverPerBp)
    
    def initMovers():
        for i in movers:
            myarray[i[0]] = 1
            myarray[i[1]] = 1
    
    def addMovers():
        for _ in xrange(numpy.random.poisson(onsetRate)):
            pos = numpy.random.randint(N-1)
            if myarray[pos:pos+2].sum() == 0:
                movers.append((pos,pos+1))
                myarray[pos:pos+2] = 1
        
    def translocateMovers():
        moved = False
        for j,mover in enumerate(movers): 
            left,right = mover             
            if left > 0:
                if myarray[left-1] == 0:                                         
                    myarray[left] = 0
                    myarray[left-1] = 1 
                    left = left - 1
                    moved = True 
            if right < N-1:
                if myarray[right + 1] == 0:                    
                    myarray[right] = 0
                    myarray[right+1] = 1
                    right = right + 1
                    moved = True
            movers[j] = (left,right)
        return moved 
    
    for _ in xrange(numSteps):
        addMovers()
        translocateMovers()
    while translocateMovers():
        pass
    return movers 
             
            
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



                
    
         
def _test():
    
    print "testing save/load"    
    a = numpy.random.random((20000,3))
    save(a,"bla", mode = "txt")
    b = load("bla")
    assert abs(b.mean() - a.mean()) < 0.00001
    
    save(a,"bla", mode = "joblib")
    b = load("bla")
    assert abs(b.mean() - a.mean()) < 0.00001
    
    save(a,"bla", mode = "h5dict")
    b = load("bla")
    assert abs(b.mean() - a.mean()) < 0.00001
    
    os.remove("bla")
    
    print "Finished testing save/load, successfull"
    
#_test()
    


