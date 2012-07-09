import mirnylib.systemutils
from mirnylib.systemutils import fmap, fmapav, setExceptionHook
import polymerutils
mirnylib.systemutils.setExceptionHook() 
from mirnylib.numutils import logbins
from mirnylib.h5dict import h5dict
import os   
import contactmaps 
from contactmaps import Cload, giveContacts
from polymerutils import load 
import cPickle  
from math import exp,sqrt,log 

import numpy

import matplotlib.pyplot as plt 

from copy import copy


  

def giveCpScaling(data, bins0, cutoff=1.1,integrate = False,ring=False,intContacts = False):
    """
    Returns contact probability scaling for a given polymer conformation
    
    Parameters
    ----------
    data : 3xN array of ints/floats
        Input polymer conformation
    bins0 : list
        Bins to calculate scaling. 
    cutoff : float, optional 
        Cutoff to calculate scaling
    integrate : bool, optional 
        if True, will return cumulative probability
    ring : bool, optional 
        If True, will calculate contacts for the ring
    intContacts : bool, optional 
        If True, will speed up calculation of contacts for a cubit lattice case.
        
    Returns
    -------
    (bins, contact probabilities) where "bins" are centers of bins in the input bins0 array.  
    
    """
                 
    N = len(data[0])
    bins0 = numpy.array(bins0)          
    bins = [(bins0[i], bins0[i + 1]) for i in xrange(len(bins0) - 1)]    
    
    
    if intContacts == False:   contacts = numpy.array(giveContacts(data,cutoff))      
    else:  contacts = contactmaps.giveIntContacts(data)           #integer contacts are faster
    
    contacts = contacts[:,1] - contacts[:,0]   #contact lengthes
    
    if ring == True:
        mask = contacts > N/2         
        contacts[mask] = N - contacts[mask]       
    scontacts = numpy.sort(contacts)   #sorted contact lengthes         
    connections = 1.*numpy.diff(numpy.searchsorted(scontacts, bins0,side = "left"))   #binned contact lengthes
    possible = numpy.diff(N * bins0 + 0.5 * bins0 - 0.5 * (bins0**2)) 
    print "average contacts per monomer:", connections.sum() / N
    
    if integrate == False:  connections /= possible
    if integrate == True: connections = numpy.cumsum(connections)/connections.sum() 

    a = [sqrt(i[0] * i[1]) for i in bins]
    print list(connections) 
    return (a, connections)


def give_distance(data, bins=None,ring = False ):
    """
    Returns end-to-end distance scaling of a given polymer conformation. 
    ..warning:: This method averages end-to-end scaling over bins to make better average

    Parameters
    ----------
    
    data: 3xN array 
    
    """
    N = len(data[0])
    if ring == True:
        data = numpy.concatenate([data,data],axis = 1)
    bins = [(bins[i], bins[i + 1]) for i in xrange(len(bins) - 1)]
    
    rads = [0. for i in xrange(len(bins))]  
    for i in xrange(len(bins)):
        oneBin = bins[i] 
        rad = 0.
        count = 0
        for j in xrange(oneBin[0],oneBin[1], (oneBin[1]-oneBin[0])/10 + 1):
            length = j 
            if ring == True: 
                rad += numpy.mean(numpy.sqrt(numpy.sum((data[:,:N]-data[:,length:length+N])**2,0)))
            else: 
                rad += numpy.mean(numpy.sqrt(numpy.sum((data[:,:-length]-data[:,length:])**2,0))) 
            count += 1
         
        rads[i] = rad/count
    bins = [sqrt(i[0] * i[1]) for i in bins] 
    return (bins, rads) 



def give_radius_scaling(data, bins=None,ring = False):
    "main working horse for radius of gyration"
    "uses dymanic programming algorithm"
    
    bins = [sqrt(bins[i]* bins[i + 1]) for i in xrange(len(bins) - 1)]    
        
    data = numpy.array(data,float)
    coms = numpy.cumsum(data,1)   #cumulative sum of locations to calculate COM
    coms2 = numpy.cumsum(data**2,1)  #cumulative sum of locations^2 to calculate RG    
    def radius_gyration(len2):
        data 
        if ring == True:
            comsadd = coms[:,:len2].copy()
            coms2add = coms2[:,:len2].copy()
            comsadd += coms[:,-1][:,None]
            coms2add += coms2[:,-1][:,None]
            comsw = numpy.concatenate([coms,comsadd],axis = 1)
            coms2w = numpy.concatenate([coms2,coms2add],axis = 1)  #for rings we extend need longer chain
        else:
            comsw = coms
            coms2w = coms2
            
        coms2d = (-coms2w[:,:-len2]+coms2w[:,len2:])/len2
        comsd = ((comsw[:,:-len2]-comsw[:,len2:])/len2)**2
        diffs = coms2d - comsd
        sums = numpy.sqrt(numpy.sum(diffs,0))
        return numpy.mean(sums)
        
    rads = [0. for i in xrange(len(bins))]
    for i in xrange(len(bins)):
        rads[i] = radius_gyration(int(bins[i]))
    return (copy(bins), rads)



def give_radius_scaling_eig(data, bins=None):
    #gives volume^^0.33 as  defined through  eigenvectors
    if bins==None: bins = [2*i for i in logbins(1,0.45*len(data[0]),1.3,40)]
    x,y,z = data[0],data[1],data[2]
    coords = [x,y,z]
    sums = [[i*j for j in coords] for i in coords]
    sums = numpy.array(sums) 
    N = len(data[0])
    coms = numpy.cumsum(data,1)
    sums = numpy.cumsum(sums,2)
    def tensor(a,b):        
        newsums = (sums[:,:,b]-sums[:,:,a])/float(b-a)
        newcoms = (coms[:,b] - coms[:,a])/float(b-a)
        tensor =  newsums - numpy.outer(newcoms,newcoms)
        return numpy.linalg.eigvalsh(tensor)
    ret = []
    for i in bins:
        av = 0.
        for j in xrange(1000):
            t = numpy.random.randint(5,N-5-i)
            res = tensor(t,t+i)
            av += sqrt(3) * (res[0] * res[1] * res[2] * 1. )**(1/6.) 
        ret.append(av/1000)
    retret = (copy(bins),ret)
    return  retret


def subchainDensityFunction(filenames,bins,normalize = "Rg",lengthmult  = 3, Nbins = 30,coverage = 1 ):
    "Calculates density function of subchains"
    
    results = []
    
    for filename in filenames: 
        newbins = zip(bins[::2],bins[1::2]) #bin start/end
         
        
        data = polymerutils.load(filename)
        N = len(data) 
         
        rgs = give_radius_scaling(data.T, bins, ring = False)[1][::2]          
        curresults = []
        labels = []
        for bin,rg in zip(newbins,rgs):
            labels.append("S = %.1lf; " % (numpy.mean(bin)) + " Rg=%.2lf" % rg)
            lengthbins = numpy.linspace(0,lengthmult * rg,Nbins)
            lengthBinMids = (lengthbins[:-1] + lengthbins[1:])*0.5
            volumes = (4./3.) * 3.141592 * (lengthbins**3) 
            volumes = numpy.diff(volumes)
            count = int(N * coverage / numpy.mean(bin) + 1)
            sphereCounts = numpy.zeros(len(volumes),float)
            
            for i in xrange(count):
                size = numpy.random.randint(bin[0],bin[1])
                start = numpy.random.randint(0,N-size)
                subchain = data[start:start+size]
                com = numpy.mean(subchain,axis = 0)
                shifted = subchain - com[None,:]
                dists = numpy.sqrt(numpy.sum(shifted**2,axis = 1))
                sphereCounts += numpy.histogram(dists, lengthbins)[0]
            sphereCounts /= (volumes * count)
            #curresults.append(numpy.array([lengthBinMids/rg,sphereCounts]))
            curresults.append(numpy.array([lengthBinMids,sphereCounts]))
        results.append(curresults)
    results = numpy.mean(numpy.array(results,float),axis = 0)
    for i,label in zip(results,labels):
        plt.plot(*i,label = label)
    
                 
    
    
    plt.legend()
      
        
            
            
        
#filenames = ["/home/magus/evo/trajectories/globule_creation/32000_SAW/crumpled%d.dat" % i for i in xrange(1,20)]
filenames = ["/home/magus/evo/topo37_256_battery/run%d/expanded%d.dat" % (i,j) for i in xrange(1,11) for j in [1000,1030,1060,1090]]


"""
filenames = ["/home/magus/evo/GO38_32k_diffusion/run%d/expanded%d.dat" % (i,j) for i in xrange(2,10) for j in [900,1000,10]]
subchainDensityFunction(filenames,[1000,1200])
filenames = ["/home/magus/evo/GO38_32k_diffusion/run%d/expanded%d.dat" % (i,j) for i in xrange(2,10) for j in [0]]
subchainDensityFunction(filenames,[1000,1200])
filenames = ["/home/magus/HiC2011/openmm_simulations/02_nechaev_corsslinks/trajectory%d/block%d.dat" % (i,j) for i in xrange(2,5) for j in range(100,200,10)]
subchainDensityFunction(filenames,[1000,1200])

filenames = ["/home/magus/evo/trajectories/globule_creation/32000_equilibrium/equilibrium%d.dat" % (i) for i in xrange(2,40)]
subchainDensityFunction(filenames,[1000,1200])


plt.show() 
"""
plt.subplot(121)
plt.title("Mixed globule")
filenames = ["/home/magus/evo/topo37_256_battery/run%d/expanded%d.dat" % (i,j) for i in xrange(1,11) for j in xrange(900,1000,10)]
subchainDensityFunction(filenames,[100,120,1000,1200,10000,12000])

plt.subplot(122)
plt.title("Equilibrium globule")
filenames = ["/home/magus/evo//topo371_256_equilibrium/run%d/expanded%d.dat" % (i,j) for i in xrange(1,11) for j in xrange(130,170,3)]
subchainDensityFunction(filenames,[100,120,1000,1200,10000,12000])

plt.show() 


plt.subplot(121)
plt.title("Mixed globule")
filenames = ["/home/magus/evo/GO38_32k_diffusion/run%d/expanded%d.dat" % (i,j) for i in xrange(2,10) for j in [1000]]
subchainDensityFunction(filenames,[10,15,100,120,1000,1200,10000,10200])

plt.subplot(122)
plt.title("Crumpled globule")
filenames = ["/home/magus/evo/GO38_32k_diffusion/run%d/expanded%d.dat" % (i,j) for i in xrange(2,10) for j in [0]]
subchainDensityFunction(filenames,[10,15,100,120,1000,1200,10000,10200])


plt.show()

exit()
        
         
        
    
    
    





def give_slices(base, tosave, slices, sliceParams, multipliers, mode = "chain", loadFunction = Cload, integrate = False,normalize=False, exceptionList = [],nproc=4):
    numpy.seterr(invalid='raise')
    
    plotsBySlice = []

    for cur_slice in slices:
         
        files = []

        def slice2D(a, b,mult=[1]):             
            tm = []
            if type(b) == tuple:
                for i in xrange(b[0], b[1]+1):
                    tm.append((i, a))
            elif type(b) == int:
                for i in xrange(1, b + 1):
                    tm.append((i, a))
            elif type(b) == list: tm=[(i,a) for i in b]
            tm2 = sorted(list(set([(i[0],int(float(i[1])*m)) for i in tm for m in mult])))
            print tm2
            return tm2  
        
        def slice3D(a, b,c,mult=[1]):
            tm = []
            for i in xrange(b[0], b[1]+1):
                for t in xrange(c[0],c[1]+1):
                    tm.append((i, a,t))
            tm2 = sorted(list(set([(i[0],int(float(i[1])*m)) for i in tm for m in mult])))
            print tm2
            return tm2   
        
        #sluces actually are defined
        runs = slice2D(cur_slice, sliceParams ,multipliers)
        #runs = slice3D(cur_slice, (1,14),(1,10),multipliers)
        
        for i in runs:            
            #filename is replaced in slices
            try: files.append(base.replace("DATA1", str(i[0])).replace("DATA2", str(i[1])).replace("DATA3",str(i[2])))
            except: files.append(base.replace("DATA1", str(i[0])).replace("DATA2", str(i[1])))
            
    
        datas = []
        
        def newload(i):
            #loads a file  
            try:                    
                data =  loadFunction(i,False)
                if len(data) !=3:
                    data = data.T
                if len(data) != 3:
                    raise StandardError("Wrong shape of data")
                data = numpy.asarray(data,order = "C", dtype = float )                
                return data
            except exceptionList: 
                print "file not found" , i
                return None
            
        
        #use this for determining the file size 
        datas = filter(lambda x:x!=None,fmap(newload,files[::len(files)/20 + 1],n=3))
        datlen = len(datas[0][0])
        
        if mode == "chain": bins2 =  logbins(4, datlen-100,1.25)        
        if mode == "parts": bins2 =  logbins(4, datlen-100,1.25)
        if (mode == "ring") or (mode == "intring"): 
            b1 = logbins(2, datlen/4-1,1.34)
            bins2 =  [2*i for i in b1]
            print bins2
        binsrg = logbins(4, datlen-100,1.25)
        
        def give_plots(i):
            data = newload(i)
            if data == None: return None                  
            i = data                  

            if (mode == "ring") or (mode == "intring"):
                b = give_radius_scaling(i,binsrg,ring = True)
            else: 
                b = give_radius_scaling(i,binsrg,ring = False)
            
            
            if (mode == "chain") : a = giveCpScaling(i, bins2,1.7,integrate)
            if (mode == "ring"): a = giveCpScaling(i, bins2,1.7,integrate,ring = True)
            if (mode == "intring"): a = giveCpScaling(i, bins2,1.7,integrate,ring = True,project=False,intContacts = True)             
            if (mode == "project"): a = giveCpScaling(i, bins2,1.450,integrate,project=True)
            
            if (mode == "ring") or (mode == "intring"):
                c = give_distance(i, bins2,ring = True)
            else: 
                c = give_distance(i, bins2,ring = False)
                
            
                     
            flatten = normalize
            if (flatten == True):
            #some weired corrections... 
                print "using normal correction!!!"
                a[1] = [a[0][i]*a[1][i] for i in xrange(len(a[0]))]
                b = list(b)
                c = list(c)
                b[1] = [(b[1][i]/(b[0][i]**0.333333)) for i in xrange(len(b[0]))]
                c[1] = [(c[1][i]/(c[0][i]**0.333333)) for i in xrange(len(c[0]))]
            elif flatten=="super":
                print "using super correction!!!"
                b = list(b)
                c = list(c)
                for i in b[0]:
                    e = (log(datlen) - log(b[0][i]))/(log(datlen) - log(1))
                    #print e
                    b[1][i] = b[1][i]  * ((2-2*e)/(2-e))**(1/3.)
                    e = (log(datlen) - log(c[0][i]))/(log(datlen) - log(1))
                    c[1][i] = c[1][i] * ((2-2*e)/(2-e))**(1/3.)                                     
            a = numpy.array(a,dtype = float)
            b = numpy.array(b,dtype = float)
            c = numpy.array(c, dtype = float)  

            
            return numpy.array([a,b,c])
                
        
        parPlots = fmap(give_plots,files,n=nproc)
        
        parPlots = filter(lambda x:x!=None,parPlots)        
        
        means = numpy.mean(parPlots,axis = 0)                
        plotsBySlice.append([means,{"slice":cur_slice}])
          
        
            

    if tosave != None: cPickle.dump(plotsBySlice,open(tosave,'wb'),-1)
    print "Finished!!!"
    return plotsBySlice


#give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run4eq/expandedDATA2.dat", 
#            tosave  = "data/DNA_conf/plots/paper_scalings/ring13eq", 
#             slices = [11000], sliceParams = (1,1), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)

#give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/runDATA1/expandedDATA2.dat",
#             tosave  = "data/DNA_conf/plots/paper_scalings/ring13", 
#             slices = [300,1000,70000], sliceParams = (3,4), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)


#give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run5eq/expandedDATA2.dat",
#             tosave  = "data/DNA_conf/plots/paper_scalings/chain32eq", 
#             slices = [7000], sliceParams = (1,1), multipliers = numpy.arange(0.5,1,0.0001), mode = "chain", loadFunction = Cload)

#give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run5/expandedDATA2.dat",
#             tosave  = "data/DNA_conf/plots/paper_scalings/chain32", 
#             slices = [300,1000,30000], sliceParams = (1,1), multipliers = numpy.arange(0.5,1,0.0001), mode = "chain", loadFunction = Cload)
#  
#give_slices(base = "/home/magus/evo/topo371_256_equilibrium/runDATA1/expandedDATA2.dat",
#             tosave  = "data/DNA_conf/plots/paper_scalings/ring_256_eq", 
#             slices = [2,3,7,13], sliceParams = (1,10), multipliers = numpy.arange(0.7,1,0.0001), mode = "intring", loadFunction = intload)

#-------- equilibrium Rings low density -----------

#give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run6_lowden/expandedDATA2.dat",
#             tosave  = "data/DNA_conf/plots/paper_scalings/ring13_lowd", 
#             slices = [6200], sliceParams = (3), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)


#give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run6_eq/expandedDATA2.dat",
#             tosave  = "data/DNA_conf/plots/paper_scalings/ring13_lowd_eq", 
#             slices = [4600], sliceParams = (3), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)

#---------------ETE distance for N=108000 ----------------
#give_slices(base = "/home/magus/evo/topo36_108_grow/runDATA1/expandedDATA2.dat",
#             tosave  = "data/DNA_conf/plots/paper_scalings/ring_108_eq", 
#             slices = [20000], sliceParams = (1,10), multipliers = numpy.arange(0.7,1,0.01), mode = "intring", loadFunction = intload)

#---------------------------Equilibrium rings, 4k ---------
"""give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run7_small/expandedDATA2.dat",
             tosave  = "/home/magus/workspace/testnucl2/data/DNA_conf/plots/paper_scalings/ring4", 
             slices = [4000], sliceParams = (3), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)

give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run7_eq/expandedDATA2.dat",
             tosave  = "data/DNA_conf/plots/paper_scalings/ring4_eq", 
             slices = [7000], sliceParams = (3), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)

give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run8_tiny/expandedDATA2.dat",
             tosave  = "data/DNA_conf/plots/paper_scalings/ring2", 
             slices = [14000], sliceParams = (3), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)


give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run8_tiny_eq/expandedDATA2.dat",
             tosave  = "data/DNA_conf/plots/paper_scalings/ring2_eq", 
             slices = [9000], sliceParams = (3), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)
"""


class h5dictLoad(object):
    """
    An experimental class to fetch h5dict values based on a fake filename. 
    
    It accepts filenames in a form
    /path-to-h5dict/h5dictKey
    
    It automatically caches h5dicts that are already open
    
    use simpleFetch method for concurrent fetching using fmap
    """
    def __init__(self):
        self.baseDict = {}
    def fetch(self,filename, dummy = True ):
        base,num = os.path.split(filename)
                    
        if base not in self.baseDict.keys():  
            self.baseDict[base] = h5dict(base,mode = 'r')
        return self.baseDict[base][num]
    def simpleFetch(self,filename,dummy = True):
        base,num = os.path.split(filename)
        h5 = h5dict(base,mode = "r")
        toret = h5[num]
        del h5
        return toret
        
#import joblib


def _test():
    def giveStraightChain(x,y):
        N = 4000
        a = numpy.zeros((N,3),float)
        a[:,0] = numpy.arange(N)
        return a 
    
    scalings = give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run8_tiny_eq/expandedDATA2.dat",
                 tosave  = None, 
                 slices = [9000], sliceParams = (3), multipliers = numpy.arange(0.5,1,0.01), mode = "chain", loadFunction = giveStraightChain)
    plt.plot(*scalings[0][0][0])
    plt.show()
    plt.plot(*scalings[0][0][1])
    plt.show()
    plt.plot(*scalings[0][0][2])
    plt.show()



#scalings = give_slices(base = "/home/magus/HiC2011/openmm_simulations/02_nechaev_corsslinks/blockDATA2.dat",
#             tosave  = None, 
#             slices = [100,500], sliceParams = (3), multipliers = numpy.arange(0.5,1,0.0001), mode = "chain", loadFunction = (lambda x,y:joblib.load(x)["data"]) )


#setExceptionHook()
#plt.plot(*scalings[0][0][0])
#plt.show()
#plt.plot(*scalings[0][0][1])
#plt.show()
#plt.plot(*scalings[0][0][2])
#plt.show()
#exit()
