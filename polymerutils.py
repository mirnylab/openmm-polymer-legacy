from mirnylib import numutils  
import numpy 
from scipy import weave
from math import sqrt
from mirnylib.systemutils import fmapred,fmap , deprecate
import joblib
import sys 
import mirnylib
import os 
from mirnylib.h5dict import h5dict

def load(filename, h5dictKey = None):
    """Universal load function for any type of data file"""
    
    if not os.path.exists(filename):
        raise IOError("File not found :( ")
    
                 
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
    


