import os
import joblib
import numpy
import sys


infile = sys.argv[1]
try:
    outfile = sys.argv[2]
except:
    outfile = sys.argv[1]
    
lines = open(infile).readlines()[1:]
data = [[float(i) for i in j.split()] for j in lines if len(j) > 3]

data2 = numpy.array(data,dtype = float)
ourdict = {"data":data2,"OldFilename":os.path.abspath(infile)}
joblib.dump(ourdict,outfile)