import os
import joblib
import numpy
import sys

"""
Sample txt file looks like this
next line starts the file
10
x1 y1 z1
x2 y2 z2
....
x10 y10 z10

So the file has the first line with a number of atoms, and then 10 lines with their coordinates
"""
infile = sys.argv[1]
try:
    outfile = sys.argv[2]
except:
    outfile = sys.argv[1]

lines = open(infile).readlines()[1:]
data = [[float(i) for i in j.split()] for j in lines if len(j) > 3]

data2 = numpy.array(data, dtype=float)
ourdict = {"data": data2, "OldFilename": os.path.abspath(infile)}
joblib.dump(ourdict, outfile, compress=3)
