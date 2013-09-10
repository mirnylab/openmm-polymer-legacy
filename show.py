#(c) 2013 Massachusetts Institute of Technology. All Rights Reserved
# Code written by: Maksim Imakaev (imakaev@mit.edu)


#!/usr/bin/env python
import numpy
import os
import tempfile
import sys
import subprocess
import joblib

def showData(data, rotate=(0, 0, 0), links=None, drawBonds=True):
    #if you want to change positions of the spheres along each segment, change these numbers
    #e.g. [0,.1, .2 ...  .9] will draw 10 spheres, and this will look better
    if drawBonds:
        shifts = [0., 0.2, 0.4, 0.6, 0.8]
    else:
        shifts = [0.]

    #determining the 95 percentile distance between particles,
    meandist = numpy.percentile(numpy.sqrt(
        numpy.sum(numpy.diff(data, axis=0) ** 2, axis=1)), 95)
    #rescaling the data, so that bonds are of the order of 1. This is because rasmol spheres are of the fixed diameter.
    data /= meandist

    #writing the rasmol script. Spacefill controls radius of the sphere.
    rascript = open(tempfile.NamedTemporaryFile().name, 'w')
    script = ("""set write on
        wireframe off
        color temperature
        spacefill 100
        background white
        rotate x {0}
        rotate y {1}
        rotate z {1}

        """.format(rotate[0], rotate[1], rotate[2]))
    rascript.write(script)
    rascript.flush()

    #creating the array, linearly chanhing from -225 to 225, to serve as an array of colors
    #(rasmol color space is -250 to 250, but it  still sets blue to the minimum color it found and red to the maximum).
    colors = numpy.array([int(
        (j * 450.) / (len(data))) - 225 for j in xrange(len(data))])

    #creating spheres along the trajectory
    #for speedup I just create a Nx4 array, where first three columns are coordinates, and fourth is the color
    newData = numpy.zeros((len(data) * len(shifts) - (len(shifts) - 1), 4))
    for i in xrange(len(shifts)):
        #filling in the array like 0,5,10,15; then 1,6,11,16; then 2,7,12,17, etc.
        #this is just very fast
        newData[i:-1:len(shifts), :3] = data[:-1] * shifts[i] + \
            data[1:] * (1 - shifts[i])
        newData[i:-1:len(shifts), 3] = colors[:-1]
    newData[-1, :3] = data[-1]
    newData[-1, 3] = colors[-1]

    #inserting links
    if not(links is None):
        if issubclass(type(links), numpy.ndarray):
            linksData = numpy.zeros((links.shape[0] * len(shifts) + 1, 4))
        elif issubclass(type(links), list) and issubclass(type(links[0]), tuple):
            linksData = numpy.zeros((len(links) * len(shifts) + 1, 4))
            links = numpy.array(links)
        else:
            raise Exception('Unknown format of the links')

        for i in xrange(len(shifts)):
            linksData[i:-1:len(shifts), :3] = (data[links[:, 0]] * shifts[i] +
                data[links[:, 1]] * (1 - shifts[i]))
            linksData[i:-1:len(shifts), 3] = colors[-1]

        newData = numpy.vstack([linksData[:-1], newData])

    towrite = open(tempfile.NamedTemporaryFile().name, 'w')
    towrite.write("%d\n\n" % (len(newData)))
        #number of atoms and a blank line after is a requirement of rasmol

    for i in newData:
        towrite.write("CA\t%lf\t%lf\t%lf\t%d\n" % tuple(i))
    towrite.flush()
    #For windows you might need to change the place where your rasmol file is
    if os.name == "posix":  # if linux
        subprocess.Popen(
            "rasmol -xyz {0} -script {1}; rm {0}; rm {1}".format(towrite.name, rascript.name),
            shell=True)
    else:  # if windows
        os.system("C:/RasWin/raswin.exe -xyz %s -script %s" % (
            towrite.name, rascript.name))
    #exit()


def load(filename):
    try:
        return joblib.load(filename)["data"]
    except (IOError, KeyError):

        line0 = open(filename).readline()
        N = int(line0)
        lines = open(filename).readlines()[1:]
        data = [[float(i) for i in j.split()] for j in lines if len(j) > 3]

        if len(data) != N:
            raise StandardError(
                "N does not correspond to the number of lines!")
        return data

if __name__ == '__main__':
    if len(sys.argv) == 3:
        print "Assuming h5dict file first"
        try:
            from mirnylib.h5dict import h5dict
            data = h5dict(path=sys.argv[1], mode="r")[sys.argv[2]]
            showData(data)
            sys.exit()
        except IOError:
            print "failed to load h5dict file, trying regular file"


    if len(sys.argv) > 2:
        showData(load(sys.argv[1]), (sys.argv[2], sys.argv[3], sys.argv[4]))
    try:
        showData(load(sys.argv[1]))
        sys.exit()

    except:
        showData(load("block%s.dat" % sys.argv[1]))



#show3D(numpy.cumsum(numpy.random.randint(-1,2,(3,10000)),axis = 1))  #an example
