import numpy as np
import joblib
import os
from math import sqrt, sin, cos

from mirnylib.numutils import isInteger, rotationMatrix

import numpy
from mirnylib.plotting import showPolymerRasmol



def load(filename, h5dictKey=None):
    """Universal load function for any type of data file"""

    if not os.path.exists(filename):
        raise IOError("File not found :( \n %s" % filename)

    try:
        "checking for a text file"
        line0 = open(filename).readline()
        try:
            N = int(line0)
        except ValueError:
            raise TypeError("Cannot read text file... reading pickle file")
        lines = open(filename).readlines()[1:]
        data = [[float(i) for i in j.split()] for j in lines if len(j) > 3]

        if len(data) != N:
            raise ValueError("N does not correspond to the number of lines!")
        return np.array(data)

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
        mydict = h5dict(path=filename, mode='r')
        if h5dictKey is None:
            keys = mydict.keys()
            if len(keys) != 1:
                raise ValueError("H5Dict has more than one key. Please specify the key.")
            h5dictKey = keys[0]
        assert h5dictKey in mydict.keys()
        data = mydict[str(h5dictKey)]
        return data
    except IOError:
        raise IOError("Failed to open file")

def save(data, filename, mode="txt", h5dictKey="1"):
    h5dictKey = str(h5dictKey)
    mode = mode.lower()

    if mode == "h5dict":
        from mirnylib.h5dict import h5dict
        mydict = h5dict(filename, mode="w")
        mydict[h5dictKey] = data
        del mydict
        return

    elif mode == "joblib":
        metadata = {}
        metadata["data"] = data
        joblib.dump(metadata, filename=filename, compress=3)
        return

    elif mode == "txt":
        lines = [str(len(data)) + "\n"]
        for particle in data:
            lines.append("".join([str(j) + " " for j in particle]) + "\n")
        with open(filename, 'w') as myfile:
            myfile.writelines(lines)
        return

    elif mode == 'pdb':
        data = data - np.min(data, axis=0)[None, :]

        N = len(data)
        retret = ""

        def add(st, n):
            if len(st) > n:
                return st[:n]
            else:
                return st + " " * (n - len(st))

        for i, line in enumerate(data):
            line = [1. * (float(j) + 1) for j in line]
            ret = add("ATOM", 7)
            ret = add(ret + "%i" % (i + 1), 13)
            ret = add(ret + "CA", 17)
            ret = add(ret + "ALA", 22)

            ret = add(ret + "%i" % (i), 30)
            ret = add(ret + ("%8.3f" % line[0]), 37)
            ret = add(ret + ("%8.3f" % line[1]), 45)
            ret = add(ret + ("%8.3f" % line[2]), 53)
            ret = add(ret + (" 1.00"), 61)
            ret = add(ret + str(float(i % 8 > 4)), 67)
            retret += (ret + "\n")

        f = open(filename, 'w')
        f.write(retret)
        f.flush()

    else:
        raise ValueError("Unknown mode : %s, use h5dict, joblib, txt or pdb" % mode)

def generateRandomLooping(length=10000, oneMoverPerBp=1000, numSteps=100):

    N = length
    myarray = np.zeros(N, int)
    movers = []
    onsetRate = length / float(oneMoverPerBp)

    def initMovers():
        for i in movers:
            myarray[i[0]] = 1
            myarray[i[1]] = 1

    def addMovers():
        for _ in xrange(np.random.poisson(onsetRate)):
            pos = np.random.randint(N - 1)
            if myarray[pos:pos + 2].sum() == 0:
                movers.append((pos, pos + 1))
                myarray[pos:pos + 2] = 1

    def translocateMovers():
        moved = False
        for j, mover in enumerate(movers):
            left, right = mover
            if left > 0:
                if myarray[left - 1] == 0:
                    myarray[left] = 0
                    myarray[left - 1] = 1
                    left = left - 1
                    moved = True
            if right < N - 1:
                if myarray[right + 1] == 0:
                    myarray[right] = 0
                    myarray[right + 1] = 1
                    right = right + 1
                    moved = True
            movers[j] = (left, right)
        return moved

    for _ in xrange(numSteps):
        addMovers()
        translocateMovers()
    while translocateMovers():
        pass
    return movers


def create_spiral(r1, r2, N):
    Pi = 3.141592
    points = []
    finished = [False]

    def rad(phi):
        return phi / (2 * Pi)

    def ang(rad):
        return 2 * Pi * rad

    def coord(phi):
        r = rad(phi)
        return (r * sin(phi), r * cos(phi))

    def fullcoord(phi, z):
        c = coord(phi)
        return [c[0], c[1], z]

    def dist(phi1, phi2):
        c1 = coord(phi1)
        c2 = coord(phi2)
        d = sqrt((c1[1] - c2[1]) ** 2 + (c1[0] - c2[0]) ** 2)
        return d

    def nextphi(phi):
        phi1 = phi
        phi2 = phi + 0.7 * Pi
        mid = phi2
        while abs(dist(phi, mid) - 1) > 0.00001:
            mid = (phi1 + phi2) / 2.
            if dist(phi, mid) > 1:
                phi2 = mid
            else:
                phi1 = mid
        return mid

    def prevphi(phi):

        phi1 = phi
        phi2 = phi - 0.7 * Pi
        mid = phi2

        while abs(dist(phi, mid) - 1) > 0.00001:
            mid = (phi1 + phi2) / 2.
            if dist(phi, mid) > 1:
                phi2 = mid
            else:
                phi1 = mid
        return mid

    def add_point(point, points=points, finished=finished):
        if (len(points) == N) or (finished[0] == True):
            points = np.array(points)
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
            return np.transpose(points)
        if forward == True:
            curphi = nextphi(curphi)
            add_point(fullcoord(curphi, z))
            if(rad(curphi) > r2):
                forward = False
                z += 1
                add_point(fullcoord(curphi, z))
        else:
            curphi = prevphi(curphi)
            add_point(fullcoord(curphi, z))
            if(rad(curphi) < r1):
                forward = True
                z += 1
                add_point(fullcoord(curphi, z))


def grow_rw(step, size, method="line"):
    numpy = np
    t = size / 2
    if method == "standart":
        a = [(t, t, t), (t, t, t + 1), (t, t + 1, t + 1), (t, t + 1, t)]
    if method == "line":
        a = []
        for i in xrange(1, size):
            a.append((t, t, i))

        for i in xrange(size - 1, 0, -1):
            a.append((t, t - 1, i))
    if method == "linear":
        a = []
        for i in xrange(0, size + 1):
            a.append((t, t, i))
        if (len(a) % 2) != (step % 2):
            a = a[:-1]


    b = numpy.zeros((size + 1, size + 1, size + 1), int)
    for i in a:
        b[i] = 1
    for i in xrange((step - len(a)) / 2):
        print len(a)
        while True:
            t = numpy.random.randint(0, len(a))
            if t != len(a) - 1:
                c = numpy.abs(numpy.array(a[t]) - numpy.array(a[t + 1]))
                t0 = numpy.array(a[t])
                t1 = numpy.array(a[t + 1])
            else:
                c = numpy.abs(numpy.array(a[t]) - numpy.array(a[0]))
                t0 = numpy.array(a[t])
                t1 = numpy.array(a[0])
            cur_direction = numpy.argmax(c)
            while True:
                direction = numpy.random.randint(0, 3)
                if direction != cur_direction:
                    break
            if numpy.random.random() > 0.5:
                shift = 1
            else:
                shift = -1
            shiftar = numpy.array([0, 0, 0])
            shiftar[direction] = shift
            t3 = t0 + shiftar
            t4 = t1 + shiftar
            if (b[tuple(t3)] == 0) and (b[tuple(t4)] == 0) and (numpy.min(t3) >= 1) and (numpy.min(t4) >= 1) and (numpy.max(t3) < size) and (numpy.max(t4) < size):
                a.insert(t + 1, tuple(t3))
                a.insert(t + 2, tuple(t4))
                b[tuple(t3)] = 1
                b[tuple(t4)] = 1
                break
        #print a
    return numpy.array(a)



def _test():

    print "testing save/load"
    a = np.random.random((20000, 3))
    save(a, "bla", mode="txt")
    b = load("bla")
    assert abs(b.mean() - a.mean()) < 0.00001

    save(a, "bla", mode="joblib")
    b = load("bla")
    assert abs(b.mean() - a.mean()) < 0.00001

    save(a, "bla", mode="h5dict")
    b = load("bla")
    assert abs(b.mean() - a.mean()) < 0.00001

    os.remove("bla")

    print "Finished testing save/load, successfull"

def createSpiralRing(N, twist, r=0, offsetPerParticle=np.pi, offset=0):
    """
    Creates a ring of length N. Then creates a spiral
    """
    if not isInteger(N * offsetPerParticle / (2 * np.pi)):
        print N * offsetPerParticle / (2 * np.pi)
        raise ValueError("offsetPerParticle*N should be multitudes of 2*Pi")
    totalTwist = twist * N
    totalTwist = np.floor(totalTwist / (2 * np.pi)) * 2 * np.pi
    alpha = np.linspace(0, 2 * np.pi, N + 1)[:-1]
    #print alpha
    twistPerParticle = totalTwist / float(N) + offsetPerParticle
    R = float(N) / (2 * np.pi)
    twist = np.cumsum(np.ones(N, dtype=float) * twistPerParticle) + offset
    #print twist
    x0 = R + r * np.cos(twist)
    z = 0 + r * np.sin(twist)
    x = x0 * np.cos(alpha)
    y = x0 * np.sin(alpha)
    return np.array(np.array([x, y, z]).T, order="C")


def getLinkingNumber(data1, data2):
    if len(data1) == 3:
        data1 = numpy.array(data1.T)
    if len(data2) == 3:
        data2 = numpy.array(data2.T)
    if len(data1[0]) != 3: raise ValueError
    if len(data2[0]) != 3: raise ValueError

    olddata = numpy.concatenate([data1, data2], axis=0)
    olddata = numpy.array(olddata, dtype=float, order="C")
    M = len(data1)
    N = len(olddata)
    returnArray = numpy.array([0])
    from scipy import weave

    support = r"""
#include <stdio.h>
#include <stdlib.h>

float *cross(float *v1, float *v2) {
    float *v1xv2 = new float[3];
    v1xv2[0]=-v1[2]*v2[1] + v1[1]*v2[2];
    v1xv2[1]=v1[2]*v2[0] - v1[0]*v2[2];
    v1xv2[2]=-v1[1]*v2[0] + v1[0]*v2[1];
    return v1xv2;
}

float *linearCombo(float *v1, float *v2, float s1, float s2) {
    float *c = new float[3];
    c[0]=s1*v1[0]+s2*v2[0];
    c[1]=s1*v1[1]+s2*v2[1];
    c[2]=s1*v1[2]+s2*v2[2];
        return c;
}

int intersectValue(float *p1, float *v1, float *p2, float *v2) {
    int x=0;
    float *v2xp2 = cross(v2,p2), *v2xp1 = cross(v2,p1), *v2xv1 = cross(v2,v1);
    float *v1xp1 = cross(v1,p1), *v1xp2 = cross(v1,p2), *v1xv2 = cross(v1,v2);
    float t1 = (v2xp2[2]-v2xp1[2])/v2xv1[2];
    float t2 = (v1xp1[2]-v1xp2[2])/v1xv2[2];
    if(t1<0 || t1>1 || t2<0 || t2>1) {
        free(v2xp2);free(v2xp1);free(v2xv1);free(v1xp1);free(v1xp2);free(v1xv2);
        return 0;
    }
    else {
        if(v1xv2[2]>=0) x=1;
        else x=-1;
    }
    float *inter1 = linearCombo(p1,v1,1,t1), *inter2 = linearCombo(p2,v2,1,t2);
    float z1 = inter1[2];
    float z2 = inter2[2];

    free(v2xp2);free(v2xp1);free(v2xv1);free(v1xp1);free(v1xp2);free(v1xv2);free(inter1);free(inter2);
    if(z1>=z2) return x;
    else return -x;
}
    """


    code = r"""
    #line 1149 "numutils.py"
    float **data = new float*[N];
    int i,j;
    for(i=0;i<N;i++) {

        data[i] = new float[3];
        data[i][0]=olddata[3*i];
        data[i][1]=olddata[3*i+1];
        data[i][2]=olddata[3*i + 2];
    }

    int L = 0;
        for(i=0;i<M;i++) {
            for(j=M;j<N;j++) {
                float *v1, *v2;
                if(i<M-1) v1 = linearCombo(data[i+1],data[i],1,-1);
                else v1 = linearCombo(data[0],data[M-1],1,-1);

                if(j<N-1) v2 = linearCombo(data[j+1],data[j],1,-1);
                else v2 = linearCombo(data[M],data[N-1],1,-1);
                L+=intersectValue(data[i],v1,data[j],v2);
                free(v1);free(v2);
            }
        }

    returnArray[0] =  L;

"""
    M, N  #Eclipse warning removal
    weave.inline(code, ['M', 'olddata', 'N', "returnArray"], extra_compile_args=['-march=native -malign-double -O3'], support_code=support)
    return returnArray[0]





#_test()
