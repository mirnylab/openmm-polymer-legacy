import numpy
from scipy import weave
import os
import os.path
from tempfile import NamedTemporaryFile
from time import sleep

folderName = os.path.split(__file__)[0]
reduceKnotFilename = os.path.join(folderName, "Reduce_knot20")


def findSimplifiedPolymer(data):
    """a weave.inline wrapper for polymer simplification code
    Calculates a simplified topologically equivalent polymer ring"""

    if len(data) != 3:
        data = numpy.transpose(data)
    if len(data) != 3:
        raise ValueError("Wrong dimensions of data")
    datax = numpy.array(data[0], float, order="C")
    datay = numpy.array(data[1], float, order="C")
    dataz = numpy.array(data[2], float, order="C")
    N = len(datax)
    ret = numpy.array([1])
    datax, datay, dataz, N  # eclipse warning removal
    code = r"""
    #line 290 "binary_search.py"
    int M = 0;
    int sum = 0;
    int t=0,s=0,k=0;
    int turn=0;
    bool breakflag;
    double maxdist=0;
    int a;
    position=vector<point>(N);
    newposition=vector<point>(N);

    for (i=0;i<N;i++)
    {
    position[i].x = datax[i] +  0.00001*(rand()%1000);
    position[i].y = datay[i] +0.00001*(rand()%1000);
    position[i].z  = dataz[i] +  0.00001*(rand()%1000);
    }
    todelete = vector <int> (N);
    for (i=0;i<N;i++) todelete[i] == -2;
    while (true)
        {
        maxdist = 0;
        for (i=0;i<N-1;i++)
        {
        if (dist(i,i+1) > maxdist) {maxdist = dist(i,i+1);}
        }
        //printf ("maxdist = %lf\n",maxdist);
        turn++;
        M=0;
        for (i=0;i<N;i++) todelete[i] = -2;
        for (int j=1;j<N-1;j++)  //going over all elements trying to delete
            {

            breakflag = false; //by default we delete thing
            for (k=0;k<N-1;k++)  //going over all triangles to check
                {
                double dd = dist(j,k);
                if (dd  < 3 * maxdist)
                {

                if (k < j-2 || k > j+1)
                    {
                    sum = intersect(position[j-1],position[j],position[
                        j+1],position[k],position[k+1]);
                    if (sum!=0)
                        {
                        //printf("intersection at %d,%d\n",j,k);
                        breakflag = true; //keeping thing
                        break;
                        }
                    }
                }
                else k+= abs((int)((float)dd/(float)maxdist )- 3);
                }
            if (breakflag ==false)
            {
            todelete[M++] = j;
            position[j] = (position[j-1] + position[j+1])* 0.5;
            //printf("%d will be deleted at %d\n",j,k);
            j++;
            //break;
            }
            }
        t = 0;//pointer for todelete
        s = 0;//pointer for newposition
        if (M==0)
            {
            break;
            }
        for (int j=0;j<N;j++)
            {
            if (todelete[t] == j)
                {
                t++;
                continue;
                }
            else
                {
                newposition[s++] = position[j];
                }
            }
        N = s;
        M = 0;
        t = 0;
        position = newposition;
        }
    ret[0] = N;

    for (i=0;i<N;i++)
    {
    datax[i]  = position[i].x;
    datay[i]  = position[i].y;
    dataz[i]  = position[i].z;
    }
    """
    support = r"""
#line 400 "binary_search.py"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <ctime>
#include <omp.h>
#include <stdio.h>
using namespace std;
struct point{
    double x,y,z;
    point operator + (const point &p) const {
        return (point) {x+p.x, y+p.y, z+p.z};
    }
    point operator - (const point &p) const {
        return (point) {x-p.x, y-p.y, z-p.z};
    }
/* cross product */
    point operator * (const point &p) const {
        return (point) {y*p.z - z*p.y,
                        z*p.x - x*p.z,
                        x*p.y - y*p.x};
    }
    point operator * (const double &d) const {
        return (point) {d*x, d*y, d*z};
    }

    point operator / (const double &d) const {
        return (point) {x/d, y/d, z/d};
    }
};

vector <point> position;
vector <point> newposition;
vector <int> todelete;
int N;
int i;
double dist(int i,int j);
double dotProduct(point a,point b);
int intersect(point t1,point t2,point t3,point r1,point r2);

inline double sqr(double x){
    return x*x;
}
inline double dist(int i,int j){
return sqrt(dotProduct((position[i]-position[j]),(position[i]-position[j])));
}

inline double dist(point a,point b){
    return sqr(a.x-b.x)+sqr(a.y-b.y)+sqr(a.z-b.z);
}

inline double dotProduct(point a,point b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

int intersect(point t1,point t2,point t3,point r1,point r2)
{
point A,B,C,D,n;
int r;
double det,t,u,v,c1,d1,d2,d3;
B = t2 - t1;
C = t3 - t1;
D = r2 - t1;
A = r2 - r1;

d1 = (B.y*C.z-C.y*B.z);
d2 = (B.x*C.z-B.z*C.x);
d3 = (B.x*C.y-C.x*B.y);
det = A.x*d1-A.y*d2+A.z*d3;
if (det == 0) return 0;
if (det >0){
t = D.x*d1-D.y*d2+D.z*d3;
if (t<0 || t>det) return 0;
u = A.x*(D.y*C.z-C.y*D.z)-A.y*(D.x*C.z-D.z*C.x)+A.z*(D.x*C.y-C.x*D.y);
if (u<0 || u>det) return 0;
v = A.x*(B.y*D.z-D.y*B.z)-A.y*(B.x*D.z-B.z*D.x)+A.z*(B.x*D.y-D.x*B.y);
if (v<0 || v>det || (u+v)>det) return 0;
//printf("\n%lf,%lf,%lf, ",t/det,u/det,v/det);
n = B*C;
c1 = dotProduct(r1-t1,n);
if (c1>0) return 1;
else return -1;
}
else{
t = D.x*d1-D.y*d2+D.z*d3;
if (t>0 || t<det) return 0;
u = A.x*(D.y*C.z-C.y*D.z)-A.y*(D.x*C.z-D.z*C.x)+A.z*(D.x*C.y-C.x*D.y);
if (u>0 || u<det) return 0;
v = A.x*(B.y*D.z-D.y*B.z)-A.y*(B.x*D.z-B.z*D.x)+A.z*(B.x*D.y-D.x*B.y);
if (v>0 || v<det || (u+v)<det) return 0;
//printf("\n%lf,%lf,%lf, ",t/det,u/det,v/det);
n = B*C;
c1 = dotProduct(r1-t1,n);
if (c1>0) return 1;
else return -1;
}
}
//DNA conformation
"""
    weave.inline(code, ['datax', 'datay', 'dataz', 'N', 'ret'],
                 extra_compile_args=['-malign-double'], support_code=support)

    data = numpy.array([datax, datay, dataz]).T

    return data[:ret[0]]


def getKnotNumber(data, evalAt= -1.1):
    data = numpy.array(data)
    if len(data) == 3:
        data = data.T

    with  NamedTemporaryFile() as newfile:
        newfile.write("t=0\n\n%d\n" % len(data))
        for j, i in enumerate(data):
            newfile.write("%d %lf %lf %lf\n" % tuple([j + 1] + list(i)))

        name = newfile.name
        newfile.flush()

        os.system("{0} {1} -p {4}  > {2}_{3}".format(reduceKnotFilename, name,
                                     name, "_output", evalAt))
        lines = open("%s_%s" % (name, "_output")).readlines()
        os.remove("%s_%s" % (name, "_output"))
        return lines


def expandPolymerRing(data, mode="auto"):
    """
    Expands polymer ring or chain using OpenMM.

    Parameters
    ----------
    data : Nx3 or 3xN array of coordinates
        Input coordinates of the polymer
    mode : str, optional
        "ring", or "chain", default - autodetect
    """

    from openmmlib import Simulation
    sim = Simulation(
        timestep=100, thermostat=0.002, velocityReinitialize=True)
    sim.setup(platform="opencl")
    sim.load(data)
    sim.randomizeData()
    if mode == "auto":
        if sim.dist(0, sim.N - 1) < 2:
            mode = "ring"
        else:
            mode = "chain"
    sim.setLayout(mode=mode)
    sim.addHarmonicPolymerBonds(wiggleDist=0.1)
    sim.addGrosbergRepulsiveForce(trunc=40)
    sim.addGrosbergStiffness(k=2)
    #sim.energyMinimization(stepsPerIteration=50)
    sim.doBlock(10)
    for _ in xrange(10):
        sim.doBlock(1000)
    data = sim.getData()
    del sim
    sleep(0.5)
    return data



def analyzeKnot(data, useOpenmm=False, evalAt= -1.1):

    data = numpy.asarray(data)
    if len(data) == 3:
        data = data.T

    t = findSimplifiedPolymer(data)
    if useOpenmm == True:
        if len(t) > 200:
            data = expandPolymerRing(data)
            t = findSimplifiedPolymer(data)


    #t = data
    print "simplified to: %d monomers" % len(t)
    number = getKnotNumber(t, evalAt=evalAt)
    print number
    num = float(number[0].split()[1])
    return num
