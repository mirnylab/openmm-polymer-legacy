"""
This file contains a bunch of method to work on contact maps of a Hi-C data. 
It uses a lot of methods from mirnylab repository. 

Save/load functions
-------------------

Use :py:func:`load` to load any polymer. 

Use :py:func:`save` to save an XYZ. Otherwise you can just use joblib.dump({"data":data},filename)

Other functions are about to be deprecated


Find contacts of a conformation
-------------------------------

To be filled in later


Find average contact maps
-------------------------

bla

"""

from mirnylab import numutils  
from array import array 
import numpy 
from scipy import weave
from math import sqrt
from mirnylab.systemutils import fmapred,fmap 
import joblib



def load(filename,center = False):
    """loads joblib or xyz polymer file.      
    
    Parameters
    ----------
    Filename : str
        Filename to load
    center : bool, optional
        Shift center of mass to zero, default = True. 
    """
    try: int(open(filename).readline())
    except: 
        data = joblib.load(filename)["data"]
        if center == True: data -= numpy.mean(data,axis = 0)[None,:]
        return data
    
    f = open(filename, 'r')
    lines = f.readlines()
    N = int(lines[0])            
    datax = array('f', [0.] * N )
    datay = array('f', [0.] * N )
    dataz = array('f', [0.] * N )
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
    datax -= datax.mean()
    datay -= datay.mean()
    dataz -= dataz.mean()
    return numpy.array([datax,datay,dataz])



def Cload(filename,center = False):
    """fast polymer loader using weave.inline
    
    ..warning:: About to be deprecated 
    """
    f = open(filename, 'r')
    line = f.readline()
    N = int(line)
    ret = numpy.zeros((3,N),order = "C",dtype = numpy.double)
    code = """
    #line 85 "binary_search.py"
    using namespace std;
    FILE *in;
    const char * myfile = filename.c_str(); 
    in = fopen(myfile,"r");
    int dummy;
    dummy = fscanf(in,"%d",&dummy);    
    for (int i = 0; i < N; i++)
    {
    dummy = fscanf(in,"%lf %lf %lf ",&ret[i],&ret[i + N],&ret[i + 2 * N]);
    if (dummy < 3){printf("Error in the file!!!");throw -1;}
    } 
    """
    support = """
    #include <math.h>  
    """
    weave.inline(code, ['filename', 'N' , 'ret' ], extra_compile_args=['-march=native -malign-double'],support_code =support )
    if center == True: ret -=  numpy.mean(ret,axis = 1)[:,None]
    return ret
    


def intload(filename,center = "N/A"):
    """
    .. warning:: About to be deprecated 
    """
    "dedicated weave.inline for loading int polymers"
    f = open(filename, 'r')
    line = f.readline()
    f.close()
    N = int(line)
    ret = numpy.zeros((3,N),order = "C",dtype = numpy.int)
    code = """
    #line 85 "binary_search.py"
    using namespace std;
    FILE *in;
    const char * myfile = filename.c_str(); 
    in = fopen(myfile,"r");
    int dummy;
    dummy = fscanf(in,"%d",&dummy);    
    for (int i = 0; i < N; i++)
    {
    dummy = fscanf(in,"%ld %ld %ld ",&ret[i],&ret[i + N],&ret[i + 2 * N]);
    if (dummy < 3){printf("Error in the file!!!");throw -1;}        
    } 
    fclose(in);
    """
    support = """
    #include <math.h>  
    """
    weave.inline(code, ['filename', 'N' , 'ret' ], extra_compile_args=['-march=native -malign-double'],support_code =support )
    
    return ret

def save(filename, data,doint=False):
    """
    Old function to save an xyz conformation. 
    .. warning:: About to be deprecated. 
    
    """
    f = open(filename, 'w')
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


def rad2(data):
    "returns Rg(N^(2/3)"
    
    if len(data) != 3: data = numpy.transpose(data)
    if len(data) != 3: raise ValueError("Wrong dimensions of data")    
    def give_radius_scaling(data):
        N = len(data[0])
        target = int(N**(2/3.))    
        coms = numpy.cumsum(data,1)
        coms2 = numpy.cumsum(data**2,1)
        def radius_gyration(len2):
            coms2d = (-coms2[:,:-len2]+coms2[:,len2:])/len2
            comsd = ((coms[:,:-len2]-coms[:,len2:])/len2)**2
            diffs = coms2d - comsd
            sums = numpy.sqrt(numpy.sum(diffs,0))
            return numpy.mean(sums)
        return radius_gyration(target)     
    return give_radius_scaling(data)



def giveIntContacts(data):    
    """give all contacts of a polymer on a cubic lattice
    Intersections are not counted as contacts. Sorry :( 
    
    Parameters
    ----------
    data : Nx3 or 3xN array of ints
    """

    data = numpy.asarray(data,dtype  = int)
    
    if len(data.shape) != 2: raise ValueError("Wrong dimensions of data")
    if 3 not in data.shape: raise ValueError("Wrong size of data: %s,%s" % data.shape)
    if data.shape[0] == 3: data = data.T
    
    data -= numpy.min(data,axis = 0)[None,:]    
    
    M = numpy.max(data)+1
    if M>1500: raise ValueError("Polymer is to big, can't create bounding box!")
    
    N = len(data)
    tocheck = numpy.zeros(M*M*M,dtype = numpy.int32) - 1
    tocheck[data[:,0] + data[:,1] * M + data[:,2] * M * M] = numpy.arange(N,dtype = numpy.int32)
    tocheck.shape = (M,M,M)
    contacts1 = numpy.concatenate([tocheck[1:,:,:].ravel(),tocheck[:,1:,:].ravel(),tocheck[:,:,1:].ravel()]) 
    contacts2 = numpy.concatenate([tocheck[:-1,:,:].ravel(),tocheck[:,:-1,:].ravel(),tocheck[:,:,:-1].ravel()])        
    mask = (contacts1 != -1) * (contacts2 != -1)
    contacts1 = contacts1[mask]
    contacts2 = contacts2[mask]
    contacts3 = numpy.minimum(contacts1,contacts2)
    contacts4 = numpy.maximum(contacts1,contacts2)    
    return numpy.concatenate([contacts3[:,None],contacts4[:,None]],1)
    



 
def giveContacts(data,cutoff=1.7,maxContacts = 40 ):
    """Returns contacts of a single polymer with a given cutoff
    
    .. warning:: Use this only to find contacts of a single polymer chain with distance between monomers of 1. 
    Multiple chains will lead to silent bugs. 
 
    Parameters
    ----------
    data : Nx3 or 3xN array
        Polymer configuration. One chaon only. 
    cutoff : float , optional
        Cutoff distance that defines contact
    maxContacts : int 
        Maximum number of contacts per monomer. 
        If total number of contacts exceeds maxContacts*N, program becomes uncontrollable. 
        
    Returns
    -------
    
    k by 2 array of contacts. Each row corresponds to a contact. 
    """
    data = numpy.asarray(data)     
    if len(data.shape) != 2: raise ValueError("Wrong dimensions of data")
    if 3 not in data.shape: raise ValueError("Wrong size of data: %s,%s" % data.shape)
    if data.shape[0] == 3: data = data.T
    data = numpy.asarray(data,float,order = "C")
    
    ####---------------------------------------------------- Finish it!     
#    dists2 = numpy.sum(numpy.diff(data,axis = 0)**2,axis = 1)
#    if dists2.max() > 3:
#        raise RuntimeError("Maximum distance for a ")
    #-------------------------------------------------------------
    N = len(data)     
    points = numpy.zeros((maxContacts*N,2),int,order = "C")
    N #Eclipse warning remover 
    code = r"""
    #line 196 "binary_search.py"
        using namespace std; 
        int counter  = 0;
          double  d;
        for (int i=0;i<N;i++)
            {
            for (int j=i+1;j<N;j++)
                {
                d = dist(data,N,i,j);
                
                if (d<= cutoff) 
                    {
                    int t = j - i;                     
                    points[2*counter] = i;
                    points[2*counter + 1 ] = j; 
                    counter ++ ;
                    //cout << i << " "   << j  << " "  << d <<  " " << counter << endl;   
                
                    }
                else if (d>4)
                    {
                    j +=  d-4; 
                    }
            }
            } 
       """
    support = """
        #include <math.h>  
        double  dist(double* data, int N, int i,int j)
        {
        return sqrt(pow((data[3 * i]-data[3 * j]),2.) + pow((data[3 * i + 1] - data[3 * j + 1]),2.) + pow((data[3 * i + 2] - data[3 * j + 2]),2.));
        }
       """
    weave.inline(code, ['data', 'N' , 'cutoff' , 'points'], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
    k = numpy.max(numpy.nonzero(points[:,1]))     
    return points[:k+1,:]

 
def giveContactsAny(data,cutoff=1.4,maxContacts = 100):
    """Returns contacts of any sets of molecules with a given cutoff. 
    Is slower than give_contacts, but can tolerate multiple chains. 
    
    Parameters
    ----------
    data : Nx3 or 3xN array
        Polymer configuration. One chaon only. 
    cutoff : float , optional
        Cutoff distance that defines contact
    maxContacts : int 
        Maximum number of contacts per monomer. 
        If total number of contacts exceeds maxContacts*N, program becomes uncontrollable. 
        
    Returns
    -------
    
    k by 2 array of contacts. Each row corresponds to a contact. 
    """    
    data = numpy.asarray(data)
    if len(data.shape) != 2: raise ValueError("Wrong dimensions of data")
    if 3 not in data.shape: raise ValueError("Wrong size of data: %s,%s" % data.shape)
    if data.shape[0] == 3: data = data.T
    data = numpy.asarray(data,float,order = "C")
    N = len(data) 
    points = numpy.zeros((maxContacts*N,2),int,order = "C")
    N #Eclipse warning remover 
    
    code = """
    #line 50 "binary_search.py"
        using namespace std; 
    
        
        int counter  = 0;
         
        double  d;
        for (int i=0;i<N;i++)
            {
            for (int j=i+1;j<N;j++)
                {
                 
                if (dist(data,N,i,j) <= CUTOFF2) 
                    {
                    int t = j - i;                     
                    points[2*counter] = i;
                    points[2*counter + 1 ] = j; 
                    counter ++ ;
                    //cout << i << " "   << j  << " "  << d <<  " " << counter << endl;                   
                    }
            }
            }
       """
    support = """
        #include <math.h>
        #define CUTOFF2  %lf
        double sq(double x) {return x*x;} 
        double  dist(double* data, int N, int i,int j)
        {
        return (pow((data[3 * i]-data[3 * j]),2.) + pow((data[3 * i + 1] - data[3 * j + 1]),2.) + pow((data[3 * i + 2] - data[3 * j + 2]),2.)); 
        }
       """ % (cutoff*cutoff) 
    #
    weave.inline(code, ['data', 'N' , 'cutoff' , 'points'], extra_compile_args=['-march=native -malign-double'],support_code =support )
    k = numpy.max(numpy.nonzero(points[:,1]))
    return points[:k+1,:]

def giveDistanceMap(data,size = 1000):
    """returns  a distance map of a polymer with a given size"""    
    toret = numpy.zeros((size,size),float,order = "C")
    
    data = numpy.asarray(data)
    if len(data.shape) != 2: raise ValueError("Wrong dimensions of data")
    if 3 not in data.shape: raise ValueError("Wrong size of data: %s,%s" % data.shape)
    if data.shape[0] == 3: data = data.T
    data = numpy.asarray(data,float,order = "C")
    
    N = len(data)         
    N #Eclipse warning remover 
    code = """
    #line 50 "binary_search.py"
        using namespace std; 
        int counter  = 0;
          double  d;
        for (int i=0;i<size;i++)
            {
            for (int j=i;j<size;j++)
                {
                d = dist(data,N,N*i/size,N*j/size);
                toret[size*i+j] = d;                 
                toret[size*j+i] = d;
            }
            }
       """ 
    support = """
        #include <math.h>  
        double  dist(double* data, int N, int i,int j)
        {
        return sqrt(pow((data[3 * i]-data[3 * j]),2.) + pow((data[3 * i + 1] - data[3 * j + 1]),2.) + pow((data[3 * i + 2] - data[3 * j + 2]),2.)); 
        }
       """
    weave.inline(code, ['data', 'N','size', 'toret'], extra_compile_args=['-march=native -malign-double'],support_code =support )
    return toret



def findSimplifiedPolymer(data):
    """a weave.inline wrapper for polymer simplification code
    Calculates a simplified topologically equivalent polymer ring"""
    
    if len(data) != 3: data = numpy.transpose(data)
    if len(data) != 3: raise ValueError("Wrong dimensions of data")        
    datax = numpy.array(data[0],float, order = "C")
    datay = numpy.array(data[1], float,order = "C")
    dataz = numpy.array(data[2], float,order = "C")
    N = len(datax)
    ret = numpy.array([1])    
    datax,datay,dataz,N ##eclipse warning removal 
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
                    sum = intersect(position[j-1],position[j],position[j+1],position[k],position[k+1]);
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
    weave.inline(code, ['datax', 'datay','dataz', 'N','ret'], extra_compile_args=['-malign-double'],support_code =support )
    return ret[0]


def rescalePoints(points, res = 100):
    "converts array of contacts to the reduced resolution contact map"   
    a = numpy.histogram2d(points[:,0],points[:,1],res)[0]
    a = a + numpy.transpose(a)
    return a 

 
def rescaledMap(data,res,cutoff = 1.4 ):
    "calculates a rescaled contact map of a structure"    
    t = giveContacts(data,cutoff)
    return rescalePoints(t,res) 
    
def pureMap(data,cutoff=1.4,contactMap = None):
    """calculates an all-by-all contact map of a single polymer chain.
    Doesn't work for multi-chain polymers!!!  
    If contact map is supplied, it just updates it 
    
    Parameters
    ----------
    data : Nx3 or 3xN array
        polymer conformation
    cutoff : float
        cutoff for contacts
    contactMap : NxN array, optional 
        contact map to update, if averaging is used 
    """
    
    if len(data) != 3: data = numpy.transpose(data)
    if len(data) != 3: raise ValueError("Wrong dimensions of data")        
    t = giveContacts(data,cutoff)
    N = len(data[0])
    if contactMap == None: contactMap = numpy.zeros((N,N),int)
    contactMap[t[:,0],t[:,1]] += 1 
    contactMap[t[:,1],t[:,0]] += 1    
    return contactMap  
        

def observedOverExpected(matrix):
    "Calculates observedOverExpected of any contact map"
    data = numpy.array(matrix, order = "C")
    N = data.shape[0]
    bins = numutils.logbins(1,N,1.2)
    bins = [(0,1)] + [(bins[i],bins[i+1]) for i in xrange(len(bins)-1)]
    bins = numpy.array(bins,order = "C")
    M = len(bins)
    M #Eclipse warning remover
    code = r"""
    #line 50 "binary_search.py"
    using namespace std;
    for (int bin = 0; bin < M; bin++)
    {
        int start = bins[2 * bin];
        int end = bins[2 * bin + 1];
        
        double ss = 0 ;
        int count   = 0 ;  
        for (int offset = start; offset < end; offset ++)
        {
            for (int j = 0; j < N - offset; j++)
            {
                ss += data[(offset + j) * N + j];
                count += 1;                            
            }
        }
        double meanss = ss / count;
        printf("%lf\n",meanss); 
        for (int offset = start; offset < end; offset ++)
        {
            for (int j = 0; j < N - offset; j++)
            {
                data[(offset + j) * N + j] /= meanss;                                             
                if (offset > 0) {data[(offset + j)  + j*N] /= meanss;}
            }
        }
    }
    
    """
    support = """
    #include <math.h>  
    """
    weave.inline(code, ['data', 'bins' , 'N' ,'M'], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
    return data
    



def cool_trunk(data):
    "somehow trunkates the globule so that the ends are on the surface"
    CUTOFF = 0.7
    CUTOFF2 = 3    
    datax, datay, dataz = data[0], data[1], data[2]
    N = len(data[0])
    found1, found2 = 0, 0
    def dist(i):
        return sqrt(datax[i] ** 2 + datay[i] ** 2 + dataz[i] ** 2)
    def sqdist(i, j):
        return sqrt((datax[i] - datax[j]) ** 2 + (datay[i] - datay[j]) ** 2 + (dataz[i] - dataz[j]) ** 2)
    def sqdist2(i, j, scale):
        return sqrt((scale * datax[i] - datax[j]) ** 2 + (scale * datay[i] - datay[j]) ** 2 + (scale * dataz[i] - dataz[j]) ** 2)
    breakflag = 0
    escapeflag = 0
    for i in xrange(N):
        exitflag = 0
        pace = 0.25 / dist(i)
        for j in numpy.arange(1, 10, pace):
            #print j
            escapeflag = 1
            for k in xrange(N):
                if k == i: continue
                escapeflag, breakflag
                
                if 0.001 < sqdist2(i, k, j) < CUTOFF:
                    breakflag = 1
                    print "breaking at", k
                    break
                if 0.001 < sqdist2(i, k, j) < CUTOFF2:                    
                    escapeflag = 0
            
            if breakflag == 1:
                breakflag = 0
                print i, dist(i), j
                break
            if escapeflag == 1:
                print i, dist(i)
                found1 = i
                exitflag
                exitflag = 1
                break
        if exitflag == 1: break
    for i in xrange(N - 1, 0, - 1):
        exitflag = 0
        pace = 0.25 / dist(i)
        for j in numpy.arange(1, 10, pace):
            #print j
            escapeflag = 1
            for k in xrange(N - 1, 0, - 1):
                if k == i: continue
                escapeflag, breakflag
                
                if 0.001 < sqdist2(i, k, j) < CUTOFF:
                    breakflag = 1
                    print "breaking at", k
                    break
                if 0.001 < sqdist2(i, k, j) < CUTOFF2:                    
                    escapeflag = 0
            
            if breakflag == 1:
                breakflag = 0
                print i, dist(i), j
                break
            if escapeflag == 1:
                print i, dist(i)
                found2 = i
                exitflag
                exitflag = 1
                break
        if exitflag == 1: break
        
    f1 = found1
    f2 = found2
    datax[f1 - 1], datay[f1 - 1], dataz[f1 - 1] = datax[f1] * 2, datay[f1] * 2, dataz[f1] * 2    
    datax[f2 + 1], datay[f2 + 1], dataz[f2 + 1] = datax[f2] * 2, datay[f2] * 2, dataz[f2] * 2
    a = (datax[f1 - 1], datay[f1 - 1], dataz[f1 - 1])
    b = (datax[f2 + 1], datay[f2 + 1], dataz[f2 + 1])
    c = (a[0] + b[0], a[1] + b[1], a[2] + b[2])
    absc = sqrt(c[0] ** 2 + c[1] ** 2 + c[2] ** 2)
    c = (c[0] * 40 / absc, c[1] * 40 / absc, c[2] * 40 / absc)
    datax[f1 - 2], datay[f1 - 2], dataz[f1 - 2] = c
    datax[f2 + 2], datay[f2 + 2], dataz[f2 + 2] = c        
    return [datax[f1 - 2:f2 + 3], datay[f1 - 2:f2 + 3], dataz[f1 - 2:f2 + 3]]


def averageContactMap(filenames, resolution = 500 ,  cutoff = 1.7, usePureMap = False, n = 4 , loadFunction = load, exceptionsToIgnore = []):
    """
    Returns an average contact map of a set of conformations. 
    Non-existing files are ignored if exceptionsToIgnore is set to IOError. 
    Can use both rescaled or pure contact map.
    example:\n
    filenames = ["/home/magus/evo/GO41_interaction_test/different_strength/expanded%d.dat" % i for i in xrange(100)] \n            
    mat_img(numpy.log(  averageContactMap(filenames,500,usePureMap = True)   +1))  \n
     
    Parameters
    ----------
    filenames : list of strings
        Filenames to average map over
    resolution : int, optional
        Resolution for rescaled map, not needed for pure map
    cutoff : float, optional
        Cutoff to calculate contacts
    usePureMap : bool, optional
        Calculate a pure (NxN) contact map. This may be slow for N>5000 
    n : int, optional 
        Number of threads to use. By default 4 to minimize RAM consumption with pure maps.
    exceptionsToIgnore : list of Exceptions
        List of exceptions to ignore when finding the contact map. 
        Put IOError there if you want it to ignore missing files. 
    
    """        
    
    Nbase = resolution
    if usePureMap == False:
        
        def action(i):   #Fetch rescale map from filename.  
            print i
            try:data = loadFunction(i)    #load filename 
            except exceptionsToIgnore:        # if file faled to load, return empty map
                print "file not found"
                return numpy.zeros((Nbase,Nbase),"float")             
            return rescaledMap(data,Nbase,cutoff)   #return rescaled map 
        
        return fmapred(action, filenames ,n=n,exceptionList= exceptionsToIgnore)   #using simple fork-map-reduce with out exception list
    else:
        """
        Now we actually need to modify our contact map by adding contacts from each new file to the contact map. 
        We do it this way because our contact map is huge (maybe a gigabyte!), 
        so we can't just add many gigabyte-sized arrays together. 
        Instead of this each worker creates an empty "average contact map", 
        and then loads files one by one and adds contacts from each file to a contact map. 
        Maps from different workers are then added together manually. 
        """
        n = min(n,len(filenames))   
        subvalues = [filenames[i::n] for i in xrange(n)]
        def myaction(values):   #our worker receives some filenames
            mysum = None   #contact map. 
            for i in values:
                try:
                    data = loadFunction(i)
                    print i
                except exceptionsToIgnore:
                    print "file not found", i
                    continue
                if data.shape[0] == 3: data = data.T
                if mysum == None: #if it's the first filename, 
                    N = len(data)#find the size of the data,
                    mysum = numpy.zeros((N,N),dtype = int) #create an empty map
                pureMap(data, cutoff, mysum)   #and pass it to pureMap to fill it in 
            return mysum

        blocks  =  fmap(myaction,subvalues,exceptionList= exceptionsToIgnore)
        blocks = [i for i in blocks if i != None]
        a = blocks[0]
        for i in blocks[1:]:
            a = a + i
        return a
    
        
