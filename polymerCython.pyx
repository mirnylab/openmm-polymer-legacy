import numpy as np 
cimport numpy as np
cimport cython   


from cpolymerCython cimport *

ctypedef fused my_type:
    cython.int
    cython.float
    cython.char
    cython.long
    cython.short
    cython.double
    cython.complex
        

    
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def fastTxtSave(np.ndarray data,filename):
    
    cdef np.ndarray[np.float64_t,ndim = 2] _data = np.array(data, dtype=float)
    cdef int i 
    cdef np.ndarray[np.float64_t,ndim = 1] particle
    
    lines = []
    lines.append(str(len(data)) + "\n")
    
    for i in range(len(data)):
        particle = _data[i]
        lines.append("%lf %lf %lf\n" % (particle[0], particle[1], particle[2]))
    
    if type(filename) == str:
        with open(filename, 'w') as myfile:
            myfile.writelines(lines)

    elif hasattr(filename, "writelines"):
        filename.writelines(lines)
    else:
        return lines


def fastMSD(data1, data2):
    """
    a faster implementation of mean square displacement. 
    """

    mean1 = np.mean(data1, axis=0)
    data1 -= mean1[None, :]
    mean2 = np.mean(data2, axis=0)
    data2 -= mean2[None, :]

    cdef np.ndarray[np.double_t, ndim = 2, mode="c"] data1_cpp = np.asanyarray(data1.T, dtype=np.double, order="C")
    cdef np.ndarray[np.double_t, ndim = 2, mode="c"] data2_cpp = np.asanyarray(data2.T, dtype=np.double, order="C")
    
    cdef int N = len(data1)
    cdef double* weight = NULL
    cdef double* data1_p[3]
    cdef double* data2_p[3]
    cdef double result    
    for i in range(3):
        data1_p[i] = &data1_cpp[i,0]
        data2_p[i] = &data2_cpp[i,0]
    
    
    cdef np.ndarray[np.double_t, ndim = 1] rot = np.zeros(9, dtype=np.double)
    
     
    result =  CalcRMSDRotationalMatrix(data1_p, data2_p, N,  &rot[0], NULL) #clib    
    return np.double(result)



