import numpy as np
import numpy.linalg
cimport numpy as np 
import sys 

cdef public double det(double ** matrix, int N):
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] myarray = np.zeros((N,N), dtype=np.double)
    cdef int i,j
    sys.stderr.write('crossings: {0}; '.format(N))
    for i in xrange(N):
        myarray[i,i] = 1.1


    for i in xrange(N):
        for j in xrange(N):
            myarray[i,j] = matrix[i+1][j+1]
    cdef double ret
    ret = numpy.linalg.slogdet(myarray)[1]
    sys.stderr.write('log determinant: {0}; '.format(ret))

    return ret
    
            

