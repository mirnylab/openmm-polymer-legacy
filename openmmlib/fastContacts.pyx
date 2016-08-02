import numpy as np
cimport numpy as np
cimport cython

cdef extern from "<vector>" namespace "std":
    cdef cppclass vector[T]:
        cppclass iterator:
            T operator*()
            iterator operator++()
            bint operator==(iterator)
            bint operator!=(iterator)
        vector()
        void push_back(T&)
        T& operator[](int)
        T& at(int)
        iterator begin()
        iterator end()

np.import_array() # initialize C API to call PyArray_SimpleNewFromData
cdef public api tonumpyarray(int* data, long long size) with gil:
    if not (data and size >= 0): raise ValueError
    cdef np.npy_intp dims = size
    #NOTE: it doesn't take ownership of `data`. You must free `data` yourself
    return np.PyArray_SimpleNewFromData(1, &dims, np.NPY_INT, <void*>data)

@cython.boundscheck(False)
@cython.wraparound(False)
def contactsCython(inArray, cutoff):
    inArray = np.asarray(inArray, dtype = np.float64, order = "C")
    cdef int N = len(inArray)
    cdef np.ndarray[np.double_t, ndim = 2] data = inArray
    cdef int j,i
    cdef double curdist
    cdef double cutoff2 = cutoff * cutoff  # IMPORTANT to avoid slow sqrt calculation
    cdef vector[int] contacts1
    cdef vector[int] contacts2
    for i in range(N):
        for j in range(i+1, N):
            curdist = (data[i,0] - data[j,0]) **2 +(data[i,1] - data[j,1]) **2 + (data[i,2] - data[j,2]) **2
            if curdist < cutoff2:
                contacts1.push_back(i)
                contacts2.push_back(j)
    cdef int M = len(contacts1)

    cdef np.ndarray[np.int32_t, ndim = 2] contacts = np.zeros((M,2), dtype = np.int32)
    for i in range(M):
        contacts[i,0] = contacts1[i]
        contacts[i,1] = contacts2[i]
    return contacts



