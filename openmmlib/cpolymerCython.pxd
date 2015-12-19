cdef extern from "qcprot.h":
    cdef double CalcRMSDRotationalMatrix(double **coords1, double **coords2, const int len, double *rot, const double *weight)

