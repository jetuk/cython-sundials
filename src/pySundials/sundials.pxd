cimport libsundials as sun
#cimport libcvode as cvode
#cimport libkinsol as kinsol
from cpython cimport bool
cimport numpy as np
from mpi4py cimport MPI

cdef class N_Vector:
    cdef sun.N_Vector _v
    cdef sun.N_Vector_Ops _ops
    cdef long LengthRealType
    cdef long LengthIntType
    
cdef class NvectorNdarrayFloat64(N_Vector):
    cdef public np.ndarray data
    cdef object shp    

cdef class NvectorMPINdarrayFloat64(N_Vector):
    cdef public np.ndarray data
    cdef object local_shp    
    cdef object global_shp    
    cdef MPI.Comm comm


cdef class NvectorMemoryViewDouble5D(N_Vector):
    cdef public double[:, :, :, :, :] data
    cdef np.ndarray _data
    cdef object shp
    cdef bint debug
    
cdef class pyDlsMat:
    cdef sun.DlsMat _m
    