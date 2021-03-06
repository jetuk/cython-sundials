cimport libsundials as sun
#cimport libcvode as cvode
#cimport libkinsol as kinsol
from cpython cimport bool
cimport numpy as np

cdef class N_Vector:
    cdef sun.N_Vector _v
    cdef sun.N_Vector_Ops _ops
    cdef long LengthRealType
    cdef long LengthIntType
    

    
cdef class NvectorNdarrayFloat64(N_Vector):
    cdef public np.ndarray data
    cdef object shp    


cdef class NvectorMemoryViewDouble5D(N_Vector):
    cdef public double[:, :, :, :, :] data
    cdef np.ndarray _data
    cdef object shp
    cdef bint debug
    
cdef class VariableNvector(N_Vector):
    cdef public dict variables
    
cdef class Variable1D(object):
    cdef public str name
    cdef public dict saved_values
    cdef list _copied_attr_names
    cdef double[:] value
    cdef dict _attrs        

cdef class Variable1DNvector(N_Vector):
    cdef public list variables
    
    
    
    
cdef class pyDlsMat:
    cdef sun.DlsMat _m
    