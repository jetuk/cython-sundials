cimport libkinsol as kinsol
from sundials cimport N_Vector


#cdef class Cvode:
#    cdef void *_cv
#    cdef int _ms
#    cdef int _it
#    cdef N_Vector y0
    
cdef class BaseKinsol:
    cdef void *_kn
    cdef int _ms
    cdef int _it
    cdef N_Vector tmpl