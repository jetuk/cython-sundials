# cython: profile=True


cimport libsundials as sun

from libc.stdlib cimport abort, malloc, free
from libc.stdint cimport uintptr_t

from cpython cimport Py_INCREF, Py_DECREF

import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

import numpy as np
cimport numpy as np
np.import_array() # initialize C API to call PyArray_SimpleNewFromData

import cython


cpdef int denseGETRF( sun.realtype[:,::1] a, long[::1] p ):
    
    cdef int m = a.shape[0]
    cdef int n = a.shape[1]
    cdef ret    
    
    ret = sun.denseGETRF( &a[0,0], m, n, &p[0])
    
    return ret
    
    
cpdef denseGETRS(sun.realtype[:,::1] a, long[::1] p ):
    
    cdef int n = a.shape[0]
    cdef sun.realtype[::1] b = np.empty( (n,) )
    
    sun.denseGETRS(&a[0,0], n, &p[0], &b[0] )
    
    return np.asarray(b)


include "Nvector.pxi"

