cimport libsundials as sun
cimport libcvodes as cvode

from sundials cimport N_Vector, pyDlsMat
from sundials import PRE_CONDITION_CONSTS, GS_TYPES

from libc.stdlib cimport abort, malloc, free
from cpython cimport Py_INCREF, Py_DECREF

import cython

import numpy as np
cimport numpy as np
np.import_array() # initialize C API to call PyArray_SimpleNewFromData



class CVodeError(Exception):
    pass

include 'denseGET.pxi'
include 'cvodes_properties.pxi'


cdef class Cvodes(BaseCvodes):
    
    def __init__(self, *args, **kwds):        
        self._nrtfn = kwds.pop('nrtfn', 0)
        
