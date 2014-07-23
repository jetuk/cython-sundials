cimport libsundials as sun
cimport libcvode as cvode

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
include 'cvode_properties.pxi'



cdef class Cvode(BaseCvode):
    
    pass