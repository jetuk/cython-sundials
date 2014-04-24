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


include "Nvector.pxi"

PRE_CONDITION_CONSTS = {
    'none': sun.PREC_NONE, 'left': sun.PREC_LEFT, 
    'right': sun.PREC_RIGHT, 'both': sun.PREC_BOTH
}

GS_TYPES = {'modified': sun.MODIFIED_GS, 'classical':sun.CLASSICAL_GS}

# Cython wrapper of DlsMat

cdef class pyDlsMat:
        
    def __setitem__(self, index, sun.realtype value ):
        cdef long int i = index[0]
        cdef long int j = index[1]
        
        if self._m.type == sun.SUNDIALS_DENSE:
            self._m.cols[j][i] = value
            #sun.DENSE_ELEM(self._m, index[0], index[1]) = value
    
    
        