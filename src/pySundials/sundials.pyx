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

