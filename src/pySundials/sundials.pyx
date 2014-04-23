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

"""
 * -----------------------------------------------------------------
 * DENSE_COL and DENSE_ELEM
 * -----------------------------------------------------------------
 *
 * DENSE_COL(A,j) references the jth column of the M-by-N dense
 * matrix A, 0 <= j < N. The type of the expression DENSE_COL(A,j) 
 * is (realtype *). After the assignment in the usage above, col_j 
 * may be treated as an array indexed from 0 to M-1. The (i,j)-th 
 * element of A is thus referenced by col_j[i].
 *
 * DENSE_ELEM(A,i,j) references the (i,j)th element of the dense 
 * M-by-N matrix A, 0 <= i < M ; 0 <= j < N.
 *
 * -----------------------------------------------------------------
"""
# These aren't really needed when using the ndarray versions of N_Vector

#define DENSE_COL(A,j) ((A->cols)[j])
#define DENSE_ELEM(A,i,j) ((A->cols)[j][i])




"""
 * -----------------------------------------------------------------
 * BAND_COL, BAND_COL_ELEM, and BAND_ELEM
 * -----------------------------------------------------------------
 *  
 * BAND_COL(A,j) references the diagonal element of the jth column 
 * of the N by N band matrix A, 0 <= j <= N-1. The type of the 
 * expression BAND_COL(A,j) is realtype *. The pointer returned by 
 * the call BAND_COL(A,j) can be treated as an array which is 
 * indexed from -(A->mu) to (A->ml).
 * 
 * BAND_COL_ELEM references the (i,j)th entry of the band matrix A 
 * when used in conjunction with BAND_COL. The index (i,j) should 
 * satisfy j-(A->mu) <= i <= j+(A->ml).
 *
 * BAND_ELEM(A,i,j) references the (i,j)th element of the M-by-N 
 * band matrix A, where 0 <= i,j <= N-1. The location (i,j) should 
 * further satisfy j-(A->mu) <= i <= j+(A->ml). 
 *
 * -----------------------------------------------------------------
"""
 
#define BAND_COL(A,j) (((A->cols)[j])+(A->s_mu))
#define BAND_COL_ELEM(col_j,i,j) (col_j[(i)-(j)])
#define BAND_ELEM(A,i,j) ((A->cols)[j][(i)-(j)+(A->s_mu)])


# Cython wrapper of DlsMat

cdef class pyDlsMat:
        
    def __setitem__(self, index, sun.realtype value ):
        cdef long int i = index[0]
        cdef long int j = index[1]

        if self._m.type == sun.SUNDIALS_DENSE:
            self._m.cols[j][i] = value
            #sun.DENSE_ELEM(self._m, index[0], index[1]) = value
        elif self._m.type == sun.SUNDIALS_BAND:
            assert j-self._m.mu <= i and i <= j+self._m.ml
            self._m.cols[j][i-j+self._m.s_mu] = value
    
    
        
