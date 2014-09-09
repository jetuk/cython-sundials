# -*- coding: utf-8 -*-


import numpy as np
from pySundials.cvode import denseGETRF, denseGETRS, denseAddIdentity

from scipy import linalg

from numpy.testing import assert_allclose

def testDenseGET():
    """
    Simple test of denseGETRF and denseGETRS to solve a linear system
    """
    N = 2
    A = np.array( [[4.0,6.0],[3.0,3.0]], dtype=np.float64, order='F')
    A_sp = np.ascontiguousarray(A)
    P = np.arange( N, dtype=np.int32 )
    b = np.ones( N, dtype=np.float64 )

    b_sp = linalg.lu_solve(linalg.lu_factor(A_sp), b)
    
    ret = denseGETRF(A, P)
    assert ret == 0
    
    denseGETRS(A, P, b)
    
    assert_allclose(b, b_sp)
    assert ret == 0
    
    
def testDenseAddIdentity():
    """
    Simple test of denseAddIdentity to check it works like the scipy
    equivalent.
    """
    
    N = 3
    I = np.identity(N)
    
    A_sp = np.random.random(size=(N,N))
    A = np.asfortranarray(A_sp)
    
    denseAddIdentity(A)
    
    assert_allclose(A, A_sp+I)
    
if __name__ == '__main__':
    
    testDenseGET()
    testDenseAddIdentity()