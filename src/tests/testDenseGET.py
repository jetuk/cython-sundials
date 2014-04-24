# -*- coding: utf-8 -*-


import numpy as np
from pySundials.cvode import denseGETRF, denseGETRS

from scipy import linalg

def testDenseGET():
    """
    Simple test of denseGETRF to ensure no memory corruption
    """
    N = 2
    A = np.array( [[4.0,6.0],[3.0,3.0]], dtype=np.float64)
    P = np.arange( N, dtype=np.int64 )
    b = np.ones( N, dtype=np.float64 )
    print linalg.lu_factor(A)
    print linalg.lu_solve(linalg.lu_factor(A.T), b)
    
    
    print A, P
    
    ret = denseGETRF(A, P)
    
    print A, P
    
    
    
    denseGETRS(A, P, b)
    
    print A, P, b, ret
    assert ret == 0
    
    
    
if __name__ == '__main__':
    
    testDenseGET()