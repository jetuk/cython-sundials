# denseGETRF and denseGETRS Python wrappers included in the 
# separate sundials solvers



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