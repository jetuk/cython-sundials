# denseGETRF and denseGETRS Python wrappers included in the 
# separate sundials solvers



cpdef int denseGETRF( sun.realtype[:,::1] a, long int[::1] p ):
    
    cdef long int m = a.shape[0]
    cdef long int n = a.shape[1]
    cdef long int ret    
    cdef long int i
    
    cdef sun.realtype **col_ptrs = <sun.realtype**>malloc(n*cython.sizeof(<sun.realtype*>NULL))
    
    for i in range(n):
        col_ptrs[i] = &a[i,0]

    ret = sun.denseGETRF( col_ptrs, m, n, &p[0])    
    free(col_ptrs)
    
    return ret
    
    
cpdef denseGETRS(sun.realtype[:,::1] a, long[::1] p, sun.realtype[::1] b ):
    
    cdef int n = a.shape[0]
        
    cdef sun.realtype **col_ptrs = <sun.realtype**>malloc(n*cython.sizeof(<sun.realtype*>NULL))
    
    for i in range(n):
        col_ptrs[i] = &a[i,0]    
    
    sun.denseGETRS(col_ptrs, n, &p[0], &b[0] )
    
    