# denseGETRF and denseGETRS Python wrappers included in the 
# separate sundials solvers

cpdef int denseGETRF( sun.realtype[::1,:] a, long int[::1] p ):
    
    cdef long int m = a.shape[0]
    cdef long int n = a.shape[1]
    cdef long int ret    
    cdef long int i
    
    cdef sun.realtype **col_ptrs = <sun.realtype**>malloc(n*cython.sizeof(<sun.realtype*>NULL))
    
    for i in range(n):
        col_ptrs[i] = &a[0,i]

    ret = sun.denseGETRF( col_ptrs, m, n, &p[0])    
    free(col_ptrs)
    
    return ret
    
    
cpdef denseGETRS(sun.realtype[::1,:] a, long[::1] p, sun.realtype[::1] b ):
    
    cdef long int m = a.shape[0]
    cdef long int n = a.shape[1]
    cdef long int ret    
    cdef long int i
    
    cdef sun.realtype **col_ptrs = <sun.realtype**>malloc(n*cython.sizeof(<sun.realtype*>NULL))
    
    for i in range(n):
        col_ptrs[i] = &a[0,i]  
    
    sun.denseGETRS(col_ptrs, n, &p[0], &b[0] )
    
    
    
cpdef densePOTRF(sun.realtype[::1,:] a):    
    """
     DensePOTRF computes the Cholesky factorization of a real symmetric
     positive definite matrix A.
    """    
    cdef long int m = a.shape[0]
    cdef long int n = a.shape[1]
    cdef long int ret    
    cdef long int i
    
    cdef sun.realtype **col_ptrs = <sun.realtype**>malloc(n*cython.sizeof(<sun.realtype*>NULL))
    
    for i in range(n):
        col_ptrs[i] = &a[0,i]   
    
    ret = sun.densePOTRF(col_ptrs, n)
    free(col_ptrs)
    return ret
    
cpdef densePOTRS(sun.realtype[:,::1] a, sun.realtype[::1] b):
    """
     DensePOTRS solves a system of linear equations A*X = B with a 
     symmetric positive definite matrix A using the Cholesky factorization
     A = L*L**T computed by DensePOTRF.
    """
    cdef long int m = a.shape[0]
    cdef long int n = a.shape[1]
    cdef long int ret    
    cdef long int i
    
    cdef sun.realtype **col_ptrs = <sun.realtype**>malloc(n*cython.sizeof(<sun.realtype*>NULL))
    
    for i in range(n):
        col_ptrs[i] = &a[0,i] 
        
    sun.densePOTRS(col_ptrs, n, &b[0])
    free(col_ptrs)
    
cpdef denseGEQRF(sun.realtype[:,::1] a, sun.realtype[::1] beta, sun.realtype[::1] v):
    """
     DenseGEQRF computes a QR factorization of a real M-by-N matrix A:
      A = Q * R (with M>= N).

     DenseGEQRF requires a temporary work vector wrk of length M.
    """
    cdef long int m = a.shape[0]
    cdef long int n = a.shape[1]
    cdef long int ret    
    cdef long int i
    
    cdef sun.realtype **col_ptrs = <sun.realtype**>malloc(n*cython.sizeof(<sun.realtype*>NULL))
    
    for i in range(n):
        col_ptrs[i] = &a[0,i]  

    ret = sun.denseGEQRF(col_ptrs, m, n, &beta[0], &v[0])
    free(col_ptrs)
    return ret
    
cpdef denseORMQR(sun.realtype[::1,:] a, sun.realtype[::1] beta, sun.realtype[::1] v, 
                 sun.realtype[::1] w, sun.realtype[::1] wrk):
    """
     DenseORMQR computes the product w = Q * v where Q is a real 
     orthogonal matrix defined as the product of k elementary reflectors
    
            Q = H(1) H(2) . . . H(k)
    
     as returned by DenseGEQRF. Q is an M-by-N matrix, v is a vector
     of length N and w is a vector of length M (with M>=N).
    
     DenseORMQR requires a temporary work vector wrk of length M.
    """
    cdef long int m = a.shape[0]
    cdef long int n = a.shape[1]
    cdef long int ret    
    cdef long int i
    
    cdef sun.realtype **col_ptrs = <sun.realtype**>malloc(n*cython.sizeof(<sun.realtype*>NULL))
    
    for i in range(n):
        col_ptrs[i] = &a[0,i]  
        
    ret = sun.denseORMQR(col_ptrs, m, n, &beta[0], &v[0], &w[0], &wrk[0])
    free(col_ptrs)
    return ret
    
cpdef denseCopy(sun.realtype[:,::1] a, sun.realtype[:,::1] b):
    """
     DenseCopy copies the contents of the M-by-N matrix A into the
     M-by-N matrix B.
     
     DenseCopy is a wrapper around denseCopy which accesses the data
     in the DlsMat A and B (i.e. the fields cols)
    """
    cdef long int m = a.shape[0]
    cdef long int n = a.shape[1]
    cdef long int ret    
    cdef long int i
    
    cdef sun.realtype **col_ptrs_a = <sun.realtype**>malloc(n*cython.sizeof(<sun.realtype*>NULL))
    cdef sun.realtype **col_ptrs_b = <sun.realtype**>malloc(n*cython.sizeof(<sun.realtype*>NULL))
    
    for i in range(n):
        col_ptrs_a[i] = &a[0,i]      
        col_ptrs_b[i] = &a[0,i]

    sun.denseCopy(col_ptrs_a, col_ptrs_b, m, n)
    free(col_ptrs_a)    
    free(col_ptrs_b)    
    
cpdef denseScale(sun.realtype c, sun.realtype[::1,:] a):
    """
     DenseScale scales the elements of the M-by-N matrix A by the
     constant c and stores the result back in A.
    
     DenseScale is a wrapper around denseScale which performs the actual
     scaling by accessing the data in the DlsMat A (i.e. the field
     cols).
    """
    
    cdef long int m = a.shape[0]
    cdef long int n = a.shape[1]
    cdef long int ret    
    cdef long int i
    
    cdef sun.realtype **col_ptrs = <sun.realtype**>malloc(n*cython.sizeof(<sun.realtype*>NULL))
    
    for i in range(n):
        col_ptrs[i] = &a[0,i] 
        
    
    sun.denseScale(c, col_ptrs, m, n)
    free(col_ptrs)
    
cpdef denseAddIdentity(sun.realtype[::1,:] a):
    """
     denseAddIdentity adds the identity matrix to the n-by-n matrix
     stored in the realtype** arrays.
    """
    cdef long int m = a.shape[0]
    cdef long int n = a.shape[1]
    cdef long int ret    
    cdef long int i
    
    cdef sun.realtype **col_ptrs = <sun.realtype**>malloc(n*cython.sizeof(<sun.realtype*>NULL))
    
    for i in range(n):
        col_ptrs[i] = &a[0,i]
        
    sun.denseAddIdentity(col_ptrs, n)
    free(col_ptrs)