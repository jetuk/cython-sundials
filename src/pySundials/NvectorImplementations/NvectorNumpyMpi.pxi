


cdef class NvectorMPINdarrayFloat64(N_Vector):
    #cdef public np.ndarray data
    #cdef object shp
    #cdef public double[:,::1] data
    
    def __init__(self, comm, local_shp, global_shp):
        N_Vector.__init__(self, )
        
        self.local_shp = local_shp
        self.global_shp = global_shp
        self.comm = comm
        
        npes = comm.Get_size()
        self.LengthRealType = np.prod(self.global_shp)
        self.LengthIntType = 2*npes
        
        self.data = np.zeros(self.local_shp)
        
    def Size(self, ):
        return np.prod(self.global_shp)
        
    def Clone(self, ):
        return NvectorMPINdarrayFloat64(self.comm, self.local_shp, self.global_shp)
        
    def CloneEmpty(self, ):
        return NvectorMPINdarrayFloat64(self.comm, self.local_shp, self.global_shp)
        
    def Destroy(self, ):
        self.data = None
        
        N_Vector.Destroy(self, )
        
    def GetArrayPointer(self, ):
        #Py_INCREF(self.data)
        return <uintptr_t>self.data.data
        
    def SetArrayPointer(self, uintptr_t ptr):
        cdef int ndims = np.PyArray_NDIM(self.data)
        cdef np.npy_intp* shp = np.PyArray_DIMS(self.data)
        Py_INCREF(self.data)
        self.data = np.PyArray_SimpleNewFromData(ndims, shp, np.NPY_DOUBLE, <void*>ptr)
        
    def LinearSum(self, sun.realtype a, sun.realtype b, NvectorMPINdarrayFloat64 y, NvectorMPINdarrayFloat64 z):
        z.data[...] = a*self.data + b*y.data
        
        
    def Constant(self, sun.realtype c):
        self.data[...] = c
        
    def Prod(self, NvectorMPINdarrayFloat64 y, NvectorMPINdarrayFloat64 z):
        z.data[...] = self.data*y.data
        
    def Div(self, NvectorMPINdarrayFloat64 y, NvectorMPINdarrayFloat64 z):
        z.data[...] = self.data/y.data
        
    def Scale(self, sun.realtype c, NvectorMPINdarrayFloat64 z):
        z.data[...] = c*self.data
        
    def Abs(self, NvectorMPINdarrayFloat64 z):
        z.data[...] = np.abs(self.data)
        
    def Inv(self, NvectorMPINdarrayFloat64 z):
        z.data[...] = 1.0/self.data
        
    def AddConst(self, sun.realtype b, NvectorMPINdarrayFloat64 z):
        z.data[...] = self.data + b
        
    def DotProd(self, NvectorMPINdarrayFloat64 y):
        lsum = np.sum(self.data*y.data)
        gsum = np.array(0.0, 'd')
        self.comm.Allreduce(lsum, gsum)
        return gsum
        
    def Min(self, ):
        lmin = np.min(self.data)
        gmin = np.array(0.0, 'd')
        self.comm.Allreduce(lmin,gmin,op=MPI.MIN)
        return gmin
        
    def MaxNorm(self, ):
        lmax = np.max(np.abs(self.data))
        gmax = np.array(0.0, 'd')
        self.comm.Allreduce(lmax,gmax,op=MPI.MAX)
        return gmax
        
    def WrmsNorm(self, NvectorMPINdarrayFloat64 w):
        lsum = np.sum( (self.data*w.data)**2)
        gsum = np.array(0.0, 'd')
        self.comm.Allreduce(lsum,gsum,op=MPI.SUM)
        return np.sqrt( gsum/self.Size() )    
        
    def WL2Norm(self, NvectorMPINdarrayFloat64 w):
        lsum = np.sum( (self.data*w.data)**2)
        gsum = np.array(0.0, 'd')
        self.comm.Allreduce(lsum,gsum,op=MPI.SUM)
        return np.sqrt( gsum )  
        
    def L1Norm(self, ):
        lsum = np.sum(np.abs(self.data))
        gsum = np.array(0.0, 'd')
        self.comm.Allreduce(lsum,gsum,op=MPI.SUM)
        return gsum
        
    def ConstrMask(self, c, m):
        """
         * N_VConstrMask
         *   Performs the operation : 
         *       m[i] = 1.0 if constraint test fails for x[i]
         *       m[i] = 0.0 if constraint test passes for x[i]
         *   where the constraint tests are as follows:
         *      If c[i] = +2.0, then x[i] must be >  0.0.
         *      If c[i] = +1.0, then x[i] must be >= 0.0.
         *      If c[i] = -1.0, then x[i] must be <= 0.0.
         *      If c[i] = -2.0, then x[i] must be <  0.0.
         *   This routine returns a boolean FALSE if any element failed
         *   the constraint test, TRUE if all passed. It also sets a
         *   mask vector m, with elements equal to 1.0 where the
         *   corresponding constraint test failed, and equal to 0.0
         *   where the constraint test passed.
         *   This routine is specialized in that it is used only for
         *   constraint checking.
        """
        m.data[...] = 0.0
        
        ind = np.where(c.data == 2.0)        
        m.data[ind] = 1.0-(self.data[ind] > 0.0)
        
        ind = np.where(c.data == 1.0)        
        m.data[ind] = 1.0-(self.data[ind] >= 0.0)
        
        ind = np.where(c.data == -1.0)        
        m.data[ind] = 1.0-(self.data[ind] <= 0.0)
        
        ind = np.where(c.data == -2.0)        
        m.data[ind] = 1.0-(self.data[ind] < 0.0)
        
        lmax  = m.data.max()
        gmax = np.array(0.0, 'd')
        self.comm.Allreduce(lmax,gmax,op=MPI.MAX)
        return gmax == 0.0
        
        
        
        