


cdef class NvectorNdarrayFloat64(N_Vector):
    #cdef public np.ndarray data
    #cdef object shp
    #cdef public double[:,::1] data
    
    def __init__(self, shp):
        #pass
        N_Vector.__init__(self, )
        
        self.shp = shp
        
        self.LengthRealType = np.prod(self.shp)
        self.LengthIntType = 0
        
        self.data = np.zeros(self.shp)
        
    def Size(self, ):
        return np.prod(self.shp)
        
    def Clone(self, ):
        return NvectorNdarrayFloat64(self.shp)        
        
    def CloneEmpty(self, ):
        return NvectorNdarrayFloat64(self.shp)
        
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
        
        
    def LinearSum(self, sun.realtype a, sun.realtype b, NvectorNdarrayFloat64 y, NvectorNdarrayFloat64 z):
        z.data[...] = a*self.data + b*y.data
        
        
    def Constant(self, sun.realtype c):
        self.data[...] = c
        
    def Prod(self, NvectorNdarrayFloat64 y, NvectorNdarrayFloat64 z):
        z.data[...] = self.data*y.data
        
    def Div(self, NvectorNdarrayFloat64 y, NvectorNdarrayFloat64 z):
        z.data[...] = self.data/y.data
        
    def Scale(self, sun.realtype c, NvectorNdarrayFloat64 z):
        z.data[...] = c*self.data
        
    def Abs(self, NvectorNdarrayFloat64 z):
        z.data[...] = np.abs(self.data)
        
    def Inv(self, NvectorNdarrayFloat64 z):
        z.data[...] = 1.0/self.data
        
    def AddConst(self, sun.realtype b, NvectorNdarrayFloat64 z):
        z.data[...] = self.data + b
        
    def DotProd(self, NvectorNdarrayFloat64 y):
        return np.sum(self.data*y.data)
        
    def Min(self, ):
        return np.min(self.data)
        
    def MaxNorm(self, ):
        return np.max(np.abs(self.data))
        
    def WrmsNorm(self, NvectorNdarrayFloat64 w):
        return np.sqrt(np.sum( (self.data*w.data)**2)/self.Size() )    
        
    def WL2Norm(self, NvectorNdarrayFloat64 w):
        return np.sqrt(np.sum( (self.data*w.data)**2 ) )
        
    def L1Norm(self, ):
        return np.sum(np.abs(self.data))
        
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
        
        return m.data.max() == 0.0
        
        
        
        