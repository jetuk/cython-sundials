


cdef class NvectorMemoryViewDouble5D(N_Vector):
    
    def __init__(self, shp, debug=False):
        #pass
        N_Vector.__init__(self, )
        
        assert len(shp) == 5
        self.shp = shp
        
        self.LengthRealType = np.prod(self.shp)
        self.LengthIntType = 1
        
        self.data = np.zeros(self.shp)
        
        self.debug = debug
        print self.data.shape[4]
        
    def Size(self, ):
        return np.prod(self.shp)
        
    def Clone(self, ):
        return NvectorMemoryViewDouble5D(self.shp, debug=self.debug) 
        
    def CloneEmpty(self, ):
        return NvectorMemoryViewDouble5D(self.shp, debug=self.debug)
        
    def Destroy(self, ):
        self.data = None
        
        N_Vector.Destroy(self, )
        
    def LinearSum(self, sun.realtype a, sun.realtype b, NvectorMemoryViewDouble5D y, NvectorMemoryViewDouble5D z):
        cdef int i0,i1,i2,i3,i4
        for i0 in range(self.data.shape[0]):
            for i1 in range(self.data.shape[1]):
                for i2 in range(self.data.shape[2]):
                    for i3 in range(self.data.shape[3]):
                        for i4 in range(self.data.shape[4]):
                            z.data[i0,i1,i2,i3,i4] = a*self.data[i0,i1,i2,i3,i4] + b*y.data[i0,i1,i2,i3,i4]
        if self.debug: print 'LinearSum', self.data[0,0,0,0,0], a, y.data[0,0,0,0,0], b, z.data[0,0,0,0,0]
        
    def Constant(self, sun.realtype c):
        self.data[...] = c
        if self.debug: print 'Constant', c
        
    def Prod(self, NvectorMemoryViewDouble5D y, NvectorMemoryViewDouble5D z):
        cdef int i0,i1,i2,i3,i4
        for i0 in range(self.data.shape[0]):
            for i1 in range(self.data.shape[1]):
                for i2 in range(self.data.shape[2]):
                    for i3 in range(self.data.shape[3]):
                        for i4 in range(self.data.shape[4]):
                            z.data[i0,i1,i2,i3,i4] = self.data[i0,i1,i2,i3,i4]*y.data[i0,i1,i2,i3,i4]        
        if self.debug: print 'Prod', self.data[0,0,0,0,0], y.data[0,0,0,0,0], z.data[0,0,0,0,0]
        
    @cython.cdivision(True)   
    def Div(self, NvectorMemoryViewDouble5D y, NvectorMemoryViewDouble5D z):
        cdef int i0,i1,i2,i3,i4
        for i0 in range(self.data.shape[0]):
            for i1 in range(self.data.shape[1]):
                for i2 in range(self.data.shape[2]):
                    for i3 in range(self.data.shape[3]):
                        for i4 in range(self.data.shape[4]):
                            z.data[i0,i1,i2,i3,i4] = self.data[i0,i1,i2,i3,i4]/y.data[i0,i1,i2,i3,i4]        
        if self.debug: print 'Div', self.data[0,0,0,0,0], y.data[0,0,0,0,0], z.data[0,0,0,0,0]
        
    def Scale(self, sun.realtype c, NvectorMemoryViewDouble5D z):
        cdef int i0,i1,i2,i3,i4
        for i0 in range(self.data.shape[0]):
            for i1 in range(self.data.shape[1]):
                for i2 in range(self.data.shape[2]):
                    for i3 in range(self.data.shape[3]):
                        for i4 in range(self.data.shape[4]):
                            z.data[i0,i1,i2,i3,i4] = c*self.data[i0,i1,i2,i3,i4]
        if self.debug: print 'Scale', self.data[0,0,0,0,0], z.data[0,0,0,0,0]                            
        
    def Abs(self, NvectorMemoryViewDouble5D z):
        cdef int i0,i1,i2,i3,i4
        for i0 in range(self.data.shape[0]):
            for i1 in range(self.data.shape[1]):
                for i2 in range(self.data.shape[2]):
                    for i3 in range(self.data.shape[3]):
                        for i4 in range(self.data.shape[4]):
                            z.data[i0,i1,i2,i3,i4] = abs(self.data[i0,i1,i2,i3,i4])
        if self.debug: print 'Abs', self.data[0,0,0,0,0], z.data[0,0,0,0,0]                            
    
    @cython.cdivision(True)    
    def Inv(self, NvectorMemoryViewDouble5D z):
        cdef int i0,i1,i2,i3,i4
        for i0 in range(self.data.shape[0]):
            for i1 in range(self.data.shape[1]):
                for i2 in range(self.data.shape[2]):
                    for i3 in range(self.data.shape[3]):
                        for i4 in range(self.data.shape[4]):
                            z.data[i0,i1,i2,i3,i4] = 1.0/self.data[i0,i1,i2,i3,i4]
        if self.debug: print 'Inv', self.data[0,0,0,0,0], z.data[0,0,0,0,0]                            
        
    def AddConst(self, sun.realtype b, NvectorMemoryViewDouble5D z):        
        cdef int i0,i1,i2,i3,i4
        for i0 in range(self.data.shape[0]):
            for i1 in range(self.data.shape[1]):
                for i2 in range(self.data.shape[2]):
                    for i3 in range(self.data.shape[3]):
                        for i4 in range(self.data.shape[4]):
                            z.data[i0,i1,i2,i3,i4] = self.data[i0,i1,i2,i3,i4] + b
        if self.debug: print 'AddConst', self.data[0,0,0,0,0], z.data[0,0,0,0,0]
        
    def DotProd(self, NvectorMemoryViewDouble5D y):
        cdef double ret = 0
        cdef int i0,i1,i2,i3,i4
        for i0 in range(self.data.shape[0]):
            for i1 in range(self.data.shape[1]):
                for i2 in range(self.data.shape[2]):
                    for i3 in range(self.data.shape[3]):
                        for i4 in range(self.data.shape[4]):
                            ret += self.data[i0,i1,i2,i3,i4]*y.data[i0,i1,i2,i3,i4]
        if self.debug: print 'DotProd', self.data[0,0,0,0,0] , y.data[0,0,0,0,0], ret                            
        return ret
        
    def Min(self, ):
        cdef double ret = self.data[0,0,0,0,0]        
        cdef int i0,i1,i2,i3,i4
        for i0 in range(self.data.shape[0]):
            for i1 in range(self.data.shape[1]):
                for i2 in range(self.data.shape[2]):
                    for i3 in range(self.data.shape[3]):
                        for i4 in range(self.data.shape[4]):
                            
                            ret = min(ret,self.data[i0,i1,i2,i3,i4])
        if self.debug: print 'Min', self.data[0,0,0,0,0], ret                            
        return ret
        
    def MaxNorm(self, ):

        cdef double ret = abs(self.data[0,0,0,0,0])
        cdef int i0,i1,i2,i3,i4
        for i0 in range(self.data.shape[0]):
            for i1 in range(self.data.shape[1]):
                for i2 in range(self.data.shape[2]):
                    for i3 in range(self.data.shape[3]):
                        for i4 in range(self.data.shape[4]):
                            ret = max(ret,abs(self.data[i0,i1,i2,i3,i4]))
        if self.debug: print 'MaxNorm', self.data[0,0,0,0,0], ret                            
        return ret
        
    def WrmsNorm(self, NvectorMemoryViewDouble5D w):
        cdef double ret = 0.0
        cdef int N = self.Size()
        cdef int i0,i1,i2,i3,i4
        for i0 in range(self.data.shape[0]):
            for i1 in range(self.data.shape[1]):
                for i2 in range(self.data.shape[2]):
                    for i3 in range(self.data.shape[3]):
                        for i4 in range(self.data.shape[4]):
                            ret += (self.data[i0,i1,i2,i3,i4]*w.data[i0,i1,i2,i3,i4])**2
        if self.debug: print 'WrmsNorm', self.data[0,0,0,0,0], np.sqrt(ret/N), N
        return np.sqrt(ret/N)
        
        
    def WL2Norm(self, NvectorMemoryViewDouble5D w):
        cdef double ret = 0.0
        cdef int N = self.Size()
        cdef int i0,i1,i2,i3,i4
        for i0 in range(self.data.shape[0]):
            for i1 in range(self.data.shape[1]):
                for i2 in range(self.data.shape[2]):
                    for i3 in range(self.data.shape[3]):
                        for i4 in range(self.data.shape[4]):
                            ret += (self.data[i0,i1,i2,i3,i4]*w.data[i0,i1,i2,i3,i4])**2
        if self.debug: print 'WL2Norm', self.data[0,0,0,0,0], w.data[0,0,0,0,0], ret
        return np.sqrt(ret)
        
    def L1Norm(self, ):
        
        cdef double ret = 0.0
        cdef int i0,i1,i2,i3,i4
        for i0 in range(self.data.shape[0]):
            for i1 in range(self.data.shape[1]):
                for i2 in range(self.data.shape[2]):
                    for i3 in range(self.data.shape[3]):
                        for i4 in range(self.data.shape[4]):
                            ret += abs(self.data[i0,i1,i2,i3,i4])
        if self.debug: print 'L1Norm', self.data[0,0,0,0,0], ret
        return ret
        