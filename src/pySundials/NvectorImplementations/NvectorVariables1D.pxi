


def common_entries(*dcts):
    #for i in set(dcts[0]).intersection(*dcts[1:]):
    #    yield (i,) + tuple(d[i] for d in dcts)
    for i in dcts[0]:
        yield (i,) + tuple(d[i] for d in dcts)


cdef class Variable1D(object):

    
    def __init__(self, name, value=None, shape=None, copied_attrs=None):
        self.name = name
        
        if not value is None:
            self.value = np.asarray(value)
        elif not shape is None: 
            self.value = np.zeros( shape )
        else:
            raise ValueError('Either keyword argument value or shape must not be None')
        
        self.saved_values = {}
        self._copied_attr_names = []
        self._attrs = {}
        if not copied_attrs is None:
            for name, val in copied_attrs.iteritems():
                #setattr(self, name, val)
                self._attrs[name] = val
                self._copied_attr_names.append(name)
                
    def __getattr__(self, name):
        if name in self._attrs:
            return self._attrs[name]
        return object.__getattribute__(self, name)

    def cvalue(self, t):
        return self.value

    def copy(self, ):
        #print self.value.shape, self.name, type(self.value)
        return Variable1D(self.name, value=self.value.copy(),
                        copied_attrs={n:getattr(self,n) for n in self._copied_attr_names})
        
    #def copyEmpty(self, ):
        
    def saveValue(self, t):
        self.saved_values[t] = self.value.copy()
        
    def savedLimits(self, ):
        mn = np.min([vals.min() for vals in self.saved_values.itervalues()])
        mx = np.max([vals.max() for vals in self.saved_values.itervalues()])
        return mn,mx
        
        
cdef class Variable1DNvector(N_Vector):
    def __init__(self, variables):
        N_Vector.__init__(self, )
        
        self.variables = variables
        
        self.LengthRealType = self.Size()
        self.LengthIntType = 1

        
    def Size(self, ):
        cdef Variable1D var
        cdef int N
        N = 0
        for var in self.variables:
            N += var.value.shape[0]
        return N
        
    def Clone(self, ):
        return Variable1DNvector([var.copy() for var in self.variables])

        
    def CloneEmpty(self, ):
        return Variable1DNvector([var.copy() for var in self.variables])
        
    def Destroy(self, ):
        pass
        
    def LinearSum(self, double a, double b, N_Vector y, N_Vector z): 
        cdef Variable1D vx, vy, vz
        cdef int I, i

        for vx, vy, vz in zip(self.variables, y.variables, z.variables):
            I = vx.value.shape[0]
            for i in range(I):
                vz.value[i] = a*vx.value[i] + b*vy.value[i]
        
    def Constant(self, c):
        cdef Variable1D vx
        cdef int I, i        
        for vx in self.variables:
            I = vx.value.shape[0]
            for i in range(I):
                vx.value[i] = c
        
    def Prod(self, y, z):
        cdef Variable1D vx, vy, vz
        cdef int I, i
        for vx, vy, vz in zip(self.variables, y.variables, z.variables):
            I = vx.value.shape[0]
            for i in range(I):            
                vz.value[i] = vx.value[i]*vy.value[i]
        
    def Div(self, y, z):
        cdef Variable1D vx, vy, vz
        cdef int I, i
        for vx, vy, vz in zip(self.variables, y.variables, z.variables):         
            I = vx.value.shape[0]
            for i in range(I):            
                vz.value[i] = vx.value[i]/vy.value[i]
        
        
    def Scale(self, double c, z):
        cdef Variable1D vx, vz
        cdef int I, i
        for vx, vz in zip(self.variables, z.variables):  
            I = vx.value.shape[0]
            for i in range(I):               
                vz.value[i] = vx.value[i]*c        
        
        
    def Abs(self, z):
        cdef Variable1D vx, vz    
        cdef int I, i
        for vx, vz in zip(self.variables, z.variables):          
            I = vx.value.shape[0]
            for i in range(I):               
                vz.value[i] = abs(vx.value[i])
        
        
    def Inv(self, z):
        cdef Variable1D vx, vz    
        cdef int I, i
        for vx, vz in zip(self.variables, z.variables):          
            I = vx.value.shape[0]
            for i in range(I):               
                vz.value[i] = 1.0/vx.value[i]
                
    def AddConst(self, double b, z):
        cdef Variable1D vx, vz    
        cdef int I, i
        for vx, vz in zip(self.variables, z.variables):          
            I = vx.value.shape[0]
            for i in range(I):               
                vz.value[i] = vx.value[i] + b
        
    def DotProd(self, y):
        cdef double ret = 0.0  
        cdef Variable1D vx, vy
        cdef int I, i
        for vx, vy in zip(self.variables, y.variables):          
            I = vx.value.shape[0]
            for i in range(I):   
                ret += vx.value[i]*vy.value[i]

        return ret
        
    def Min(self, ):
        #print 'Min', np.min([var.value.min() for var in self.variables.itervalues()])
        cdef Variable1D var
        return np.min([np.asarray(var.value).min() for var in self.variables])
        
        
    def MaxNorm(self, ):
        cdef Variable1D var
        return np.max([np.abs(np.asarray(var.value)).max() for var in self.variables])
        
    def WrmsNorm(self, w):        
        cdef double ret = 0.0

        cdef Variable1D vx, vw
        cdef int I, i
        for vx, vw in zip(self.variables, w.variables):         
            I = vx.value.shape[0]
            for i in range(I):               
                ret += (vx.value[i]*vw.value[i])**2 
        #print 'WrmsNorm', ret, np.sqrt(ret/self.Size()), type(np.sqrt(ret/self.Size()))
        return np.sqrt(ret/self.Size())