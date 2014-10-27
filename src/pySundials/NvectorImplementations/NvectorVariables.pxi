


def common_entries(*dcts):
    #for i in set(dcts[0]).intersection(*dcts[1:]):
    #    yield (i,) + tuple(d[i] for d in dcts)
    for i in dcts[0]:
        yield (i,) + tuple(d[i] for d in dcts)


cdef class Variable(object):
    cdef public str name
    cdef public dict saved_values
    cdef list _copied_attr_names
    cdef public np.ndarray value
    cdef dict _attrs    
    
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
        return Variable(self.name, value=self.value.copy(),
                        copied_attrs={n:getattr(self,n) for n in self._copied_attr_names})
        
    #def copyEmpty(self, ):
        
    def saveValue(self, t):
        self.saved_values[t] = self.value.copy()
        
    def savedLimits(self, ):
        mn = np.min([vals.min() for vals in self.saved_values.itervalues()])
        mx = np.max([vals.max() for vals in self.saved_values.itervalues()])
        return mn,mx
        
        
cdef class VariableNvector(N_Vector):
    def __init__(self, variables={}):
        N_Vector.__init__(self, )
        
        self.variables = variables
        
        self.LengthRealType = self.Size()
        self.LengthIntType = 1
        
        
    def addVariable(self, var):       
        self.variables[var.name] = var
        self.LengthRealType = self.Size()
        self.LengthIntType = 1
        
    def Size(self, ):
        cdef Variable var
        return np.sum([var.value.size for var in self.variables.itervalues()])
        
    def Clone(self, ):
        return VariableNvector(dict([(name,var.copy()) for name,var in self.variables.iteritems()]))

        
    def CloneEmpty(self, ):
        return VariableNvector(dict([(name,var.copy()) for name,var in self.variables.iteritems()]))
        
    def Destroy(self, ):
        pass
        
    def LinearSum(self, double a, double b, N_Vector y, N_Vector z):
        #for vname, vx, vy, vz in common_entries(self.variables, y.variables, z.variables):
        cdef str key    
        cdef Variable vx, vy, vz
        for key in self.variables:
            vx = self.variables[key]
            vy = y.variables[key]
            vz = z.variables[key]            
            vz.value[...] = a*vx.value[...] + b*vy.value[...]
        
    def Constant(self, c):
        cdef Variable vx
        for vx in self.variables.itervalues():
            vx.value[...] = c
        
    def Prod(self, y, z):
        cdef str key   
        cdef Variable vx, vy, vz
        for key in self.variables:
        #for vname, vx, vy, vz in common_entries(self.variables, y.variables, z.variables):
            vx = self.variables[key]
            vy = y.variables[key]
            vz = z.variables[key]
            vz.value[...] = vx.value*vy.value
        
    def Div(self, y, z):
#        for vname, vx, vy, vz in common_entries(self.variables, y.variables, z.variables):
        cdef str key       
        cdef Variable vx, vy, vz
        for key in self.variables:
        #for vname, vx, vy, vz in common_entries(self.variables, y.variables, z.variables):
            vx = self.variables[key]
            vy = y.variables[key]
            vz = z.variables[key]            
            vz.value[...] = vx.value/vy.value       
        
        
    def Scale(self, double c, VariableNvector z):
        #print 'Scale', c
        #for vname, vx, vz in common_entries(self.variables, z.variables):
        cdef str key   
        cdef Variable vx, vz
        for key in self.variables:
            vx = self.variables[key]
            vz = z.variables[key] 
            vz.value[...] = vx.value*c        
        
        
    def Abs(self, z):
        #print 'Abs'
        #for vname, vx, vz in common_entries(self.variables, z.variables):
        cdef str key   
        cdef Variable vx, vz    
        for key in self.variables:
            vx = self.variables[key]
            vz = z.variables[key]             
            vz.value[...] = np.abs(vx.value)
        
        
    def Inv(self, z):
        #print 'Inv'
        #for vname, vx, vz in common_entries(self.variables, z.variables):
        cdef str key   
        cdef Variable vx, vz    
        for key in self.variables:
            vx = self.variables[key]
            vz = z.variables[key]             
            vz.value[...] = 1.0/vx.value
                
    def AddConst(self, b, z):
        #print 'AddConst', b
        #for vname, vx, vz in common_entries(self.variables, z.variables):
        cdef str key   
        cdef Variable vx, vz    
        for key in self.variables:
            vx = self.variables[key]
            vz = z.variables[key]             
            vz.value[...] = vx.value + b
        
    def DotProd(self, y):
        cdef float ret = 0.0
#        for vname, vx, vy in common_entries(self.variables, y.variables):
        cdef str key   
        cdef Variable vx, vy        
        for key in self.variables:
            vx = self.variables[key]
            vy = y.variables[key]
 
            ret += np.dot(vx.value,vy.value)
        #print 'DotProd', ret
        return ret
        
    def Min(self, ):
        #print 'Min', np.min([var.value.min() for var in self.variables.itervalues()])
        cdef Variable var
        return np.min([var.value.min() for var in self.variables.itervalues()])
        
        
    def MaxNorm(self, ):
        cdef Variable var
        return np.max([np.max(np.abs(var.value)) for var in self.variables.itervalues()])
        
    def WrmsNorm(self, VariableNvector w):        
        cdef float ret = 0.0
        #for vname, vx, vw in common_entries(self.variables, w.variables):
        cdef str key
        cdef Variable vx, vw
        for key in self.variables:
            vx = self.variables[key]
            vw = w.variables[key]            
            ret += np.sum( (vx.value*vw.value)**2 )
        #print 'WrmsNorm', ret, np.sqrt(ret/self.Size()), type(np.sqrt(ret/self.Size()))
        return np.sqrt(ret/self.Size())