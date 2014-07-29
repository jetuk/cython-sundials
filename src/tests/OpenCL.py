# -*- coding: utf-8 -*-
"""
An example of using OpenCL with Sundials
"""



from pySundials.sundials import N_Vector
import numpy as np
import pyopencl as cl
import pyopencl.array
from pyopencl.elementwise import ElementwiseKernel






class NvectorCL(N_Vector):
    
    def __init__(self, ctx, cqa, shape, allocator=None):
        N_Vector.__init__(self,)
        
        self.cqa = cqa
        self.shape = shape
        self.allocator = allocator
        self.array = cl.array.empty(cqa, shape, dtype=np.float32, allocator=allocator)
        self.ctx = ctx
        self._create_kernels(ctx)
        
    def _create_kernels(self,ctx):
        
        self.linear_sum = ElementwiseKernel(ctx,
            "float a, float b, float *a_x, float *a_y, float *a_z",
            "a_z[i] = a * a_x[i] + b * a_y[i]",
            "lin_sum"
        )
        
    def Size(self, ):
        return np.prod(self.shp)
        
    def Clone(self, ):
        return NvectorCL(self.ctx, self.cqa, self.shape, self.allocator) 
        
    def CloneEmpty(self, ):
        return NvectorCL(self.ctx, self.cqa, self.shape, self.allocator) 
        
    def Destroy(self, ):
        self.data = None
        N_Vector.Destroy(self, )
        
    def LinearSum(self, a,  b, y, z):
        self.linear_sum(a,b,self.array, y.array, z.array)
        #z.data[...] = a*self.data + b*y.data
        
        
    def Constant(self, c):
        self.array.fill(c)
#        
#    def Prod(self, NvectorNdarrayFloat64 y, NvectorNdarrayFloat64 z):
#        z.data[...] = self.data*y.data
#        
#    def Div(self, NvectorNdarrayFloat64 y, NvectorNdarrayFloat64 z):
#        z.data[...] = self.data/y.data
#        
#    def Scale(self, sun.realtype c, NvectorNdarrayFloat64 z):
#        z.data[...] = c*self.data
#        
#    def Abs(self, NvectorNdarrayFloat64 z):
#        z.data[...] = np.abs(self.data)
#        
#    def Inv(self, NvectorNdarrayFloat64 z):
#        z.data[...] = 1.0/self.data
#        
#    def AddConst(self, sun.realtype b, NvectorNdarrayFloat64 z):
#        z.data[...] = self.data + b
#        
#    def DotProd(self, NvectorNdarrayFloat64 y):
#        return np.sum(self.data*y.data)
#        
#    def Min(self, ):
#        return np.min(self.data)
#        
#    def MaxNorm(self, ):
#        return np.max(np.abs(self.data))
#        
#    def WrmsNorm(self, NvectorNdarrayFloat64 w):
#        return np.sqrt(np.sum( (self.data*w.data)**2)/self.Size() )    
#        
#    def WL2Norm(self, NvectorNdarrayFloat64 w):
#        return np.sqrt(np.sum( (self.data*w.data)**2 ) )
#        
#    def L1Norm(self, ):
#        return np.sum(np.abs(self.data))
        
        
        
if __name__ == '__main__':
    
    
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    x = NvectorCL(ctx, queue, 10000,)
    y = x.Clone()
    z = x.Clone()
    
    x.Constant(10.0)
    y.Constant(5.0)

    x.LinearSum(1.0,2.0,y,z)    
    
    print z.array.get().shape
    