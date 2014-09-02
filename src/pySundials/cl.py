# -*- coding: utf-8 -*-
"""
OpenCL implementation of NVector


"""

from sundials import N_Vector

import numpy as np
import pyopencl as cl
import pyopencl.array
from pyopencl.elementwise import ElementwiseKernel
from pyopencl.reduction import ReductionKernel



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
            "nv_lin_sum"
        )
        
        self.prod = ElementwiseKernel(ctx,
            "float *a_x, float *a_y, float *a_z",
            "a_z[i] = a_x[i] * a_y[i]",
            "nv_prod"
        )
        self.div = ElementwiseKernel(ctx,
            "float *a_x, float *a_y, float *a_z",
            "a_z[i] = a_x[i] / a_y[i]",
            "nv_div"
        )        
        self.scale = ElementwiseKernel(ctx,
            "float c, float *a_x, float *a_z",
            "a_z[i] = c*a_x[i]",
            "nv_scale"
        )        
        self.abs = ElementwiseKernel(ctx,
            "float *a_x, float *a_z",
            "a_z[i] = fabs(a_x[i])",
            "nv_abs"
        )                                             
        
        self.inv = ElementwiseKernel(ctx, 
            "float *a_x, float *a_z",
            "a_z[i] = 1.0/a_x[i]",
            "nv_inv"
        )   
        
        self.add_const = ElementwiseKernel(ctx,
            "float b, float *a_x, float *a_z",
            "a_z[i] = b+a_x[i]",
            "nv_add_const"
        )   
        
        
        self.sum_sq = ReductionKernel(ctx, np.float32,
                neutral="0.0",
                reduce_expr="a+b", map_expr="x[i]*x[i]*y[i]*y[i]",
                arguments="float *x, float *y",
                name = "nv_sum_sq",
        )
        
        self.sum_abs = ReductionKernel(ctx, np.float32,
                neutral="0.0",
                reduce_expr="a+b", map_expr="fabs(x[i])",
                arguments="float *x",
                name = "nv_sum_abs",
        )        
        
    def Size(self, ):
        return np.prod(self.shape)
        
    def Clone(self, ):
        return NvectorCL(self.ctx, self.cqa, self.shape, self.allocator) 
        
    def CloneEmpty(self, ):
        return NvectorCL(self.ctx, self.cqa, self.shape, self.allocator) 
        
    def Destroy(self, ):
        self.data = None
        N_Vector.Destroy(self, )
        
    def LinearSum(self, a,  b, y, z):
        self.linear_sum(a,b,self.array, y.array, z.array)       
        
    def Constant(self, c):
        self.array.fill(c)
        
    def Prod(self, y, z):
        self.prod(self.array, y.array, z.array)

    def Div(self, y, z):
        self.div(self.array, y.array, z.array)
        
    def Scale(self, c, z):
        self.scale(c, self.array, z.array)
        
    def Abs(self, z):
        self.abs(self.array, z.array)
        
    def Inv(self, z):
        self.inv(self.array, z.array)

    def AddConst(self, b, z):
        self.add_const(b, self.array, z.array)
        
    def DotProd(self, y):
        return pyopencl.array.dot(self.array, y.array).get()
        
    def Min(self, ):
        return pyopencl.array.min(self.array).get()
        
    def MaxNorm(self, ):
        mn = pyopencl.array.min(self.array).get()
        mx = pyopencl.array.max(self.array).get()
        return max( [abs(mn), abs(mx)] )
        
    def WrmsNorm(self, w):
        ssq = self.sum_sq(self.array, w.array).get()
        return np.sqrt(ssq/self.Size() )    
        
    def WL2Norm(self, w):
        ssq = self.sum_sq(self.array, w.array).get()
        return np.sqrt( ssq )
        
    def L1Norm(self, ):
        return self.sum_abs(self.array).get()