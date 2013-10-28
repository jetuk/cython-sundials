
from pySundials import sundials as csun

import numpy as np

class MyNvector(csun.N_Vector):
    def __init__(self, shp):
        #pass
        csun.N_Vector.__init__(self, )
        #
        self.shp = shp
        
        self.LengthRealType = np.prod(shp)
        self.LengthIntType = 1
        
        self.data = np.zeros(shp)
        
    def Size(self, ):
        return np.prod(self.shp)
        
    def Clone(self, ):
        nv = MyNvector(self.shp)
        return nv
        
    def CloneEmpty(self, ):
        return MyNvector(self.shp)
        
    def Destroy(self, ):
        pass
        
    def LinearSum(self, a, b, y, z):
        z.data = a*self.data + b*y.data
        
    def Constant(self, c):
        self.data[...] = c
        
    def Prod(self, y, z):
        z.data = self.data*y.data
        
    def Div(self, y, z):
        z.data = self.data/y.data
        
    def Scale(self, c, z):
        z.data = c*self.data
        
    def Abs(self, z):
        z.data = np.abs(self.data)
        
    def Inv(self, z):
        z.data = 1.0/self.data
        
    def AddConst(self, b, z):
        z.data = self.data + b
        
    def DotProd(self, y):
        return np.dot(self.data, y.data)
        
    def Min(self, ):
        return np.min(self.data)
        
    def MaxNorm(self, ):
        return np.max(np.abs(self.data))
        
    def WrmsNorm(self, w):
        return np.sqrt(np.sum( (self.data*w.data)**2)/self.Size() )          
    
class Motion(csun.Cvode):
    def RhsFn(self, t, y, ydot):
        #print 'Rhs', y.data, t
        ydot.data[0,:] = y.data[1,:]
        ydot.data[1,0] = -9.81
        ydot.data[1,1] = 0.0
        return 0
#nv = csun.N_Vector( )


nv = MyNvector( (2,2) )

v0 = 100
theta = np.pi/5

nv.data[0,0] = 0.0
nv.data[0,1] = 0.0
nv.data[1,0] = v0*np.sin(theta)
nv.data[1,1] = v0*np.cos(theta)

cv = Motion()
print cv
tout = np.linspace(0,10,50)

res = cv.Setup(nv, tout[0])


z = [nv.data[0,0], ]
x = [nv.data[0,1], ]
zd = [nv.data[1,0], ]
xd = [nv.data[1,1], ]


for t in tout[1:]:
    yout = MyNvector( 2 )
    res = cv.Solve(t, yout)
    z.append( yout.data[0,0] )
    x.append( yout.data[0,1] )
    zd.append( yout.data[1,0] )
    xd.append( yout.data[1,1] )
    
    
import pylab

pylab.subplot(2,1,1)
pylab.plot(x,z,'o')
pylab.grid()

pylab.subplot(2,1,2)
pylab.plot(tout, zd, '-o')
pylab.plot(tout, xd, '-o')
pylab.grid()

pylab.show()
