# -*- coding: utf-8 -*-
"""
An example of using OpenCL with Sundials
"""



from pySundials.cl import NvectorCL
from pySundials.sundials import NvectorNdarrayFloat64
from pySundials.cvode import Cvode
import numpy as np
import pyopencl as cl
import pyopencl.array
from pyopencl.elementwise import ElementwiseKernel

from examples.cvode.testCvAdvDiff_bnd import AdvDiff

class AdvDiffVec(AdvDiff):
    
    def RhsFn(self, t, u, udot):
        """
        f routine. Compute f(t,u). 
         
        """
        # Extract needed constants from data 
        hordc = self.hdcoef
        horac = self.hacoef
        verdc = self.vdcoef
        
        mx,my = u.data.shape
        
        udot.data[...] = 0.0
        
        udot.data[:,1:] += (hordc-horac)*u.data[:,:-1]
        udot.data[:,:] += -hordc*TWO*u.data[:,:]
        udot.data[:,:-1] += (hordc+horac)*u.data[:,1:]
        
        udot.data[1:,:] += verdc*u.data[:-1,:]
        udot.data[:,:] += -verdc*TWO*u.data[:,:]
        udot.data[:-1,:] += verdc*u.data[1:,:]        
            
#        # Loop over all grid points. 
#        for j in range(my):
#            
#            for i in range(mx):
#                # Extract u at x_i, y_j and four neighboring points 
#                uij = u.data[i,j]
#                
#                udn = ZERO if j == 0 else u.data[i,j-1]
#                uup = ZERO if j == my-1 else u.data[i,j+1]
#                ult = ZERO if i == 0 else u.data[i-1,j]
#                urt = ZERO if i == mx-1 else u.data[i+1,j]
#                
#
#                # Set diffusion and advection terms and load into udot 
#
#                hdiff = hordc*(ult - TWO*uij + urt)
#                hadv = horac*(urt - ult)
#                vdiff = verdc*(uup - TWO*uij + udn)
#                udot.data[i,j] = hdiff + hadv + vdiff
  

        return 0 
        
class AdvDiffCl(AdvDiff):
    """
    UserData (contains grid constants) 
    
    """
    
    def __init__(self, ctx, mx, my, dx, dy, hdcoef, hacoef, vdcoef, **kwds):
        AdvDiff.__init__(self, dx, dy, hdcoef, hacoef, vdcoef, **kwds)
            
        self.mx = mx
        self.my = my

        
        self.advdiff = ElementwiseKernel(ctx, 
            "int MX, int MY, float hdcoef, float hacoef, float vdcoef, float *u, float *udot",
            """        
            float const ZERO = 0.0;
            
            int ix = i / MY;
            int iy = i % MY;
            
            float udn = iy == 0 ? ZERO : u[i-1];
            float uup = iy == MY-1 ? ZERO : u[i+1];
            float ult = ix == 0 ? ZERO : u[i-MY];
            float urt = ix == MX-1 ? ZERO : u[i+MY];
            
            udot[i] = hdcoef*(ult - 2.0*u[i] + urt) + hacoef*(urt - ult) + vdcoef*(uup - 2.0*u[i] + udn);
            
            """,
            "advdiff"
        )        
        
        self.initial = ElementwiseKernel(ctx, 
            "int MX, int MY, float coef, float dx, float dy, float *u",
            """
            int ix = i / MY;
            int iy = i % MY;
            
            float x = (ix+1)*dx;
            float y = (iy+1)*dy;
            
            
            u[i] = x*((MX+1)*dx - x)*y*((MY+1)*dy - y)*exp(coef*x*y)
            """,
            "initial"
        )
    
    def RhsFn(self, t, u, udot):
        """
        f routine. Compute f(t,u). 
         
        """
        #print t, u.array.reshape(self.mx, self.my)[0,:]
        self.advdiff(self.mx, self.my, self.hdcoef, self.hacoef, 
                     self.vdcoef, u.array, udot.array)
                     
        return 0
                     
    def SetIC(self, u, coef=5.0):
        """
        Load initial profile into u vector 
        """
        self.initial(self.mx, self.my, coef, self.dx, self.dy, u.array)
                

def PrintHeader(reltol, abstol, umax):
    """
    Print first lines of output (problem description) 
    """

    print "\n2-D Advection-Diffusion Equation"
    print "Mesh dimensions = %d X %d" % (MX, MY)
    print "Total system size = %d" % NEQ
    print "Tolerance parameters: reltol = %g   abstol = %g\n" % (reltol, abstol)
    print "At t = %g      max.norm(u) =%14.6e " %(T0, umax)                     
        
if __name__ == '__main__':
    import time
    
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)
    
    XMAX  =2.0    # domain boundaries         
    YMAX  =1.0
    MX    =1000             # mesh dimensions           
    MY    =50
    NEQ   =MX*MY          # number of equations       
    ATOL  =1.0e-5 # scalar absolute tolerance 
    T0    =0.0    # initial time              
     
    DTOUT =0.01    # output time increment     
    T1    =DTOUT    # first output time        
    NOUT  =10             # number of output times    
    
    ZERO =0.0
    HALF =0.5
    ONE  =1.0
    TWO  =2.0
    FIVE =5.0
        
    reltol = ZERO  # Set the tolerances 
    abstol = ATOL    
    
    dx = XMAX/(MX+1);  # Set grid coefficients in data 
    dy = YMAX/(MY+1);

    hdcoef = ONE/(dx*dx)
    hacoef = HALF/(TWO*dx)
    vdcoef = ONE/(dy*dy)
    
    
    cvode_py = AdvDiff(dx, dy, hdcoef, hacoef, vdcoef,
                        multistep='bdf', iteration='newton')
    u_py = NvectorNdarrayFloat64( (MX,MY) )
    
    cvode_vec = AdvDiffVec(dx, dy, hdcoef, hacoef, vdcoef,
                        multistep='bdf', iteration='newton')
    u_vec = NvectorNdarrayFloat64( (MX,MY) )    
    
    u_cl = NvectorCL(ctx, queue, MX*MY,)
    # Create CL instance 
    cvode_cl = AdvDiffCl(ctx, MX, MY, dx, dy, ONE/(dx*dx), HALF/(TWO*dx),  ONE/(dy*dy),
                        multistep='bdf', iteration='newton')
                        
                        
                        
    for cvode_mem, u in ((cvode_vec, u_vec), (cvode_cl, u_cl)):
        t1 = time.time()
        cvode_mem.SetIC(u, FIVE) # Initialize u vector 
    
        # Call CVodeInit to initialize the integrator memory and specify the
        # user's right hand side function in u'=f(t,u), the inital time T0, and
        # the initial dependent variable vector u. 
        cvode_mem.initSolver( T0, u )
        # Call CVodeSStolerances to specify the scalar relative tolerance
        #  and scalar absolute tolerance 
        cvode_mem.setTolerances(reltol, abstol)
        
        # Call CVBand to specify the CVBAND band linear solver 
        # Set the user-supplied Jacobian routine Jac 
        #cvode_mem.setupBandLinearSolver( NEQ, MY, MY, user_jac=True)
        cvode_mem.setupIndirectLinearSolver(solver='spgmr', prectype='none', user_jac=False,
                                            gstype='modified', epslin=0.05)  
                                            
        #udot = u.Clone()
        #cvode_mem.RhsFn(0.0, u, udot)
        #print udot.array.reshape(MX,MY)
        
          
        # In loop over output points: call CVode, print results, test for errors 
        umax = u.MaxNorm()

        PrintHeader(reltol, abstol, umax);
        tout = T1
        for iout in range(NOUT):
            flag, t = cvode_mem.Solve(tout, u)
            umax = u.MaxNorm()
            
            cvode_mem.PrintOutput(t, umax,)
            
            tout += DTOUT
            
        cvode_mem.PrintFinalStats()  # Print some final statistics   
        
        t2 = time.time()
        
        print 'Solving took %0.3f s' % (t2-t1)
    