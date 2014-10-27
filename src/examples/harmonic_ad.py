# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 16:44:44 2014

@author: TOML7828
"""

import numpy as np
from pySundials.kinsol import Kinsol
from pySundials.cvode import Cvode
from pySundials.sundials import NvectorNdarrayFloat64
import pandas as pd

ZERO = 0.0
ATOL  =1.0e-5 # scalar absolute tolerance 
XMAX  =2.0    # domain boundaries         
MX    =20   
NEQ = MX
T0    =0.0    # initial time              
T1    =0.1    # first output time         
DTOUT =0.1    # output time increment     
NOUT  =100            # number of output times  

class kinAdvDiff(Kinsol):
    """
    UserData (contains grid constants) 
    
    """
    
    def __init__(self, mx, dx, hdcoef, hacoef, freq, **kwds):
        Kinsol.__init__(self, **kwds)
        
        self.mx = mx
        self.dx = dx
        self.hdcoef = hdcoef
        self.hacoef = hacoef
        self.freq = freq


    def RhsFn(self, u, udot):
        """
        f routine. Compute f(t,u). 
         
        """
        ud = u.data
        # Extract needed constants from data 
        hordc = self.hdcoef
        horac = self.hacoef        
        mx = self.mx
        freq = self.freq
            
        #bval = np.sin(freq*t)
        
        #bval = np.sin(freq*t)
        bval = 0.0 - 1j
        
        ucplx = ud[...,0] + 1j*ud[...,1]
        udotcplx = np.empty_like( ucplx )

        for i in range(mx):
            # Extract u at x_i, y_j and four neighboring points 
            uij = ucplx[i]
            ult = bval if i == 0 else ucplx[i-1]
            urt = ZERO if i == mx-1 else ucplx[i+1]

            # Set diffusion and advection terms and load into udot 

            hdiff = hordc*(ult - 2.0*uij + urt)
            hadv = horac*(urt - ult)
            udotcplx[i] = hdiff + hadv
  
        
        udotcplx -= 1j*freq*ucplx  
        
        print ucplx, udotcplx
        udot.data[:,0] = udotcplx.real
        udot.data[:,1] = udotcplx.imag
        return 0

class AdvDiff(Cvode):
    """
    UserData (contains grid constants) 
    
    """
    
    def __init__(self, mx, dx, hdcoef, hacoef, freq, **kwds):
        Cvode.__init__(self, **kwds)
        
        self.mx = mx
        self.dx = dx
        self.hdcoef = hdcoef
        self.hacoef = hacoef
        self.freq = freq


    def RhsFn(self, t, u, udot):
        """
        f routine. Compute f(t,u). 
         
        """
        ud = u.data
        # Extract needed constants from data 
        hordc = self.hdcoef
        horac = self.hacoef        
        mx = self.mx
        freq = self.freq
            
        bval = np.sin(freq*t)
        
        for i in range(mx):
            # Extract u at x_i, y_j and four neighboring points 
            uij = ud[i]
            ult = bval if i == 0 else ud[i-1]
            urt = ZERO if i == mx-1 else ud[i+1]

            # Set diffusion and advection terms and load into udot 

            hdiff = hordc*(ult - 2.0*uij + urt)
            hadv = horac*(urt - ult)
            udot.data[i] = hdiff + hadv

        return 0 
        
    def PrintOutput(self, t, umax):
        """
        Print current value 
        """
        
        nst = self.numSteps
        print "At t = %4.2f   max.norm(u) =%14.6e   nst = %4ld" % (t, umax, nst)


    def PrintFinalStats(self, ):
        """
        Get and print some final statistics 
        """

        print "\nFinal Statistics:" 
        print "nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld" % (
    	   self.numSteps, self.numRhsEvals, self.numLinSolvSetups, self.numRhsEvals, self.dlsNumJacEvals)
        print "nni = %-6ld ncfn = %-6ld netf = %ld\n " % (
    	   self.numNonlinSolvIters, self.numNonlinSolvConvFails, self.numErrTestFails)



def PrintHeader(reltol, abstol, umax):
    """
    Print first lines of output (problem description) 
    """

    print "\n2-D Advection-Diffusion Equation"
    print "Mesh dimensions = %d" % MX
    print "Total system size = %d" % NEQ
    print "Tolerance parameters: reltol = %g   abstol = %g\n" % (reltol, abstol)
    print "At t = %g      max.norm(u) =%14.6e " %(T0, umax)

if __name__ == '__main__':
    

    
    reltol = ZERO  # Set the tolerances 
    abstol = ATOL
    
    dx = XMAX/(MX+1);  # Set grid coefficients in data 
    
    freq = 10.0/(2.0*np.pi)
    
    
    hfit = kinAdvDiff(MX, dx, 2.0/dx, 0.5/(2.0*dx),  freq,)
    SHP = (MX,2)
    u = NvectorNdarrayFloat64(SHP)
    s = NvectorNdarrayFloat64(SHP)
    s.Constant(1.0) # /* no scaling */
    
    hfit.initSolver(u)
    hfit.setupDenseLinearSolver(MX*2)
    
    hfit.Solve(u,s,s )    
    
    F = u.data[:,0] + 1j*u.data[:,1]
    
    # Create a serial vector 
    u = NvectorNdarrayFloat64( (MX,) )
    

    # Call CVodeCreate to create the solver memory and specify the 
    # Backward Differentiation Formula and the use of a Newton iteration 
    cvode_mem = AdvDiff(MX, dx, 2.0/dx, 0.5/(2.0*dx),  freq,
                        multistep='bdf', iteration='newton')
    
    # Call CVodeInit to initialize the integrator memory and specify the
    # user's right hand side function in u'=f(t,u), the inital time T0, and
    # the initial dependent variable vector u. 
    cvode_mem.initSolver( T0, u )
    # Call CVodeSStolerances to specify the scalar relative tolerance
    #  and scalar absolute tolerance 
    cvode_mem.setTolerances(reltol, abstol)
    
    # Call CVBand to specify the CVBAND band linear solver 
    # Set the user-supplied Jacobian routine Jac 
    cvode_mem.setupIndirectLinearSolver()
    #cvode_mem.setupBandLinearSolver( NEQ, MY, MY, user_jac=True)
  

    from matplotlib import pyplot
    import matplotlib.animation as animation

    # In loop over output points: call CVode, print results, test for errors 
    umax = u.MaxNorm()

    PrintHeader(reltol, abstol, umax);
    tout = T1
    
    x = np.arange(MX)*dx            # x-array
    fig = pyplot.figure()
    pyplot.grid()
    
    ims = []
    u_all = []
    t_all = [tout,]
    
    for iout in range(NOUT):
        flag, t = cvode_mem.Solve(tout, u)
        umax = u.MaxNorm()
        
        cvode_mem.PrintOutput(t, umax,)
        
        lines1, = pyplot.plot(x,u.data,'b-o')
        lines2, = pyplot.plot(x,(F*np.exp(1j*freq*tout)).real, 'r-o')
        
        print lines1
        
        ims.append([lines1,lines2])
        u_all.append(  u.data.copy() )
        t_all.append( tout )
        
        tout += DTOUT
        
    cvode_mem.PrintFinalStats()  # Print some final statistics           
    
    im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000,
    blit=True)    
    
    u_all = np.array(u_all)
    t_all = np.array(t_all)
    
    u_freq = F*np.exp(1j*freq*t)    
    
    fig2 = pyplot.figure()
    
    pyplot.plot(u_all, 'o-')
    pyplot.grid()
    
    pyplot.show()