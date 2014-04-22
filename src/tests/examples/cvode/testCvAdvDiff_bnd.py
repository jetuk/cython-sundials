# -*- coding: utf-8 -*-
"""
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2010/12/01 22:51:32 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem with a banded Jacobian,
 * with the program for its solution by CVODE.
 * The problem is the semi-discrete form of the advection-diffusion
 * equation in 2-D:
 *   du/dt = d^2 u / dx^2 + .5 du/dx + d^2 u / dy^2
 * on the rectangle 0 <= x <= 2, 0 <= y <= 1, and the time
 * interval 0 <= t <= 1. Homogeneous Dirichlet boundary conditions
 * are posed, and the initial condition is
 *   u(x,y,t=0) = x(2-x)y(1-y)exp(5xy).
 * The PDE is discretized on a uniform MX+2 by MY+2 grid with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX*MY.
 * This program solves the problem with the BDF method, Newton
 * iteration with the CVBAND band linear solver, and a user-supplied
 * Jacobian routine.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .1, .2, ..., 1.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------
"""

from pySundials.cvode import Cvode
from pySundials.sundials import NvectorNdarrayFloat64
import numpy as np

# Problem Constants 

XMAX  =2.0    # domain boundaries         
YMAX  =1.0
MX    =10             # mesh dimensions           
MY    =5
NEQ   =MX*MY          # number of equations       
ATOL  =1.0e-5 # scalar absolute tolerance 
T0    =0.0    # initial time              
T1    =0.1    # first output time         
DTOUT =0.1    # output time increment     
NOUT  =10             # number of output times    

ZERO =0.0
HALF =0.5
ONE  =1.0
TWO  =2.0
FIVE =5.0

# User-defined vector access macro IJth 

"""
 IJth is defined in order to isolate the translation from the
   mathematical 2-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. 
   IJth(vdata,i,j) references the element in the vdata array for
   u at mesh point (i,j), where 1 <= i <= MX, 1 <= j <= MY.
   The vdata array is obtained via the macro call vdata = NV_DATA_S(v),
   where v is an N_Vector. 
   The variables are ordered by the y index j, then by the x index i. 
"""

class AdvDiff(Cvode):
    """
    UserData (contains grid constants) 
    
    """
    
    def __init__(self, dx, dy, hdcoef, hacoef, vdcoef, **kwds):
        Cvode.__init__(self, **kwds)
        
        self.dx = dx
        self.dy = dy
        self.hdcoef = hdcoef
        self.hacoef = hacoef
        self.vdcoef = vdcoef
    
    def RhsFn(self, t, u, udot):
        """
        f routine. Compute f(t,u). 
         
        """
        # Extract needed constants from data 
        hordc = self.hdcoef
        horac = self.hacoef
        verdc = self.vdcoef
        
            
        # Loop over all grid points. 
        for j in range(MY):
            
            for i in range(MX):
                # Extract u at x_i, y_j and four neighboring points 
                uij = u.data[i,j]
                
                udn = ZERO if j == 0 else u.data[i,j-1]
                uup = ZERO if j == MY-1 else u.data[i,j+1]
                ult = ZERO if i == 0 else u.data[i-1,j]
                urt = ZERO if i == MX-1 else u.data[i+1,j]
                

                # Set diffusion and advection terms and load into udot 

                hdiff = hordc*(ult - TWO*uij + urt)
                hadv = horac*(urt - ult)
                vdiff = verdc*(uup - TWO*uij + udn)
                udot.data[i,j] = hdiff + hadv + vdiff
  

        return 0 

        
    def RootFn(self, t, y, gout):
        """
        g routine. Compute functions g_i(t,y) for i = 0,1. 
        
        """
        y1 = y.data[0]; y3 = y.data[2]
        
        gout[0] = y1 - 0.0001
        gout[1] = y3 - 0.01
        
        return 0
        
    def DlsBandJacFn(self, N, mupper, mlower, t, y, fy, J, tmp1, tmp2, tmp3):     
        """
        Jacobian routine. Compute J(t,u). 
        
        The components of f = udot that depend on u(i,j) are
        f(i,j), f(i-1,j), f(i+1,j), f(i,j-1), f(i,j+1), with
          df(i,j)/du(i,j) = -2 (1/dx^2 + 1/dy^2)
          df(i-1,j)/du(i,j) = 1/dx^2 + .25/dx  (if i > 1)
          df(i+1,j)/du(i,j) = 1/dx^2 - .25/dx  (if i < MX)
          df(i,j-1)/du(i,j) = 1/dy^2           (if j > 1)
          df(i,j+1)/du(i,j) = 1/dy^2           (if j < MY)
  
        """

        hordc = self.hdcoef
        horac = self.hacoef
        verdc = self.vdcoef
  
        for j in range(MY):
            
            for i in range(MX):
                k = j + i*MY
                
                              
                
                #kthCol = BAND_COL(J,k);
            
                # set the kth column of J 
            
                #BAND_COL_ELEM(kthCol,k,k) = -TWO*(verdc+hordc);
                #if (i != 1)  BAND_COL_ELEM(kthCol,k-MY,k) = hordc + horac;
                #if (i != MX) BAND_COL_ELEM(kthCol,k+MY,k) = hordc - horac;
                #if (j != 1)  BAND_COL_ELEM(kthCol,k-1,k)  = verdc;
                #if (j != MY) BAND_COL_ELEM(kthCol,k+1,k)  = verdc;:
        return 0

        
    def SetIC(self, u):
        """
        Load initial profile into u vector 
        """
        for j in range(MY):
            y = (j+1)*self.dy
            for i in range(MX):
                x = (i+1)*self.dx
                u.data[i,j] = x*(XMAX - x)*y*(YMAX - y)*np.exp(FIVE*x*y)
                
                
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
    print "Mesh dimensions = %d X %d" % (MX, MY)
    print "Total system size = %d" % NEQ
    print "Tolerance parameters: reltol = %g   abstol = %g\n" % (reltol, abstol)
    print "At t = %g      max.norm(u) =%14.6e " %(T0, umax)



if __name__ == '__main__':
    

    
    # Create a serial vector 
    u = NvectorNdarrayFloat64( (MX,MY) )
    
    
    reltol = ZERO  # Set the tolerances 
    abstol = ATOL
    
    dx = XMAX/(MX+1);  # Set grid coefficients in data 
    dy = YMAX/(MY+1);
    # Call CVodeCreate to create the solver memory and specify the 
    # Backward Differentiation Formula and the use of a Newton iteration 
    cvode_mem = AdvDiff(dx, dy, ONE/(dx*dx), HALF/(TWO*dx),  ONE/(dy*dy),
                        multistep='bdf', iteration='newton')
    cvode_mem.SetIC(u) # Initialize u vector 

    # Call CVodeInit to initialize the integrator memory and specify the
    # user's right hand side function in u'=f(t,u), the inital time T0, and
    # the initial dependent variable vector u. 
    cvode_mem.initSolver( T0, u )
    # Call CVodeSStolerances to specify the scalar relative tolerance
    #  and scalar absolute tolerance 
    cvode_mem.setTolerances(reltol, abstol)
    
    # Call CVBand to specify the CVBAND band linear solver 
    # Set the user-supplied Jacobian routine Jac 
    cvode_mem.setupBandLinearSolver( NEQ, MY, MY, user_jac=False)
  
      
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




