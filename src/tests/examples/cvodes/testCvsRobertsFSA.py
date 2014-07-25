# -*- coding: utf-8 -*-
"""
#
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2010/12/01 22:51:32 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 * 
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODE. The problem is from
 * chemical kinetics, and consists of the following three rate
 * equations:         
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
 *    dy3/dt = 3.e7*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = 1e-4 or at which
 * y3 = 0.01. This program solves the problem with the BDF method,
 * Newton iteration with the CVDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance. Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------
 
"""

from pySundials.cvodes import Cvodes
from pySundials.sundials import NvectorNdarrayFloat64
import numpy as np

# Problem Constants 

NEQ   = 3                # number of equations  
NS = 3
Y1    = 1.0             # initial y components 
Y2    = 0.0
Y3    = 0.0
RTOL  = 1.0e-4   # scalar relative tolerance            
ATOL1 = 1.0e-8   # vector absolute tolerance components 
ATOL2 = 1.0e-14
ATOL3 = 1.0e-6
T0    = 0.0      # initial time           
T1    = 0.4      # first output time      
TMULT = 10.0     # output time factor     
NOUT  = 12               # number of output times 


class Roberts(Cvodes):
    
    def __init__(self, p, **kwds):
        Cvodes.__init__(self, **kwds)
        
        self.p = p
    
    def RhsFn(self, t, y, ydot):
        """
         f routine. Compute function f(t,y).        
         
        """
        
        y1 = y.data[0]; y2 = y.data[1]; y3 = y.data[2]
        p1,p2,p3 = self.p
        
        yd1 = ydot.data[0] = -p1*y1 + p2*y2*y3
        yd3 = ydot.data[2] = p3*y2*y2
        ydot.data[1] = -yd1 - yd3
        
        return 0
        
    def RootFn(self, t, y, gout):
        """
        g routine. Compute functions g_i(t,y) for i = 0,1. 
        
        """
        y1 = y.data[0]; y3 = y.data[2]
        
        gout[0] = y1 - 0.0001
        gout[1] = y3 - 0.01
        
        return 0
        
    def DlsDenseJacFn(self, N, t, u, fu, J, tmp1, tmp2, tmp3):
        """
        Jacobian routine. Compute J(t,y) = df/dy
        
        """
        y1 = y.data[0]; y2 = y.data[1]; y3 = y.data[2]
        
        J[0,0] = -0.04
        J[0,1] = 1.0e4*y3          
        J[0,2] = 1.0e4*y2
        J[1,0] = 0.04
        J[1,1] = -1.0e4*y3 - 6.0e7*y2
        J[1,2] = -1.0e4*y2
        J[2,1] = 6.0e7*y2

        return 0        
        
    def SensRhsFn1(self, Ns, t, y, ydot, iS, yS, ySdot, tmp1, tmp2):
        
        y1 = y.data[0]; y2 = y.data[1]; y3 = y.data[2]
        s1 = yS.data[0]; s2 = yS.data[1]; s3 = yS.data[2]
        p1,p2,p3 = self.p

        sd1 = -p1*s1 + p2*y3*s2 + p2*y2*s3;
        sd3 = 2*p3*y2*s2;
        sd2 = -sd1-sd3;
        
        if iS == 0:
            sd1 += -y1
            sd2 +=  y1
        elif iS == 1:    
            sd1 +=  y2*y3
            sd2 += -y2*y3
        elif iS == 2:
            sd2 += -y2*y2
            sd3 +=  y2*y2
    
        ySdot.data[0] = sd1
        ySdot.data[1] = sd2
        ySdot.data[2] = sd3
        
        return 0


 
    def ErrWeightFn(self, y, w,):         
        atol = [ATOL1, ATOL2, ATOL3]
     
        for i in range(3):
           yy = y.data[i]
           ww = RTOL * abs(yy) + atol[i]
           if ww <= 0.0: return -1
           w.data[i] = 1.0/ww
         
        return 0


 
# Private helper functions
  

def PrintOutput( t, y1, y2, y3):
    print "At t = %0.4e      y =%14.6e  %14.6e  %14.6e" % ( t, y1, y2, y3)


def PrintOutputS(uS):
    
    for i in range(len(uS)):
        sdata = uS[i].data
        print "                  Sensitivity %d  %12.4e %12.4e %12.4e" %(i+1, sdata[0], sdata[1], sdata[2]);



def PrintRootInfo( root_f1, root_f2):
    print "    rootsfound[] = %3d %3d" % (root_f1, root_f2)


def PrintFinalStats(cv):
    """
 * Get and print some final statistics
    """
    print "\nFinal Statistics:"
    print "nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld"  % (
	   cv.numSteps, cv.numRhsEvals, cv.numLinSolvSetups,
         cv.dlsNumRhsEvals, cv.dlsNumJacEvals)
    print "nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n" % (
	   cv.numNonlinSolvIters, cv.numNonlinSolvConvFails, 
         cv.numErrTestFails, cv.numGEvals)



if __name__ == '__main__':
    # Main Program
    # Create serial vector of length NEQ for I.C. and abstol 

    y = NvectorNdarrayFloat64(NEQ)
    abstol = NvectorNdarrayFloat64(NEQ)
  

    # Initialize y 
    y.data[:] = [Y1,Y2,Y3]
    
    p = np.array([0.04, 1.0e4, 3.0e7])
    # Set the scalar relative tolerance 
    reltol = RTOL;
    # Set the vector absolute tolerance 
    abstol.data[:] = [ATOL1,ATOL2,ATOL3]
    
  
    # Call CVodeCreate to create the solver memory and specify the 
    # Backward Differentiation Formula and the use of a Newton iteration 
    # Call CVodeRootInit to specify the root function g with 2 components   
    cvode_mem = Roberts(p, multistep='bdf', iteration='newton', nrtfn=2)
  
    # Call CVodeInit to initialize the integrator memory and specify the
    # user's right hand side function in y'=f(t,y), the inital time T0, and
    # the initial dependent variable vector y. 
    cvode_mem.initSolver(T0, y)


    # Call CVodeSVtolerances to specify the scalar relative tolerance
    # and vector absolute tolerances 
    cvode_mem.setTolerances(None, None)

    # Call CVDense to specify the CVDENSE dense linear solver   
    # Set the Jacobian routine to Jac (user-supplied) 
    cvode_mem.setupDenseLinearSolver( NEQ, user_jac=True )
  

    yS = [NvectorNdarrayFloat64(NEQ) for i in range(NS)]
    cvode_mem.initSens(NS, 'simultaneous', yS, 'one')
    
    cvode_mem.setSensTolerances()
    cvode_mem.sensErrCon(True)
    
    
    cvode_mem.sensParams(None, p, None)

    # In loop, call CVode, print results, and test for error.
    # Break out of loop when NOUT preset output times have been reached.  
    print " \n3-species kinetics problem\n\n"
    

    iout = 0;  tout = T1;
    while True:
        flag, t = cvode_mem.Solve(tout, y)        
        PrintOutput(t, y.data[0], y.data[1], y.data[2] )
        
        if flag == 'CV_ROOT_RETURN':
            rootsfound = cvode_mem.rootInfo()
            PrintRootInfo(rootsfound[0], rootsfound[1])
            
        elif flag == 'CV_SUCCESS':
            iout += 1
        else:
            break

        flag, t = cvode_mem.SolveSens( yS )

        PrintOutputS( yS )
        
        print "-----------------------------------------"
        print "------------------------------"        
        
        if iout == NOUT: break
        tout *= TMULT
        
    # Print some final statistics 
    PrintFinalStats(cvode_mem)        


