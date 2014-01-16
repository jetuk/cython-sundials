"""
/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2008/12/17 19:38:48 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example (serial):
 *
 * This example solves a nonlinear system from.
 *
 * Source: "Handbook of Test Problems in Local and Global Optimization",
 *             C.A. Floudas, P.M. Pardalos et al.
 *             Kluwer Academic Publishers, 1999.
 * Test problem 4 from Section 14.1, Chapter 14: Ferraris and Tronconi
 * 
 * This problem involves a blend of trigonometric and exponential terms.
 *    0.5 sin(x1 x2) - 0.25 x2/pi - 0.5 x1 = 0
 *    (1-0.25/pi) ( exp(2 x1)-e ) + e x2 / pi - 2 e x1 = 0
 * such that
 *    0.25 <= x1 <=1.0
 *    1.5 <= x2 <= 2 pi
 * 
 * The treatment of the bound constraints on x1 and x2 is done using
 * the additional variables
 *    l1 = x1 - x1_min >= 0
 *    L1 = x1 - x1_max <= 0
 *    l2 = x2 - x2_min >= 0
 *    L2 = x2 - x2_max >= 0
 * 
 * and using the constraint feature in KINSOL to impose
 *    l1 >= 0    l2 >= 0
 *    L1 <= 0    L2 <= 0
 * 
 * The Ferraris-Tronconi test problem has two known solutions.
 * The nonlinear system is solved by KINSOL using different 
 * combinations of globalization and Jacobian update strategies 
 * and with different initial guesses (leading to one or the other
 * of the known solutions).
 *
 *
 * Constraints are imposed to make all components of the solution
 * positive.
 * -----------------------------------------------------------------
 */
"""
import numpy as np
from pySundials.kinsol import Kinsol
from pySundials.sundials import NvectorNdarrayFloat64

from nose.tools import assert_almost_equal


class FerTron(Kinsol):
    """
    Subclass for the predator-prey system.     
    """
    def __init__(self, lb, ub, PT25=0.25, PT5=0.5, ONE=1.0, TWO=2.0):
        self.lb = lb
        self.ub = ub
        
        self.PT25 = PT25
        self.PT5 = PT5
        self.ONE = ONE
        self.TWO = TWO

    def RhsFn(self, u, f):
        """
        System function for predator-prey system 
        """
        lb = self.lb
        ub = self.ub
        
        x1 = u.data[0]
        x2 = u.data[1]
        l1 = u.data[2]
        L1 = u.data[3]
        l2 = u.data[4]
        L2 = u.data[5]
        
        f.data[0] = self.PT5 * np.sin(x1*x2) - self.PT25 * x2 / np.pi - self.PT5 * x1
        f.data[1] = (self.ONE - self.PT25/np.pi)*(np.exp(self.TWO*x1)-np.e) + np.e*x2/np.pi - self.TWO*np.e*x1
        f.data[2] = l1 - x1 + lb[0]
        f.data[3] = L1 - x1 + ub[0]
        f.data[4] = l2 - x2 + lb[1]
        f.data[5] = L2 - x2 + ub[1]
        
        return 0


class testKinFerTron_dns(object):
    @classmethod
    def setUpAll(self, ):
        """
        Class level setup of constants
        """
        
        # /* Problem Constants */

        self.NVAR = 2
        self.NEQ = 3*self.NVAR
        
        self.FTOL=1.e-5 #  /* function tolerance */
        self.STOL=1.e-5 # /* step tolerance */
        
        self.ZERO=0.0
        self.PT25=0.25
        self.PT5=0.5
        self.ONE=1.0
        self.ONEPT5=1.5
        self.TWO=2.0
        
        # Print out the problem size, solution parameters, initial guess.
        PrintHeader('none', self.FTOL, self.STOL)
        
        self.solution = {1: (0.299449,2.83693),
                         2: (0.5, 3.14159) }
        
        
    def setUp(self, ):     
        
        # temporary local variables for convenience
        NEQ = self.NEQ
        ZERO = self.ZERO
        ONE = self.ONE
         
        self.data = FerTron( [self.PT25,self.ONEPT5], [ONE,self.TWO*np.pi] )
        
        #/* Create serial vectors of length NEQ */
        
        
        self.u = NvectorNdarrayFloat64(NEQ) # temporary for init
        self.s = NvectorNdarrayFloat64(NEQ)
        self.c = NvectorNdarrayFloat64(NEQ)
        
        self.s.Constant(ONE) # /* no scaling */
        
        self.c.data[0] =  ZERO   #/* no constraint on x1 */
        self.c.data[1] =  ZERO   #/* no constraint on x2 */
        self.c.data[2] =  ONE    #/* l1 = x1 - x1_min >= 0 */
        self.c.data[3] = -ONE    #/* L1 = x1 - x1_max <= 0 */
        self.c.data[4] =  ONE    #/* l2 = x2 - x2_min >= 0 */
        self.c.data[5] = -ONE    #/* L2 = x2 - x22_min <= 0 */
          

        # Seutp dense solver and tolerances
        self.data.Setup(self.u, normtol=self.FTOL,scsteptol=self.STOL,dense=True)
        self.data.SetConstraints(self.c)
        
        
    def testSolve1(self, ):
        """
        Generator for all of the tests in this problem
        """
        msets = (1,0)
        glstrs = ('none','linesearch')
        

        u1 = self.SetInitialGuess1(self.data)
        print "\n------------------------------------------"
        print "\nInitial guess on lower bounds"
        print "  [x1,x2] = "
        PrintOutput(u1)
        for mset,glstr in ((m,g) for m in msets for g in glstrs):
            u1.Scale(self.ONE, self.u)
            self.SolveIt(glstr, mset, self.solution[1])
                    
                    
    def testSolve2(self, ):
        """
        Generator for all of the tests in this problem
        """
        msets = (1,0)
        glstrs = ('none','linesearch')
        

        u2 =  self.SetInitialGuess2(self.data)
        print "\n------------------------------------------"
        print "\nInitial guess on in middle of feasible region"
        print "  [x1,x2] = "
        PrintOutput(u2)
        for mset,glstr in ((m,g) for m in msets for g in glstrs):
            u2.Scale(self.ONE, self.u)
            self.SolveIt(glstr, mset, self.solution[2])
                    
        
        
    def SetInitialGuess1(self, data):
        """
        Initial guess 1
        
        There are two known solutions for this problem    
        this init. guess should take us to (0.29945; 2.83693)
        """
        u = NvectorNdarrayFloat64(self.NEQ)
    
        lb = data.lb
        ub = data.ub
    

        x1 = lb[0]
        x2 = lb[1]
    
        u.data[0] = x1;
        u.data[1] = x2;
        u.data[2] = x1 - lb[0]
        u.data[3] = x1 - ub[0]
        u.data[4] = x2 - lb[1]
        u.data[5] = x2 - ub[1]
        
        return u
        
    def SetInitialGuess2(self, data):
        """
        Initial guess 2
        
        There are two known solutions for this problem    
        this init. guess should take us to (0.5; 3.1415926)        
        """        
        u = NvectorNdarrayFloat64(self.NEQ)
    
        lb = data.lb
        ub = data.ub
    

        x1 = self.PT5 * (lb[0] + ub[0])
        x2 = self.PT5 * (lb[1] + ub[1])
    
        u.data[0] = x1
        u.data[1] = x2
        u.data[2] = x1 - lb[0]
        u.data[3] = x1 - ub[0]
        u.data[4] = x2 - lb[1]
        u.data[5] = x2 - ub[1]
        
        return u



    def SolveIt(self, glstr, mset, solution):
    
        if mset==1:
            print "Exact Newton"
        else:
            print "Modified Newton"
    
        if glstr == 'none':
            print ""
        else:
            print " with line search"
            
            
        self.data.SetMaxSetupCalls(mset)    
        self.data.Solve(self.u,self.s,self.s,strategy=glstr,print_level=0)
    

        print "Solution:\n  [x1,x2] = "
        PrintOutput(self.u)
        
        assert_almost_equal(solution[0], self.u.data[0], places=4, msg='x1')
        assert_almost_equal(solution[1], self.u.data[1], places=4, msg='x2')
    
        PrintFinalStats(self.data)
    
        return 0


def PrintHeader(globalstrategy, fnormtol, scsteptol):

    print "\nFerraris and Tronconi test problem"
    print "Tolerance parameters:"
    print "  fnormtol  = {}\n  scsteptol = {}".format(fnormtol, scsteptol)



def PrintOutput(u):
    print " {}  {}".format(u.data[0], u.data[1])


def PrintFinalStats(kmem):
    kmem.PrintFinalStats()


if __name__ == '__main__':
    
    
    t = testKinFerTron_dns()
    t.setUpAll()
    t.setUp()
    t.testSolve1()
