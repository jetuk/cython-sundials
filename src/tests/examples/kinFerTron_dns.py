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

# /* Problem Constants */

NVAR = 2
NEQ = 3*NVAR

FTOL=1.e-5 #  /* function tolerance */
STOL=1.e-5 # /* step tolerance */

ZERO=0.0
PT25=0.25
PT5=0.5
ONE=1.0
ONEPT5=1.5
TWO=2.0

#define PI     RCONST(3.1415926)
#define E      RCONST(2.7182818)

class FerTron(Kinsol):
    def __init__(self, lb, ub):
        self.lb = lb
        self.ub = ub

    
    #/* 
    # * System function for predator-prey system 
    # */

    def RhsFn(self, u, f):
        lb = self.lb
        ub = self.ub
        
        x1 = u.data[0]
        x2 = u.data[1]
        l1 = u.data[2]
        L1 = u.data[3]
        l2 = u.data[4]
        L2 = u.data[5]
        
        f.data[0] = PT5 * np.sin(x1*x2) - PT25 * x2 / np.pi - PT5 * x1
        f.data[1] = (ONE - PT25/np.pi)*(np.exp(TWO*x1)-np.e) + np.e*x2/np.pi - TWO*np.e*x1
        f.data[2] = l1 - x1 + lb[0]
        f.data[3] = L1 - x1 + ub[0]
        f.data[4] = l2 - x2 + lb[1]
        f.data[5] = L2 - x2 + ub[1]
        
        #print 'u',u,u.data
        #print 'f',f,f.data
        
        return 0


#/* Functions Called by the KINSOL Solver */
#static int func(N_Vector u, N_Vector f, void *user_data);
#
#/* Private Helper Functions */
#static void SetInitialGuess1(N_Vector u, UserData data);
#static void SetInitialGuess2(N_Vector u, UserData data);
#static int SolveIt(void *kmem, N_Vector u, N_Vector s, int glstr, int mset);
#static void PrintHeader(int globalstrategy, realtype fnormtol,
#                        realtype scsteptol);
#static void PrintOutput(N_Vector u);
#static void PrintFinalStats(void *kmem);
#static int check_flag(void *flagvalue, char *funcname, int opt);
#
#/*
# *--------------------------------------------------------------------
# * MAIN PROGRAM
# *--------------------------------------------------------------------
# */

def main():

    #UserData data;
    #realtype fnormtol, scsteptol;
    #N_Vector u1, u2, u, s, c;
    #int glstr, mset, flag;
    #void *kmem;
    
    #u1 = u2 = u = NULL;
    #s = c = NULL;
    #kmem = NULL;
    #data = NULL;
    #glstr = KIN_NONE;
      
    #/* User data */
     
    data = FerTron( [PT25,ONEPT5], [ONE,TWO*np.pi] )
    
    #/* Create serial vectors of length NEQ */
    u1 = NvectorNdarrayFloat64(NEQ)
    #if (check_flag((void *)u1, "N_VNew_Serial", 0)) return(1);
    
    u2 = NvectorNdarrayFloat64(NEQ)
    #if (check_flag((void *)u2, "N_VNew_Serial", 0)) return(1);
    
    u = NvectorNdarrayFloat64(NEQ)
    #if (check_flag((void *)u, "N_VNew_Serial", 0)) return(1);
    
    s = NvectorNdarrayFloat64(NEQ)
    #if (check_flag((void *)s, "N_VNew_Serial", 0)) return(1);
    
    c = NvectorNdarrayFloat64(NEQ)
    #if (check_flag((void *)c, "N_VNew_Serial", 0)) return(1);
    
    SetInitialGuess1(u1,data)
    SetInitialGuess2(u2,data)
     
    s.Constant(ONE) # /* no scaling */
    
    c.data[0] =  ZERO   #/* no constraint on x1 */
    c.data[1] =  ZERO   #/* no constraint on x2 */
    c.data[2] =  ONE    #/* l1 = x1 - x1_min >= 0 */
    c.data[3] = -ONE    #/* L1 = x1 - x1_max <= 0 */
    c.data[4] =  ONE    #/* l2 = x2 - x2_min >= 0 */
    c.data[5] = -ONE    #/* L2 = x2 - x22_min <= 0 */
          
    #kmem = KINCreate();
    #if (check_flag((void *)kmem, "KINCreate", 0)) return(1);
    
    #flag = KINSetUserData(kmem, data);
    #if (check_flag(&flag, "KINSetUserData", 1)) return(1);
    #flag = KINSetConstraints(kmem, c);
    #if (check_flag(&flag, "KINSetConstraints", 1)) return(1);
    #flag = KINSetFuncNormTol(kmem, fnormtol);
    #if (check_flag(&flag, "KINSetFuncNormTol", 1)) return(1);
    #flag = KINSetScaledStepTol(kmem, scsteptol);
    #if (check_flag(&flag, "KINSetScaledStepTol", 1)) return(1);
    
    #flag = KINInit(kmem, func, u);
    #if (check_flag(&flag, "KINInit", 1)) return(1);
    data.Setup(u, normtol=FTOL,scsteptol=STOL,dense=True)
    data.SetConstraints(c)
    

    #/* Call KINDense to specify the linear solver */

#flag = KINDense(kmem, NEQ);
#if (check_flag(&flag, "KINDense", 1)) return(1);

    #/* Print out the problem size, solution parameters, initial guess. */
    PrintHeader('none', FTOL, STOL)

    #/* --------------------------- */

    print "\n------------------------------------------"
    print "\nInitial guess on lower bounds"
    print "  [x1,x2] = "
    PrintOutput(u1);
    
    msets = (1,0)
    glstrs = ('none','linesearch')
    
    for mset,glstr in ((m,g) for m in msets for g in glstrs):
        u1.Scale(ONE,u)  
        SolveIt(data, u, s, glstr, mset)

    #/* --------------------------- */

    print "\n------------------------------------------"
    print "\nInitial guess in middle of feasible region"
    print "  [x1,x2] = "
    PrintOutput(u2)

    for mset,glstr in ((m,g) for m in msets for g in glstrs):
        u2.Scale(ONE,u)
        SolveIt(data, u, s, glstr, mset)


    #/* Free memory */

    #N_VDestroy_Serial(u);
    #N_VDestroy_Serial(s);
    #N_VDestroy_Serial(c);
    #KINFree(&kmem);
    #free(data);

    return 0



def SolveIt(kmem, u, s, glstr, mset):

    if mset==1:
        print "Exact Newton"
    else:
        print "Modified Newton"

    if glstr == 'none':
        print ""
    else:
        print " with line search"
        
        
    kmem.SetMaxSetupCalls(mset)

    #flag = KINSetMaxSetupCalls(kmem, mset);
    #if (check_flag(&flag, "KINSetMaxSetupCalls", 1)) return(1);

    kmem.Solve(u,s,s,strategy=glstr,print_level=0)

    #flag = KINSol(kmem, u, glstr, s, s);
    #if (check_flag(&flag, "KINSol", 1)) return(1);

    print "Solution:\n  [x1,x2] = "
    PrintOutput(u)

    PrintFinalStats(kmem)

    return 0



#/*
# *--------------------------------------------------------------------
# * FUNCTIONS CALLED BY KINSOL
# *--------------------------------------------------------------------
# */




#/*
# *--------------------------------------------------------------------
# * PRIVATE FUNCTIONS
# *--------------------------------------------------------------------
# */

#/*
# * Initial guesses
# */

def SetInitialGuess1(u, data):


    lb = data.lb
    ub = data.ub

    #/* There are two known solutions for this problem */

    #/* this init. guess should take us to (0.29945; 2.83693) */
    x1 = lb[0]
    x2 = lb[1]

    u.data[0] = x1;
    u.data[1] = x2;
    u.data[2] = x1 - lb[0]
    u.data[3] = x1 - ub[0]
    u.data[4] = x2 - lb[1]
    u.data[5] = x2 - ub[1]


def SetInitialGuess2(u, data):

    lb = data.lb
    ub = data.ub

    #/* There are two known solutions for this problem */

    #/* this init. guess should take us to (0.5; 3.1415926) */
    x1 = PT5 * (lb[0] + ub[0])
    x2 = PT5 * (lb[1] + ub[1])

    u.data[0] = x1
    u.data[1] = x2
    u.data[2] = x1 - lb[0]
    u.data[3] = x1 - ub[0]
    u.data[4] = x2 - lb[1]
    u.data[5] = x2 - ub[1]


#/* 
# * Print first lines of output (problem description)
# */

def PrintHeader(globalstrategy, fnormtol, scsteptol):

    print "\nFerraris and Tronconi test problem"
    print "Tolerance parameters:"
    print "  fnormtol  = {}\n  scsteptol = {}".format(fnormtol, scsteptol)


#/* 
# * Print solution
# */

def PrintOutput(u):
    print " {}  {}".format(u.data[0], u.data[1])


#/* 
# * Print final statistics contained in iopt 
# */

def PrintFinalStats(kmem):
    kmem.PrintFinalStats()


if __name__ == '__main__':
    main()
