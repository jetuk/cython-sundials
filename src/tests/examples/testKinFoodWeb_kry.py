"""
/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2010/12/01 23:08:49 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example (serial):
 *
 * This example solves a nonlinear system that arises from a system
 * of partial differential equations. The PDE system is a food web
 * population model, with predator-prey interaction and diffusion
 * on the unit square in two dimensions. The dependent variable
 * vector is the following:
 * 
 *       1   2         ns
 * c = (c , c ,  ..., c  )     (denoted by the variable cc)
 * 
 * and the PDE's are as follows:
 *
 *                    i       i
 *         0 = d(i)*(c     + c    )  +  f  (x,y,c)   (i=1,...,ns)
 *                    xx      yy         i
 *
 *   where
 *
 *                   i             ns         j
 *   f  (x,y,c)  =  c  * (b(i)  + sum a(i,j)*c )
 *    i                           j=1
 *
 * The number of species is ns = 2 * np, with the first np being
 * prey and the last np being predators. The number np is both the
 * number of prey and predator species. The coefficients a(i,j),
 * b(i), d(i) are:
 *
 *   a(i,i) = -AA   (all i)
 *   a(i,j) = -GG   (i <= np , j >  np)
 *   a(i,j) =  EE   (i >  np,  j <= np)
 *   b(i) = BB * (1 + alpha * x * y)   (i <= np)
 *   b(i) =-BB * (1 + alpha * x * y)   (i >  np)
 *   d(i) = DPREY   (i <= np)
 *   d(i) = DPRED   ( i > np)
 *
 * The various scalar parameters are set using define's or in
 * routine InitUserData.
 *
 * The boundary conditions are: normal derivative = 0, and the
 * initial guess is constant in x and y, but the final solution
 * is not.
 *
 * The PDEs are discretized by central differencing on an MX by
 * MY mesh.
 * 
 * The nonlinear system is solved by KINSOL using the method
 * specified in local variable globalstrat.
 *
 * The preconditioner matrix is a block-diagonal matrix based on
 * the partial derivatives of the interaction terms f only.
 *
 * Constraints are imposed to make all components of the solution
 * positive.
 * -----------------------------------------------------------------
 * References:
 *
 * 1. Peter N. Brown and Youcef Saad,
 *    Hybrid Krylov Methods for Nonlinear Systems of Equations
 *    LLNL report UCRL-97645, November 1987.
 *
 * 2. Peter N. Brown and Alan C. Hindmarsh,
 *    Reduced Storage Matrix Methods in Stiff ODE systems,
 *    Lawrence Livermore National Laboratory Report  UCRL-95088,
 *    Rev. 1, June 1987, and  Journal of Applied Mathematics and
 *    Computation, Vol. 31 (May 1989), pp. 40-91. (Presents a
 *    description of the time-dependent version of this test
 *    problem.)
 * -----------------------------------------------------------------
 */
"""
import numpy as np
from pySundials.kinsol import Kinsol
from pySundials.sundials import NvectorNdarrayFloat64, DenseGETRF

#/* Problem Constants */

NUM_SPECIES = 6#  /* must equal 2*(number of prey or predators)
               #               number of prey = number of predators       */ 

#define PI       RCONST(3.1415926535898)   /* pi */ 

MX          =8              #/* MX = number of x mesh points */
MY          =8              #/* MY = number of y mesh points */
NSMX        =NUM_SPECIES * MX
NEQ         =NSMX * MY    #/* number of equations in the system */
AA          =1.0          #/* value of coefficient AA in above eqns */
EE          =10000.       #/* value of coefficient EE in above eqns */
GG          =0.5e-6       #/* value of coefficient GG in above eqns */
BB          =1.0          #/* value of coefficient BB in above eqns */
DPREY       =1.0          #/* value of coefficient dprey above */
DPRED       =0.5          #/* value of coefficient dpred above */
ALPHA       =1.0          #/* value of coefficient alpha above */
AX          =1.0          #/* total range of x variable */
AY          =1.0          #/* total range of y variable */
FTOL        =1.e-7        #/* ftol tolerance */
STOL        =1.e-13       #/* stol tolerance */
THOUSAND    =1000.0       #/* one thousand */
ZERO        =0.           #/* 0. */
ONE         =1.0          #/* 1. */
TWO         =2.0          #/* 2. */
PREYIN      =1.0          #/* initial guess for prey concentrations. */
PREDIN      =30000.0 #/* initial guess for predator concs.      */

#/* User-defined vector access macro: IJ_Vptr */

#/* IJ_Vptr is defined in order to translate from the underlying 3D structure
#   of the dependent variable vector to the 1D storage scheme for an N-vector.
#   IJ_Vptr(vv,i,j) returns a pointer to the location in vv corresponding to 
#   indices is = 0, jx = i, jy = j.    */

#define IJ_Vptr(vv,i,j)   (&NV_Ith_S(vv, i*NUM_SPECIES + j*NSMX))

#/* Type : UserData 
#   contains preconditioner blocks, pivot arrays, and problem constants */

class FoodWeb(Kinsol):
    def __init__(self, mx, my, num_species, ax, ay, uround ):
        self.mx = mx
        self.my = my
        self.ns = num_species
        self.np = num_species/2
        self.ax = ax
        self.ay = ay
        self.dx = ax/(mx-1)
        self.dy = ay/(my-1)
        self.uround = uround
        self.sqruround = np.sqrt(uround)
        
        #/*
        # * Allocate memory for data structure of type UserData 
        # */        
        self.P = np.empty( (mx,my,num_species,num_species), dtype=np.float64 )
        self.pivot = np.empty( (mx,my,num_species), dtype=np.int64 )
        self.acoef = np.empty( (num_species,num_species) )
        self.bcoef = np.empty( num_species )
        self.cox = np.empty( num_species )
        self.coy = np.empty( num_species )
        self.rates = NvectorNdarrayFloat64((mx,my,num_species))
                
        # /* Set up the coefficients a and b plus others found in the equations */
        
        for i in range(self.np):
            a1 = (i,self.np)
            a2 = (i+np,0)
            
            for j in range(self.np):
                self.acoef[i,self.np+j]         = -GG
                self.acoef[i+self.np,j]         =  EE       
                self.acoef[i,j]                 = ZERO
                self.acoef[i+self.np,self.np+j] = ZERO
                
            self.acoef[i,i] = -AA
            self.acoef[i+self.np,i+self.np] = - AA
            
            self.bcoef[i] = BB
            self.bcoef[i+self.np] = -BB
            
            self.cox[i] = DPREY/dx2
            self.cox[i+self.np] = DPRED/dx2
            
            self.coy[i] = DPREY/dy2
            self.coy[i+self.np] = DPRED/dy2
            
    def SetInitialProfiles(self, cc, sc):
        
        ctemp = np.empty(self.ns)
        stemp = np.empty(self.ns)
        
        ctemp[:self.np] = PREYIN
        stemp[:self.np] = ONE
        
        ctemp[self.np:] = PREDIN
        stemp[self.np:] = 0.00001

        for i in range(self.ns):
            cc.data[:,:,i] = ctemp[i]
            sc.data[:,:,i] = stemp[i]
            
            
    def WebRate(xx, yy, jx, jy, cc, rates ):        

        for i in range(self.ns):
            rates[ix,iy,i] = np.dot(cc[ix,iy,:], self.acoef[i,:])
            
        fac = ONE + ALPHA * xx * yy
        
        for i in range(self.ns):
            rates[ix,iy,i] = cc[ix,iy,i] * (self.bcoef[i]*fac+rates[ix,iy,i])
     
            
    def RhsFn(self, cc, fval):
        """        
        /* 
         * System function for predator-prey system 
         */
        """
        delx = self.dx
        dely = self.dy

  
        #/* Loop over all mesh points, evaluating rate array at each point*/
        for jy in range(self.my):
            yy = dely*jy

            #/* Set lower/upper index shifts, special at boundaries. */
            idyl = 1 if jy != 0 else -1
            idyu = 1 if jy != (MY-1) else -1
            
            for jx in range(self.mx):
                xx = delx*jx

                #/* Set left/right index shifts, special at boundaries. */
                idxl = 1 if jx != 0 else -1
                idxr = 1 if jx != (MX-1) else -1

                
                #/* Get species interaction rate array at (xx,yy) */
                self.WebRate(xx, yy, jx, jy, cc, self.rates);
                
                dcyli = cc.data[ix,iy,:] - cc.data[ix,iy-idyl,:]
                dcyui = cc.data[ix,iy+idyu,:] - cc.data[ix,iy,:]
                
                dcxli = cc.data[ix,iy,:] - cc.data[ix-idxl,iy,:]
                dcxri = cc.data[ix+idxr,iy,:] - cc.data[ix,iy,:]
                
                
                #/* Compute the total rate value at (xx,yy) */
                fval.data[ix,iy,:] = self.coy*(dcyui-dcyli)  + self.cox*(dcxri-dcxli) + self.rates[ix,iy,:]
        return 0
        
    def SpilsPrecSetup(self, cc, cscale, fval, fscale,
                                       vtemp1, vtemp2):
        """
        /*
         * Preconditioner setup routine. Generate and preprocess P. 
         */
        """

  #realtype r, r0, uround, sqruround, xx, yy, delx, dely, csave, fac;
  #realtype *cxy, *scxy, **Pxy, *ratesxy, *Pxycol, perturb_rates[NUM_SPECIES];
  #long int i, j, jx, jy, ret;
  #UserData data;
  

        delx = self.dx
        dely = self.dy
  
        uround = self.uround
        sqruround = self.sqruround
        fac = fval.WL2Norm( fscale )
        
        peturb_rates = np.empty_like( self.rates )
        
        r0 = THOUSAND * uround * fac * NEQ
        if r0 == ZERO: r0 = ONE
  
        #/* Loop over spatial points; get size NUM_SPECIES Jacobian block at each */
        for jy in range(self.my):
            yy = jy*dely
            
            for jx in range(self.mx):
                xx = jx*delx
                
      #          Pxy = self.P[ix,iy]
                
      #Pxy = (data->P)[jx][jy];
      #cxy = IJ_Vptr(cc,jx,jy);
      #scxy= IJ_Vptr(cscale,jx,jy);
      #ratesxy = IJ_Vptr((data->rates),jx,jy);
                #/* Compute difference quotients of interaction rate fn. */
                for j in range(self.ns):
                    
                    csave = cc.data[ix,iy,j] # /* Save the j,jx,jy element of cc */
                    r = max(sqruround*abs(csave), r0/cscale.data[ix,iy,j])
                    cc.data[ix,iy,j] += r # /* Perturb the j,jx,jy element of cc */
                    fac = ONE/r
      
        
                    self.WebRate(xx, yy, jx, jy, cc, perturb_rates )
        
                    #/* Restore j,jx,jy element of cc */
                    cc.data[ix,iy,j] = csave
        
                    #/* Load the j-th column of difference quotients */
                    for i in range(self.ns):
                        self.P[ix,iy,j,i] = (perturb_rates[ix,iy,i] - self.rates[ix,iy,i])*fac
            
                #/* Do LU decomposition of size NUM_SPECIES preconditioner block */
            ret = denseGETRF(self.P[ix,iy,:,:], NUM_SPECIES, NUM_SPECIES, self.pivot[jx,jy])
            if ret != 0: return 1
  
  
            return 0
        
        
        def PrecSolveBD(cc, cscale, fval, fscale, vv, ftem):
            """
            /*
             * Preconditioner solve routine 
             */            
            """

          
            for jx in range(self.mx):        

                for jy in range(self.my):

                      #For each (jx,jy), solve a linear system of size NUM_SPECIES.
                      #vxy is the address of the corresponding portion of the vector vv;
                      #Pxy is the address of the corresponding block of the matrix P;
                      #piv is the address of the corresponding block of the array pivot.
                      vxy = vv.data[jx,jy]
                      Pxy = self.P[jx,jy,:,:]
                      piv = self.pivot[jx,jy]
      
                      denseGETRS(Pxy, NUM_SPECIES, piv, vxy);

            return 0



#/* Functions Called by the KINSOL Solver */
#
#static int func(N_Vector cc, N_Vector fval, void *user_data);
#
#static int PrecSetupBD(N_Vector cc, N_Vector cscale,
#                       N_Vector fval, N_Vector fscale,
#                       void *user_data,
#                       N_Vector vtemp1, N_Vector vtemp2);
#
#static int PrecSolveBD(N_Vector cc, N_Vector cscale, 
#                       N_Vector fval, N_Vector fscale, 
#                       N_Vector vv, void *user_data,
#                       N_Vector ftem);
#
#/* Private Helper Functions */
#
#static UserData AllocUserData(void);
#static void InitUserData(UserData data);
#static void FreeUserData(UserData data);
#static void SetInitialProfiles(N_Vector cc, N_Vector sc);
#static void PrintHeader(int globalstrategy, int maxl, int maxlrst, 
#                        realtype fnormtol, realtype scsteptol);
#static void PrintOutput(N_Vector cc);
#static void PrintFinalStats(void *kmem);
#static void WebRate(realtype xx, realtype yy, realtype *cxy, realtype *ratesxy, 
#                    void *user_data);
#static realtype DotProd(long int size, realtype *x1, realtype *x2);
#static int check_flag(void *flagvalue, char *funcname, int opt);

#/*
# *--------------------------------------------------------------------
# * MAIN PROGRAM
# *--------------------------------------------------------------------
# */

def main():

    #int globalstrategy;
    #realtype fnormtol, scsteptol;
    #N_Vector cc, sc, constraints;
    #UserData data;
    #int flag, maxl, maxlrst;
    #void *kmem;

        
    data = self.FoodWeb(MX,MY, NUM_SPECIES, AX,AY, 1E-10)

    #/* Create serial vectors of length NEQ */
    cc = NvectorNdarrayFloat64((MX,MY, NUM_SPECIES))    
    sc = NvectorNdarrayFloat64((MX,MY, NUM_SPECIES))
    

    constraints = NvectorNdarrayFloat64((MX,MY, NUM_SPECIES))
    constraints.Constant(TWO)
  

    data.SetInitialProfiles(cc, sc);

    fnormtol=FTOL; scsteptol=STOL;

    #/* Call KINCreate/KINInit to initialize KINSOL.
       #A pointer to KINSOL problem memory is returned and stored in kmem. */
  
  

    #/* Vector cc passed as template vector. */
    data.SetConstraints(constraints)

    #/* Call KINSpgmr to specify the linear solver KINSPGMR with preconditioner
       #routines PrecSetupBD and PrecSolveBD. */
    maxl = 15
    maxlrst = 2
    data.Setup(u, normtol=FTOL,scsteptol=STOL,dense=False,
               linear_solv='spgmr',user_Presolve=True,
               max_restarts=maxrst)

    #/* Print out the problem size, solution parameters, initial guess. */
    #PrintHeader(globalstrategy, maxl, maxlrst, fnormtol, scsteptol);

    #/* Call KINSol and print output concentration profile */   
    flag = data.Solve(cc,sc,sc)
#  flag = KINSol(kmem,           /* KINSol memory block */
#                cc,             /* initial guess on input; solution vector */
#                globalstrategy, /* global stragegy choice */
#                sc,             /* scaling vector, for the variable cc */
#                sc);            /* scaling vector for function values fval */


    print "\n\nComputed equilibrium species concentrations:\n"
    PrintOutput(cc);

    #/* Print final statistics and free memory */  
    data.PrintFinalStats()

    return 0






def PrintHeader(globalstrategy, maxl, maxlrst, fnormtol, scsteptol):
    """
    /* 
     * Print first lines of output (problem description)
     */    
    """
    print "\nPredator-prey test problem --  KINSol (serial version)\n"
    print "Mesh dimensions = %d X %d"%(MX, MY)
    print "Number of species = %d"% NUM_SPECIES
    print "Total system size = %d\n"% NEQ
    print "Flag globalstrategy = %d (0 = None, 1 = Linesearch)"%globalstrategy
    print "Linear solver is SPGMR with maxl = %d, maxlrst = %d"%(maxl, maxlrst)
    print "Preconditioning uses interaction-only block-diagonal matrix"
    print "Positivity constraints imposed on all components"

    print "Tolerance parameters:  fnormtol = %g   scsteptol = %g"%(fnormtol, scsteptol)
    print "\nInitial profile of concentration"
    print "At all mesh points:  %g %g %g   %g %g %g"%(PREYIN, PREYIN, PREYIN,
                                                      PREDIN, PREDIN, PREDIN)





def PrintOutput(cc):
    """
    /* 
     * Print sampled values of current cc 
     */    
    """
    
    
    jy = 0
    jx = 0
    print "\nAt bottom left:"

    #/* Print out lines with up to 6 values per line */
    for _is in range(NUM_SPECIES):
        if (_is%6)*6 == _is: print ""
        print " %g"% cc.data[jx,jy,_is]

  
    jy = MY-1
    jx = MX-1
    print "\n\nAt top right:"

    #/* Print out lines with up to 6 values per line */
    for _is in range(NUM_SPECIES):
        if (_is%6)*6 == _is: print ""
        print " %g"% cc.data[jx,jy,_is]
  
    print "\n"


/* 
 * Print final statistics contained in iopt 
 */

static void PrintFinalStats(void *kmem)
{
  long int nni, nfe, nli, npe, nps, ncfl, nfeSG;
  int flag;
  
  flag = KINGetNumNonlinSolvIters(kmem, &nni);
  check_flag(&flag, "KINGetNumNonlinSolvIters", 1);
  flag = KINGetNumFuncEvals(kmem, &nfe);
  check_flag(&flag, "KINGetNumFuncEvals", 1);
  flag = KINSpilsGetNumLinIters(kmem, &nli);
  check_flag(&flag, "KINSpilsGetNumLinIters", 1);
  flag = KINSpilsGetNumPrecEvals(kmem, &npe);
  check_flag(&flag, "KINSpilsGetNumPrecEvals", 1);
  flag = KINSpilsGetNumPrecSolves(kmem, &nps);
  check_flag(&flag, "KINSpilsGetNumPrecSolves", 1);
  flag = KINSpilsGetNumConvFails(kmem, &ncfl);
  check_flag(&flag, "KINSpilsGetNumConvFails", 1);
  flag = KINSpilsGetNumFuncEvals(kmem, &nfeSG);
  check_flag(&flag, "KINSpilsGetNumFuncEvals", 1);

  printf("Final Statistics.. \n");
  printf("nni    = %5ld    nli   = %5ld\n", nni, nli);
  printf("nfe    = %5ld    nfeSG = %5ld\n", nfe, nfeSG);
  printf("nps    = %5ld    npe   = %5ld     ncfl  = %5ld\n", nps, npe, ncfl);
  
}

/*
 * Check function return value...
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns a flag so check if
 *             flag >= 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer 
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); 
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  return(0);
}
