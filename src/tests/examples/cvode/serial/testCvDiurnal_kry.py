# -*- coding: utf-8 -*-
#
"""
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2010/12/01 22:51:32 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 *
 * An ODE system is generated from the following 2-species diurnal
 * kinetics advection-diffusion PDE system in 2 space dimensions:
 *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)
 *                 + Ri(c1,c2,t)      for i = 1,2,   where
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
 *   Kv(y) = Kv0*exp(y/5) ,
 * Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
 * vary diurnally. The problem is posed on the square
 *   0 <= x <= 20,    30 <= y <= 50   (all in km),
 * with homogeneous Neumann boundary conditions, and for time t in
 *   0 <= t <= 86400 sec (1 day).
 * The PDE system is treated by central differences on a uniform
 * 10 x 10 mesh, with simple polynomial initial profiles.
 * The problem is solved with CVODE, with the BDF/GMRES
 * method (i.e. using the CVSPGMR linear solver) and the
 * block-diagonal part of the Newton matrix as a left
 * preconditioner. A copy of the block-diagonal part of the
 * Jacobian is saved and conditionally reused within the Precond
 * routine.
 * -----------------------------------------------------------------
"""

from pySundials.cvode import Cvode, denseGETRF, denseGETRS
from pySundials.sundials import NvectorNdarrayFloat64
import numpy as np

# Problem Constants 

ZERO =0.0
ONE  =1.0
TWO  =2.0

NUM_SPECIES  =2                 # number of species         
KH           =4.0e-6    # horizontal diffusivity Kh 
VEL          =0.001     # advection velocity V      
KV0          =1.0e-8    # coefficient in Kv(y)      
Q1           =1.63e-16  # coefficients q1, q2, c3    
Q2           =4.66e-16
C3           =3.7e16
A3           =22.62     # coefficient in expression for q3(t) 
A4           =7.601     # coefficient in expression for q4(t) 
C1_SCALE     =1.0e6     # coefficients in initial profiles    
C2_SCALE     =1.0e12

T0           =ZERO                 # initial time 
NOUT         =12                   # number of output times 
TWOHR        =7200.0       # number of seconds in two hours  
HALFDAY      =4.32e4       # number of seconds in a half day 
PI       =3.1415926535898  # pi  

XMIN         =ZERO                 # grid boundaries in x  
XMAX         =20.0           
YMIN         =30.0         # grid boundaries in y  
YMAX         =50.0
XMID         =10.0         # grid midpoints in x,y           
YMID         =40.0

MX           =10             # MX = number of x mesh points 
MY           =10             # MY = number of y mesh points 
NSMX         =20             # NSMX = NUM_SPECIES*MX 
MM           =(MX*MY)        # MM = MX*MY 

# CVodeInit Constants 

RTOL    =1.0e-5    # scalar relative tolerance 
FLOOR   =100.0     # value of C1 or C2 at which tolerances 
                                  # change from relative to absolute      
ATOL    =(RTOL*FLOOR)      # scalar absolute tolerance 
NEQ     =(NUM_SPECIES*MM)  # NEQ = number of equations 

# User-defined vector and matrix accessor macros: IJKth, IJth 

"""
   IJKth is defined in order to isolate the translation from the
   mathematical 3-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. IJth is defined in order to
   write code which indexes into small dense matrices with a (row,column)
   pair, where 1 <= row, column <= NUM_SPECIES.   
   
   IJKth(vdata,i,j,k) references the element in the vdata array for
   species i at mesh point (j,k), where 1 <= i <= NUM_SPECIES,
   0 <= j <= MX-1, 0 <= k <= MY-1. The vdata array is obtained via
   the macro call vdata = NV_DATA_S(v), where v is an N_Vector. 
   For each mesh point (j,k), the elements for species i and i+1 are
   contiguous within vdata.

   IJth(a,i,j) references the (i,j)th entry of the small matrix realtype **a,
   where 1 <= i,j <= NUM_SPECIES. The small matrix routines in sundials_dense.h
   work with matrices stored by column in a 2-dimensional array. In C,
   arrays are indexed starting at 0, not 1. 
"""   

#define IJKth(vdata,i,j,k) (vdata[i-1 + (j)*NUM_SPECIES + (k)*NSMX])
#define IJth(a,i,j)        (a[j-1][i-1])


class Diurnal(Cvode):
    """
# Type : UserData 
   contains preconditioner blocks, pivot arrays, and problem constants 
    """
    
    def __init__(self, *args, **kwds):
        Cvode.__init__(self,*args,**kwds)
    
        # Allocate memory for data structure of type UserData 
    
        self.P = np.empty( (MX,MY,NUM_SPECIES, NUM_SPECIES) )
        self.Jbd = np.empty( (MX,MY,NUM_SPECIES, NUM_SPECIES) )
        self.pivot = np.empty( (MX,MY,NUM_SPECIES), dtype=np.int64 )


        # Load problem constants in data 

        self.om = PI/HALFDAY
        self.dx = (XMAX-XMIN)/(MX-1)
        self.dy = (YMAX-YMIN)/(MY-1)
        self.hdco = KH/self.dx**2
        self.haco = VEL/(TWO*self.dx)
        self.vdco = (ONE/self.dy**2)*KV0

 
    def SetInitialProfiles(self, u):
        """        
        Set initial conditions in u 
        """

        # Set pointer to data array in vector u. 

        udata = u.data

        # Load initial profiles of c1 and c2 into u vector 

        for jy in range(MY):
  
            y = YMIN + jy*self.dy
            cy = (0.1*(y - YMID))**2
            cy = ONE - cy + 0.5*cy**2
            
            for jx in range(MX):
    
                x = XMIN + jx*self.dx
                cx = (0.1*(x - XMID))**2
                cx = ONE - cx + 0.5*cx**2
                udata[jx,jy,0] = C1_SCALE*cx*cy
                udata[jx,jy,1] = C2_SCALE*cx*cy
                
    def PrintOutput(self, u, t):
        """
        Print current t, step count, order, stepsize, and sampled c1,c2 values         
        """

        mxh = MX/2 - 1; myh = MY/2 - 1; mx1 = MX - 1; my1 = MY - 1;

        udata = u.data
        
        nst = self.numSteps
        qu = self.lastOrder
        hu = self.lastStep

        print "t = %.2e   no. steps = %ld   order = %d   stepsize = %.2e" % (t, nst, qu, hu)
        print "c1 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e" % (udata[0,0,0], udata[mxh,myh,0], udata[mx1,my1,0])
        print "c2 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e\n" % (udata[0,0,1], udata[mxh,myh,1], udata[mx1,my1,1])

    def PrintFinalStats(self, ):
        """
        Get and print final statistics         
        
        """
        
        lenrw, leniw = self.workSpace
        nst = self.numSteps
        nfe = self.numRhsEvals
        nsetups = self.numLinSolvSetups
        netf = self.numErrTestFails
        nni = self.numNonlinSolvIters
        ncfn = self.numNonlinSolvConvFails
        
        lenrwLS, leniwLS = self.spilsWorkSpace
        nli = self.spilsNumLinIters
        npe = self.spilsNumPrecEvals
        nps = self.spilsNumPrecSolves
        ncfl = self.spilsNumConvFails
        nfeLS = self.spilsNumRhsEvals


        print "\nFinal Statistics.. \n"
        print "lenrw   = %5ld     leniw   = %5ld" %(lenrw, leniw)
        print "lenrwLS = %5ld     leniwLS = %5ld" % (lenrwLS, leniwLS)
        print "nst     = %5ld"                  % nst
        print "nfe     = %5ld     nfeLS   = %5ld" % ( nfe, nfeLS)
        print "nni     = %5ld     nli     = %5ld" % (nni, nli);
        print "nsetups = %5ld     netf    = %5ld" % (nsetups, netf)
        print "npe     = %5ld     nps     = %5ld" % (npe, nps)
        print "ncfn    = %5ld     ncfl    = %5ld\n"% (ncfn, ncfl)


     #-------------------------------
     # Functions called by the solver
     #*-------------------------------
   
    def RhsFn(self, t, u, udot):
        """
         f routine. Compute RHS function f(t,u). 
         
        """
        udata = u.data
        dudata = udot.data


        # Set diurnal rate coefficients. 

        s = np.sin(self.om*t)
        if s > ZERO:
            q3 = np.exp(-A3/s)
            self.q4 = np.exp(-A4/s)
        else:
            q3 = ZERO
            self.q4 = ZERO


        # Make local copies of problem variables, for efficiency. 

        q4coef = self.q4
        dely = self.dy
        verdco = self.vdco
        hordco  = self.hdco
        horaco  = self.haco

        # Loop over all grid points. 
        for jy in range(MY):

            # Set vertical diffusion coefficients at jy +- 1/2 

            ydn = YMIN + (jy - 0.5)*dely
            yup = ydn + dely
            cydn = verdco*np.exp(0.2*ydn)
            cyup = verdco*np.exp(0.2*yup)
            idn = 1 if jy == 0 else -1
            iup = -1 if jy == MY-1 else 1
            
            for jx in range(MX):
    
                # Extract c1 and c2, and set kinetic rate terms. 
                c1 = udata[jx,jy,0]
                c2 = udata[jx,jy,1]

                qq1 = Q1*c1*C3
                qq2 = Q2*c1*c2
                qq3 = q3*C3
                qq4 = q4coef*c2
                rkin1 = -qq1 - qq2 + TWO*qq3 + qq4
                rkin2 = qq1 - qq2 - qq4

                # Set vertical diffusion terms. 

                c1dn = udata[jx,jy+idn,0]
                c2dn = udata[jx,jy+idn,1]
                c1up = udata[jx,jy+iup,0]
                c2up = udata[jx,jy+iup,1]
                vertd1 = cyup*(c1up - c1) - cydn*(c1 - c1dn)
                vertd2 = cyup*(c2up - c2) - cydn*(c2 - c2dn)

                # Set horizontal diffusion and advection terms. 
                
                ileft = 1 if jx == 0 else -1
                iright = -1 if jx == MX-1 else 1
                c1lt = udata[jx+ileft,jy,0]
                c2lt = udata[jx+ileft,jy,1]
                c1rt = udata[jx+iright,jy,0]
                c2rt = udata[jx+iright,jy,1]
                hord1 = hordco*(c1rt - TWO*c1 + c1lt)
                hord2 = hordco*(c2rt - TWO*c2 + c2lt)
                horad1 = horaco*(c1rt - c1lt)
                horad2 = horaco*(c2rt - c2lt)

                # Load all terms into udot. 

                dudata[jx, jy, 0] = vertd1 + hord1 + horad1 + rkin1
                dudata[jx, jy, 1] = vertd2 + hord2 + horad2 + rkin2
                
        return 0

    def SpilsJacTimesVec(self, v, Jv, t, u, fu, tmp):
        """
        Jacobian-times-vector routine. 
        """

        udata = u.data
        vdata = v.data
        Jvdata = Jv.data

        # Set diurnal rate coefficients. 

        s = np.sin(self.om*t)
        if s > ZERO:
            self.q4 = np.exp(-A4/s)
        else:
            self.q4 = ZERO


        # Make local copies of problem variables, for efficiency. 

        q4coef = self.q4
        dely = self.dy
        verdco = self.vdco
        hordco  = self.hdco
        horaco  = self.haco

        # Loop over all grid points. 

        for jy in range(MY):
            
            # Set vertical diffusion coefficients at jy +- 1/2 

            ydn = YMIN + (jy - 0.5)*dely
            yup = ydn + dely

            cydn = verdco*np.exp(0.2*ydn)
            cyup = verdco*np.exp(0.2*yup)

            idn = 1 if jy == 0 else -1
            iup = -1 if jy == MY-1 else 1

            for jx in range(MX):

                Jv1 = ZERO
                Jv2 = ZERO

                # Extract c1 and c2 at the current location and at neighbors 

                c1 = udata[jx,jy,0] 
                c2 = udata[jx,jy,1]

                v1 = vdata[jx,jy,0] 
                v2 = vdata[jx,jy,1]

                c1dn = udata[jx,jy+idn,0]
                c2dn = udata[jx,jy+idn,1]
                c1up = udata[jx,jy+iup,0]
                c2up = udata[jx,jy+iup,1]

                v1dn = vdata[jx,jy+idn,0]
                v2dn = vdata[jx,jy+idn,1]
                v1up = vdata[jx,jy+iup,0]
                v2up = vdata[jx,jy+iup,1]

                ileft = 1 if jx == 0 else -1
                iright = -1 if jx == MX-1 else 1

                c1lt = udata[jx+ileft,jy,0] 
                c2lt = udata[jx+ileft,jy,1]
                c1rt = udata[jx+iright,jy,0]
                c2rt = udata[jx+iright,jy,1]

                v1lt = vdata[jx+ileft,jy,0] 
                v2lt = vdata[jx+ileft,jy,1]
                v1rt = vdata[jx+iright,jy,0]
                v2rt = vdata[jx+iright,jy,1]

                # Set kinetic rate terms. 

                #rkin1 = -Q1*C3 * c1 - Q2 * c1*c2 + q4coef * c2  + TWO*C3*q3;
                #rkin2 =  Q1*C3 * c1 - Q2 * c1*c2 - q4coef * c2;

                Jv1 += -(Q1*C3 + Q2*c2) * v1  +  (q4coef - Q2*c1) * v2
                Jv2 +=  (Q1*C3 - Q2*c2) * v1  -  (q4coef + Q2*c1) * v2

                # Set vertical diffusion terms. 

                #vertd1 = -(cyup+cydn) * c1 + cyup * c1up + cydn * c1dn;
                #vertd2 = -(cyup+cydn) * c2 + cyup * c2up + cydn * c2dn;

                Jv1 += -(cyup+cydn) * v1  +  cyup * v1up  +  cydn * v1dn
                Jv2 += -(cyup+cydn) * v2  +  cyup * v2up  +  cydn * v2dn

                # Set horizontal diffusion and advection terms. 

                #hord1 = hordco*(c1rt - TWO*c1 + c1lt);
                #hord2 = hordco*(c2rt - TWO*c2 + c2lt);

                Jv1 += hordco*(v1rt - TWO*v1 + v1lt)
                Jv2 += hordco*(v2rt - TWO*v2 + v2lt)

                #horad1 = horaco*(c1rt - c1lt)
                #horad2 = horaco*(c2rt - c2lt)

                Jv1 += horaco*(v1rt - v1lt)
                Jv2 += horaco*(v2rt - v2lt)

                # Load two components of J*v 

                #IJKth(dudata, 1, jx, jy) = vertd1 + hord1 + horad1 + rkin1; 
                #IJKth(dudata, 2, jx, jy) = vertd2 + hord2 + horad2 + rkin2;

                Jvdata[jx, jy, 0] = Jv1
                Jvdata[jx, jy, 1] = Jv2
                
        return 0
    
    def SpilsPrecSetup(self, tn, u, fu, jok, gamma, vtemp1, vtemp2, vtemp3):
        """
        Preconditioner setup routine. Generate and preprocess P.         
        """
        
        # Make local copies of pointers in user_data, and of pointer to u's data 

        P = self.P
        Jbd = self.Jbd
        pivot = self.pivot  
        udata = u.data
        
        if jok:   
            # jok = TRUE: Copy Jbd to P 
            for jy in range(MY):
                for jx in range(MX):
                    P[jx,jy,...] = Jbd[jx,jy,...]
                    
            jcurPtr = False
        else:
            # jok = FALSE: Generate Jbd from scratch and copy to P 
    
            # Make local copies of problem variables, for efficiency. 
    
            q4coef = self.q4;
            dely = self.dy;
            verdco = self.vdco;
            hordco  = self.hdco;
    
            # Compute 2x2 diagonal Jacobian blocks (using q4 values 
            # computed on the last f call).  Load into P. 
    
            for jy in range(MY):
                ydn = YMIN + (jy - 0.5)*dely;
                yup = ydn + dely;
                cydn = verdco*np.exp(0.2*ydn);
                cyup = verdco*np.exp(0.2*yup);
                diag = -(cydn + cyup + TWO*hordco);
                for jx in range(MX):
                    c1 = udata[jx,jy,0]
                    c2 = udata[jx,jy,1]

                    Jbd[jx,jy,0,0] = (-Q1*C3 - Q2*c2) + diag
                    # Transposed in Python
                    Jbd[jx,jy,1,0] = -Q2*c1 + q4coef
                    Jbd[jx,jy,0,1] = Q1*C3 - Q2*c2
                    Jbd[jx,jy,1,1] = (-Q2*c1 - q4coef) + diag
                    P[jx,jy,...] = Jbd[jx,jy,...]
            jcurPtr = True
            
        # Scale by -gamma 
            
        for jy in range(MY):
            for jx in range(MX):
                P[jx,jy,...] *= -gamma
  
  
        # Add identity matrix and do LU decompositions on blocks in place. 
        for jy in range(MY):
            for jx in range(MX):
                P[jx,jy,...] += np.eye(NUM_SPECIES)
                ier = denseGETRF(P[jx,jy,...], pivot[jx,jy,...] )
               # print ier, P[jx,jy,...]
                if ier != 0: return 1, jcurPtr
                

        return 0, jcurPtr

    def SpilsPrecSolve(self, t, u, fu, r, z, gamma, delta, lr, tmp):
        """
        Preconditioner solve routine 
        """

        # Extract the P and pivot arrays from user_data. 

        P = self.P
        pivot = self.pivot
        zdata = z.data
        
        r.Scale(ONE, z)

        
  
        # Solve the block-diagonal system Px = r using LU factors stored
        # in P and pivot data in pivot, and return the solution in z. 
  
        for jy in range(MY):
            for jx in range(MX):
                #v = 
                #print jx, jy, P[jx,jy,...], pivot[jx,jy,...], zdata[jx,jy,:]
                denseGETRS(P[jx,jy,...], pivot[jx,jy,...], zdata[jx,jy,:])
                #print jx, jy, P[jx,jy,...], pivot[jx,jy,...], zdata[jx,jy,:]
                
        return 0


# Main Program
 
if __name__ == '__main__':
    
    # Call CVodeCreate to create the solver memory and specify the 
    # Backward Differentiation Formula and the use of a Newton iteration 
    
    cvode_mem = Diurnal( multistep='bdf', iteration='newton')

    # Allocate memory, and set problem data, initial values, tolerances  
    u = NvectorNdarrayFloat64((MX,MY,NUM_SPECIES))
    
    cvode_mem.SetInitialProfiles(u,)
    abstol=ATOL 
    reltol=RTOL

    # Call CVodeInit to initialize the integrator memory and specify the
    # user's right hand side function in u'=f(t,u), the inital time T0, and
    # the initial dependent variable vector u. 
    cvode_mem.initSolver(T0, u)
    # Call CVodeSStolerances to specify the scalar relative tolerance
    # and scalar absolute tolerances     
    cvode_mem.setTolerances(reltol, abstol)

    # Call CVSpgmr to specify the linear solver CVSPGMR 
    # with left preconditioning and the maximum Krylov dimension maxl 
    # set the JAcobian-times-vector function 
    # Set modified Gram-Schmidt orthogonalization 
    # Set the preconditioner solve and setup functions 
    cvode_mem.setupIndirectLinearSolver(solver='spgmr', prectype='none', user_jac=False,
                                        gstype='modified')


    # In loop over output points, call CVode, print results, test for error 
    print " \n2-species diurnal advection-diffusion problem\n"
    tout = TWOHR
    
    for iout in range(NOUT):
        ret, t = cvode_mem.Solve(tout, u)
        cvode_mem.PrintOutput(u, t)
        tout += TWOHR


    cvode_mem.PrintFinalStats()





#
 






