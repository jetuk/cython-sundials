# -*- coding: utf-8 -*-
"""
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2010/12/01 22:51:32 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * --------------------------------------------------------------------
 * Demonstration program for CVODE - Krylov linear solver.
 * ODE system from ns-species interaction PDE in 2 dimensions.
 * 
 * This program solves a stiff ODE system that arises from a system
 * of partial differential equations. The PDE system is a food web
 * population model, with predator-prey interaction and diffusion on
 * the unit square in two dimensions. The dependent variable vector is:
 *
 *        1   2        ns
 *  c = (c , c , ..., c  )
 *
 * and the PDEs are as follows:
 *
 *    i               i      i
 *  dc /dt  =  d(i)*(c    + c   )  +  f (x,y,c)  (i=1,...,ns)
 *                    xx     yy        i   
 *
 * where
 *
 *                 i          ns         j
 *  f (x,y,c)  =  c *(b(i) + sum a(i,j)*c )
 *   i                       j=1                                         
 *                                                                       
 * The number of species is ns = 2*np, with the first np being prey
 * and the last np being predators. The coefficients a(i,j), b(i),
 * d(i) are:
 *
 *  a(i,i) = -a  (all i)
 *  a(i,j) = -g  (i <= np, j > np)
 *  a(i,j) =  e  (i > np, j <= np)
 *  b(i) =  b*(1 + alpha*x*y)  (i <= np)
 *  b(i) = -b*(1 + alpha*x*y)  (i > np)
 *  d(i) = Dprey  (i <= np)
 *  d(i) = Dpred  (i > np)
 *
 * The spatial domain is the unit square. The final time is 10.
 * The boundary conditions are: normal derivative = 0.
 * A polynomial in x and y is used to set the initial conditions.
 *
 * The PDEs are discretized by central differencing on an MX by MY mesh.
 *
 * The resulting ODE system is stiff.
 *
 * The ODE system is solved using Newton iteration and the CVSPGMR
 * linear solver (scaled preconditioned GMRES).
 *
 * The preconditioner matrix used is the product of two matrices:
 * (1) A matrix, only defined implicitly, based on a fixed number
 * of Gauss-Seidel iterations using the diffusion terms only.
 * (2) A block-diagonal matrix based on the partial derivatives
 * of the interaction terms f only, using block-grouping (computing
 * only a subset of the ns by ns blocks).
 *
 * Four different runs are made for this problem.
 * The product preconditoner is applied on the left and on the
 * right. In each case, both the modified and classical Gram-Schmidt
 * options are tested.
 * In the series of runs, CVodeInit and CVSpgmr are called only
 * for the first run, whereas CVodeReInit and CVReInitSpgmr are
 * called for each of the remaining three runs.
 *
 * A problem description, performance statistics at selected output
 * times, and final statistics are written to standard output.
 * On the first run, solution values are also printed at output
 * times. Error and warning messages are written to standard error,
 * but there should be no such messages.
 *
 * Note: This program requires the dense linear solver functions
 * newDenseMat, newLintArray, denseAddIdentity, denseGETRF, denseGETRS, 
 * destroyMat and destroyArray.
 *
 * Note: This program assumes the sequential implementation for the
 * type N_Vector and uses the NV_DATA_S macro to gain access to the
 * contiguous array of components of an N_Vector.
 * --------------------------------------------------------------------
 * Reference: Peter N. Brown and Alan C. Hindmarsh, Reduced Storage
 * Matrix Methods in Stiff ODE Systems, J. Appl. Math. & Comp., 31
 * (1989), pp. 40-91.  Also available as Lawrence Livermore National
 * Laboratory Report UCRL-95088, Rev. 1, June 1987.
 * --------------------------------------------------------------------
"""

from pySundials.cvode import Cvode, denseAddIdentity, denseGETRF, denseGETRS
from pySundials.sundials import NvectorNdarrayFloat64
import numpy as np

# Constants 

ZERO =0.0
ONE  =1.0

# Problem Specification Constants 

AA    =ONE               # AA = a 
EE    =1.0e4     # EE = e 
GG    =0.5e-6    # GG = g 
BB    =ONE               # BB = b 
DPREY =ONE
DPRED =0.5
ALPH  =ONE
NP    =1
NS    =(2*NP)

# Method Constants 

MX    =6
MY    =6
MXNS  =(MX*NS)
AX    =ONE
AY    =ONE
DX    =(AX/float(MX-1))
DY    =(AY/float(MY-1))
MP    =NS
MQ    =(MX*MY)
MXMP  =(MX*MP)
NGX   =2
NGY   =2
NGRP  =(NGX*NGY)
ITMAX =5

# CVodeInit Constants 

NEQ  =(NS*MX*MY)
SHP = (NS,MX,MY)
T0   =ZERO
RTOL =1.0e-5
ATOL =1.0e-5

# CVSpgmr Constants 

MAXL =0     # => use default = MIN(NEQ, 5)            
DELT =ZERO  # => use default = 0.05                   

# Output Constants 

T1        =1.0e-8
TOUT_MULT =10.0
DTOUT     =ONE
NOUT      =18

UNIT_ROUNDOFF = 1e-16

# Note: The value for species i at mesh point (j,k) is stored in 
# component number (i-1) + j*NS + k*NS*MX of an N_Vector,        
# where 1 <= i <= NS, 0 <= j < MX, 0 <= k < MY.                  

# Structure for user data 

class WebData(Cvode):
    """
    
    """
    def __init__(self, *args, **kwds):
        """
        Initialise User Data
        """
        Cvode.__init__(self, **kwds)
        
        ns = self.ns = NS
        
        self.P = np.empty( (NGX, NGY, NS,NS) )
        self.pivot = np.empty( NS, dtype=np.int64 )
        
        self.rewt = NvectorNdarrayFloat64( SHP )
        
        self.jgx = np.empty( NGX+1, dtype=np.int64 )
        self.jgy = np.empty( NGY+1, dtype=np.int64 )
        self.jigx = np.empty( MX, dtype=np.int64 )
        self.jigy = np.empty( MY, dtype=np.int64 )
        
        self.jxr = np.empty( NGX, dtype=np.int64 )
        self.jyr = np.empty( NGY, dtype=np.int64 )
        
        self.acoef = np.empty( (NS,NS) )
        self.bcoef = np.empty( NS )
        self.diff = np.empty( NS )
        self.cox = np.empty( NS ); self.coy = np.empty( NS )
        self.fsave = np.empty( SHP )
        
        for j in range(NS):
            for i in range(NS):
                self.acoef[i,j] = 0.
                
        for j in range(NP):
            for i in range(NP):
                self.acoef[NP+i][j] = EE
                self.acoef[i][NP+j] = -GG
            self.acoef[j,j] = -AA
            self.acoef[NP+j,NP+j] = -AA
            self.bcoef[j] = BB
            self.bcoef[NP+j] = -BB
            self.diff[j] = DPREY
            self.diff[NP+j] = DPRED
        
        # Set remaining problem parameters 
        
        self.mxns = MXNS
        dx = self.dx = DX
        dy = self.dy = DY
        for i in range(ns):
            self.cox[i] = self.diff[i]/np.sqrt(dx)
            self.coy[i] = self.diff[i]/np.sqrt(dy)
          

        # Set remaining method parameters 

        self.mp = MP
        self.mq = MQ
        self.mx = MX
        self.my = MY
        self.srur = np.sqrt(UNIT_ROUNDOFF)
        self.mxmp = MXMP
        self.ngrp = NGRP
        self.ngx = NGX
        self.ngy = NGY
        self.SetGroups(MX, NGX, self.jgx, self.jigx, self.jxr)
        self.SetGroups(MY, NGY, self.jgy, self.jigy, self.jyr)
        
    def SetGroups(self, m, ng, jg, jig, jr):
        """
         This routine sets arrays jg, jig, and jr describing
         a uniform partition of (0,1,2,...,m-1) into ng groups.
         The arrays set are:
           jg    = length ng+1 array of group boundaries.
                   Group ig has indices j = jg[ig],...,jg[ig+1]-1.
           jig   = length m array of group indices vs node index.
                   Node index j is in group jig[j].
           jr    = length ng array of indices representing the groups.
                   The index for group ig is j = jr[ig].
        """
        mper = m//ng
        
        
        jg[:-1] = np.arange(ng)*mper
        jg[ng] = m
        
        ngm1 = ng - 1;
        len1 = ngm1*mper;
        
        jig[:len1] = np.arange(len1)/mper;
        jig[len1:m] = ngm1 

        jr[:ngm1] = ((2*np.arange(ngm1)+1)*mper-1)/2
        jr[ngm1] = (ngm1*mper+m-1)/2
        
    def PrintIntro(self, ):
        
        print "\n\nDemonstration program for CVODE - CVSPGMR linear solver\n"
        print "Food web problem with ns species, ns = %d" % self.ns
        print "Predator-prey interaction and diffusion on a 2-D square\n"

        print "Matrix parameters: a = %.2g   e = %.2g   g = %.2g" % (
         AA, EE, GG)
        print "b parameter = %.2g" % BB
        print "Diffusion coefficients: Dprey = %.2g   Dpred = %.2g" % (
         DPREY, DPRED)
        print "Rate parameter alpha = %.2g\n" % ALPH

        print "Mesh dimensions (mx,my) are %d, %d.  " % (MX, MY)
        print "Total system size is neq = %d \n" % NEQ

        print "Tolerances: reltol = %.2g, abstol = %.2g \n" % (
         RTOL, ATOL)

        print "Preconditioning uses a product of:\n" + \
              "  (1) Gauss-Seidel iterations with " + \
              "itmax = %d iterations, and\n" % ITMAX + \
              "  (2) interaction-only block-diagonal matrix " + \
              "with block-grouping\n" + \
              "  Number of diagonal block groups = ngrp = %d" % NGRP + \
              "  (ngx by ngy, ngx = %d, ngy = %d)\n" % (NGX, NGY) + \
              "\n\n--------------------------------------------------------------" + \
              "--------------\n"


    def CInit(self, c,):
        """
        This routine computes and loads the vector of initial values. 
        """

        #int jx, jy, ns, mxns, ioff, iyoff, i, ici;
        #realtype argx, argy, x, y, dx, dy, x_factor, y_factor, *cdata;
  
        cdata = c.data
        ns = self.ns;
        #mxns = self.mxns;
        dx = self.dx;
        dy = self.dy;

        x_factor = 4.0/np.sqrt(AX)
        y_factor = 4.0/np.sqrt(AY)
        for jy in range(MY):
            y = jy*dy
            argy = np.sqrt(y_factor*y*(AY-y));
            #iyoff = mxns*jy
            for jx in range(MX):
                x = jx*dx
                argx = np.sqrt(x_factor*x*(AX-x))
                #ioff = iyoff + ns*jx
                for i in range(1,ns+1):
                    #ici = ioff + i-1
                    cdata[i-1,jx,jy] = 10.0 + i*argx*argy;



    def PrintHeader(self, jpre, gstype):
        print "\n\nPreconditioner type is           jpre = PREC_%s" % jpre.upper()
        print "\nGram-Schmidt method type is    gstype = %s_GS\n\n" % gstype.upper()
  
    def PrintAllSpecies(self, c, t):
        
        ns = self.ns
        #msns = self.mxns

  
        cdata = c.data
        print "c values at t = %g:\n" % t

        #for (i=1; i <= ns; i++) {
        for i in range(1,ns+1):
            print "Species %d" % i
            for jy in range(MY):
                for jx in range(MX):
                    print "%-10.6g" % cdata[(i-1),jx,-(jy+1)], 
                print ""
            print ""
            
    def PrintOutput(self, t):
        
        nst = self.numSteps
        nfe = self.numRhsEvals
        nni = self.numNonlinSolvIters
        qu = self.lastOrder
        hu = self.lastStep
    

        print "t = %10.2e  nst = %ld  nfe = %ld  nni = %ld" % (t, nst, nfe, nni),
        print "  qu = %d  hu = %11.2e" % (qu, hu)

    def PrintFinalStats(self, ):
        
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

        print "\n\n Final statistics for this run:\n"
        print " CVode real workspace length           = %4ld " % lenrw
        print " CVode integer workspace length        = %4ld " % leniw
        print " CVSPGMR real workspace length         = %4ld " % lenrwLS
        print " CVSPGMR integer workspace length      = %4ld " % leniwLS
        print " Number of steps                       = %4ld " % nst
        print " Number of f-s                         = %4ld " % nfe
        print " Number of f-s (SPGMR)                 = %4ld " % nfeLS
        print " Number of f-s (TOTAL)                 = %4ld " % (nfe + nfeLS)
        print " Number of setups                      = %4ld " % nsetups
        print " Number of nonlinear iterations        = %4ld " % nni
        print " Number of linear iterations           = %4ld " % nli
        print " Number of preconditioner evaluations  = %4ld " % npe
        print " Number of preconditioner solves       = %4ld " % nps
        print " Number of error test failures         = %4ld " % netf
        print " Number of nonlinear conv. failures    = %4ld " % ncfn
        print " Number of linear convergence failures = %4ld " % ncfl
        avdim = nli/nni if nni > 0 else ZERO
        print " Average Krylov subspace dimension     = %.3f " % avdim

        print "\n\n----------------------------------------------------------------------------"
        print     "----------------------------------------------------------------------------"

    
    def RhsFn(self, t, c, cdot,):
        """
        This routine computes the right-hand side of the ODE system and
        returns it in cdot. The interaction rates are computed by calls to WebRates,
        and these are saved in fsave for use in preconditioning.
        """
        cdata = c.data
        cdotdata = cdot.data
  
        #mxns = self.mxns
        ns = self.ns
        fsave = self.fsave
        cox = self.cox
        coy = self.coy
        #mxns = self.mxns
        dx = self.dx
        dy = self.dy
  

        for jy  in range(MY):
            y = jy*dy
            #iyoff = mxns*jy
            idyu = -1 if (jy == MY-1) else 1
            idyl = -1 if (jy == 0) else 1

            for jx in range(MX):
                x = jx*dx
                #ic = iyoff + ns*jx
                # Get interaction rates at one point (x,y). 
                self.WebRates(x, y, t, cdata[:,jx,jy], fsave[:,jx,jy])
                idxu = -1 if (jx == MX-1) else 1
                idxl = -1 if (jx == 0) else 1
                for i in range(ns):
                    #ici = ic + i-1;
                    # Do differencing in y. 
                    dcyli = cdata[i,jx,jy] - cdata[i,jx,jy-idyl]
                    dcyui = cdata[i,jx,jy+idyu] - cdata[i,jx,jy]
                    # Do differencing in x. 
                    dcxli = cdata[i,jx,jy] - cdata[i,jx-idxl,jy]
                    dcxui = cdata[i,jx+idxu,jy] - cdata[i,jx,jy]
                    # Collect terms and load cdot elements. 
                    cdotdata[i,jx,jy] = coy[i]*(dcyui - dcyli) + cox[i]*(dcxui - dcxli) + fsave[i,jx,jy]

        #print i, cdotdata[...], cdata[...], fsave[...]
        #print dcyli, dcyui, dcxli, dcxui
        return 0



    def WebRates(self, x, y, t, c, rate, ):
        """
        This routine computes the interaction rates for the species
        c_1, ... ,c_ns (stored in c[0],...,c[ns-1]), at one spatial point 
        and at time t.
        """
        ns = self.ns
        acoef = self.acoef
        bcoef = self.bcoef
        
  
        for i in range(ns):
            rate[i] = ZERO
              
        for j in range(ns):
            for i in range(ns):
                rate[i] += c[j] * acoef[i,j]  
        fac = ONE + ALPH*x*y
        
        for i in range(ns):
            #print c[i], bcoef[i], fac, rate[i]
            rate[i] = c[i]*(bcoef[i]*fac + rate[i])



    def SpilsPrecSetup(self, t, c, fc, jok, jcurPtr, gamma, vtemp1, vtemp2, vtemp3):
        """
        This routine generates the block-diagonal part of the Jacobian
        corresponding to the interaction rates, multiplies by -gamma, adds
        the identity matrix, and calls denseGETRF to do the LU decomposition of
        each diagonal block. The computation of the diagonal blocks uses
        the preset block and grouping information. One block per group is
        computed. The Jacobian elements are generated by difference
        quotients using calls to the routine fblock.

        This routine can be regarded as a prototype for the general case
        of a block-diagonal preconditioner. The blocks are of size mp, and
        there are ngrp=ngx*ngy blocks computed in the block-grouping scheme.
        """
  
        cdata = c.data
        rewt = self.rewt
        self.errWeights(rewt)
  
        rewtdata = rewt.data

        uround = UNIT_ROUNDOFF

        P = self.P
        pivot = self.pivot
        jxr = self.jxr
        jyr = self.jyr
        mp = self.mp
        srur = self.srur
        #ngrp = self.ngrp
        ngx = self.ngx
        ngy = self.ngy
        #mxmp = self.mxmp
        fsave = self.fsave
        """
        Make mp calls to fblock to approximate each diagonal block of Jacobian.
        Here, fsave contains the base value of the rate vector and 
        r0 is a minimum increment factor for the difference quotient. 
        """
  
        f1 = vtemp1.data
  
        fac = fc.WrmsNorm(rewt)
        r0 = 1000.0*np.abs(gamma)*uround*NEQ*fac
        if (r0 == ZERO): r0 = ONE
  

        for igy in range(ngy):
            jy = jyr[igy]
            #if00 = jy*mxmp
            for igx in range(ngx):
                jx = jxr[igx]
                #if0 = if00 + jx*mp
                #ig = igx + igy*ngx
                # Generate ig-th diagonal block 
                for j in range(mp):
                    # Generate the jth column as a difference quotient 
                    #jj = if0 + j
                    save = cdata[:,igx,igy]
                    r = np.max(srur*np.abs(save),r0/rewtdata[:,igx,igy])
                    cdata[:,igx,igy] += r
                    fac = -gamma/r
                    self.fblock (t, cdata, jx, jy, f1)
                    for i in range(mp):        
                        P[igx,igy,j,i] = (f1[i] - fsave[i,jx,jy])*fac;
    
                    cdata[j,jx,jy] = save;
      
  
        # Add identity matrix and do LU decompositions on blocks. 
  
        for igy in range(ngy):
            for igx in range(ngx):
                denseAddIdentity(P[igx,igy,...]);
                ier = denseGETRF(P[igx,igy,...], pivot[igx,igy,...])
                if (ier != 0): return 1
                  
            jcurPtr = True
              
        return 0
        
    def fblock(self, t, cdata, jx, jy, cdotdata, ):
        """
        This routine computes one block of the interaction terms of the
        system, namely block (jx,jy), for use in preconditioning.
        Here jx and jy count from 0.
        """

  
        #iblok = jx + jy*(wdata->mx);
        y = jy*self.dy
        x = jx*self.dx
        #ic = (wdata->ns)*(iblok);
        self.WebRates(x, y, t, cdata[:,jx,jy], cdotdata[:,jx,jy])
        
        
    def SpilsPrecSolve(self, tn, c, fc, r, z, gamma, delta, lr, vtemp):
        """
        This routine applies two inverse preconditioner matrices
        to the vector r, using the interaction-only block-diagonal Jacobian
        with block-grouping, denoted Jr, and Gauss-Seidel applied to the
        diffusion contribution to the Jacobian, denoted Jd.
        It first calls GSIter for a Gauss-Seidel approximation to
        ((I - gamma*Jd)-inverse)*r, and stores the result in z.
        Then it computes ((I - gamma*Jr)-inverse)*z, using LU factors of the
        blocks in P, and pivot information in pivot, and returns the result in z.
        """
  
  
        r.Scale(ONE, z)
  
        # call GSIter for Gauss-Seidel iterations 
  
        self.GSIter(gamma, z, vtemp)
  
        # Do backsolves for inverse of block-diagonal preconditioner factor 
  
        P = self.P
        pivot = self.pivot
        mx = self.mx
        my = self.my
        #ngx = self.ngx
        mp = self.mp
        jigx = self.jigx
        jigy = self.jigy
  
        iv = 0
        for jy in range(my):
            igy = jigy[jy]
            for jx in range(mx):
                igx = jigx[jx]
                #ig = igx + igy*ngx
                denseGETRS(P[igx,igy,...], mp, pivot[igx,igy,:], z.data[iv])
                iv += mp;
  
        return 0

    def GSIter(self, gamma, z, x):
        """
        This routine performs ITMAX=5 Gauss-Seidel iterations to compute an
        approximation to (P-inverse)*z, where P = I - gamma*Jd, and
        Jd represents the diffusion contributions to the Jacobian.
        The answer is stored in z on return, and x is a temporary vector.
        The dimensions below assume a global constant NS >= ns.
        Some inner loops of length ns are implemented with the small
        vector kernels v_sum_prods, v_prod, v_inc_by_prod.
        """

        xd = x.data
        zd = z.data
        ns = self.ns
        mx = self.mx
        my = self.my
        mxns = self.mxns
        cox = self.cox
        coy = self.coy
        
        beta = np.zeros(NS)
        beta2 = np.zeros(NS)
        cof1 = np.zeros(NS)
        gam = np.zeros(NS)
        gam2 = np.zeros(NS)
  
        # Write matrix as P = D - L - U.
        # Load local arrays beta, beta2, gam, gam2, and cof1. 
  
        for i in range(ns):
  
            temp = ONE/(ONE + 2.0*gamma*(cox[i] + coy[i]))
            beta[i] = gamma*cox[i]*temp
            beta2[i] = 2.0*beta[i]
            gam[i] = gamma*coy[i]*temp
            gam2[i] = 2.0*gam[i]
            cof1[i] = temp
  
  
        # Begin iteration loop.
        # Load vector x with (D-inverse)*z for first iteration. 
        
        for jy in range(my):
            #iyoff = mxns*jy
            for jx in range(mx):
                #ic = iyoff + ns*jx
                v_prod(xd[:,jx,jy], cof1, zd[:,jx,jy]) # x[ic+i] = cof1[i]z[ic+i] 
    
        z.Const(ZERO)
  
        # Looping point for iterations. 
        
        for iter in range(1,ITMAX+1):

            # Calculate (D-inverse)*U*x if not the first iteration. 
    
            if iter > 1:
                for jy in range(my):
                    #iyoff = mxns*jy;
                    for jx in range(mx): # order of loops matters 
                        #ic = iyoff + ns*jx
                        x_loc = 0 if jx == 0 else 2 if jx == mx-1 else 1
                        y_loc = 0 if jy == 0 else 2 if jx == my-1 else 1
                        
                        switch = 3*y_loc+x_loc
                        
                        if switch == 0:
                            v_sum_prods(xd[jx,jy,...], beta2, xd[jx+1,jy,...], gam2, xd[jx,jy+1,...], ns)
                        elif switch == 1:
                            v_sum_prods(xd[jx,jy,...], beta, xd[jx+1,jy,...], gam2, xd[jx,jy,...], ns)
                        elif switch == 2:
                            v_prod(xd[jx,jy,...], gam2, xd[jx,jy+1,...], ns)
                        elif switch == 3:
                            v_sum_prods(xd[jx,jy,...], beta2, xd[jx+1,jy,...], gam, xd[jx,jy+1,...], ns)
                        elif switch == 4:
                            v_sum_prods(xd[jx,jy,...], beta, xd[jx+1,jy,...], gam, xd[jx,jy+1,...], ns);
                        elif switch == 5:
                            v_prod(xd[jx,jy,...], gam, xd[jx,jy+1,...], ns);
                        elif switch == 6:
                            v_prod(xd[jx,jy,...], beta2, xd[jx+1,jy,...], ns);
                        elif switch == 7:
                            v_prod(xd[jx,jy,...], beta, xd[jx+1,jy,...], ns);
                        elif switch == 8:
                            v_zero(xd[jx,jy,...], ns)
                        else:
                            return 1
                            

            # end if (iter > 1) 
    
            # Overwrite x with [(I - (D-inverse)*L)-inverse]*x. 
    
            for jy in range(my):
                i#yoff = mxns*jy;
                for jx in range(mx): # order of loops matters 
                    #ic = iyoff + ns*jx;
                    x_loc = 0 if jx == 0 else 2 if jx == mx-1 else 1
                    y_loc = 0 if jy == 0 else 2 if jx == my-1 else 1
                    switch = 3*y_loc+x_loc
                    
                    if switch == 0:
                        pass
                    elif switch == 1:
                        # 1 <= jx <= mx-2, jy == 0 
                        # x[ic+i] += beta[i]x[ic-ns+i] 
                        v_inc_by_prod(xd[jx,jy,...], beta, xd[jx-1,jy,...], ns);
                    elif switch == 2 : 
                        # jx == mx-1, jy == 0 
                        # x[ic+i] += beta2[i]x[ic-ns+i] 
                        v_inc_by_prod(xd[jx,jy,...], beta2, xd[jx-1,jy,...], ns);
                    elif switch == 3 : 
                        # jx == 0, 1 <= jy <= my-2 
                        # x[ic+i] += gam[i]x[ic-mxns+i] 
                        v_inc_by_prod(xd[jx,jy,...], gam, xd[jx,jy-1,...], ns);
                    elif switch == 4:
                        # 1 <= jx <= mx-2, 1 <= jy <= my-2 
                        # x[ic+i] += beta[i]x[ic-ns+i] + gam[i]x[ic-mxns+i] 
                        v_inc_by_prod(xd[jx,jy,...], beta, xd[jx-1,jy,...], ns);
                        v_inc_by_prod(xd[jx,jy,...], gam, xd[jx,jy-1,...], ns);
                    elif switch == 5 : 
                        # jx == mx-1, 1 <= jy <= my-2 
                        # x[ic+i] += beta2[i]x[ic-ns+i] + gam[i]x[ic-mxns+i] 
                        v_inc_by_prod(xd[jx,jy,...], beta2, xd[jx-1,jy,...], ns);
                        v_inc_by_prod(xd[jx,jy,...], gam, xd[jx,jy-1,...], ns);
                    elif switch == 6 : 
                        # jx == 0, jy == my-1 
                        # x[ic+i] += gam2[i]x[ic-mxns+i] 
                        v_inc_by_prod(xd[jx,jy,...], gam2, xd[jx,jy-1,...], ns);
                    elif switch == 7:
                        # 1 <= jx <= mx-2, jy == my-1 
                        # x[ic+i] += beta[i]x[ic-ns+i] + gam2[i]x[ic-mxns+i] 
                        v_inc_by_prod(xd[jx,jy,...], beta, xd[jx-1,jy,...], ns);
                        v_inc_by_prod(xd[jx,jy,...], gam2, xd[jx,jy-1,...], ns);
                    elif switch == 8 : 
                        # jx == mx-1, jy == my-1 
                        # x[ic+i] += beta2[i]x[ic-ns+i] + gam2[i]x[ic-mxns+i] 
                        v_inc_by_prod(xd[jx,jy,...], beta2, xd[jx-1,jy,...], ns);
                        v_inc_by_prod(xd[jx,jy,...], gam2, xd[jx,jy-1,...], ns);
                        
            # Add increment x to z : z <- z+x 
    
            z.LinearSum(ONE, ONE, x, z);
    
    
    
def v_inc_by_prod(u,v,w,n):
    u[...] += v*w
  
  
def v_sum_prods(u,p,q,v,w,n):
    u[...] = p*q + v*w

def v_prod(u, v, w, n):
    u[...] = w=v*w

def v_zero(u, n):
    u[...] = ZERO
#
#typedef struct {
#  realtype **P[NGRP];
#  long int *pivot[NGRP];
#  int ns, mxns;
#  int mp, mq, mx, my, ngrp, ngx, ngy, mxmp;
#  int jgx[NGX+1], jgy[NGY+1], jigx[MX], jigy[MY];
#  int jxr[NGX], jyr[NGY];
#  realtype acoef[NS][NS], bcoef[NS], diff[NS];
#  realtype cox[NS], coy[NS], dx, dy, srur;
#  realtype fsave[NEQ];
#  N_Vector rewt;
#  void *cvode_mem;
#} *WebData;


# Implementation 

def main():
    abstol=ATOL; reltol=RTOL
    
    

    # Initializations 
    c = NvectorNdarrayFloat64(SHP)
  
    wdata = WebData()

    ns = wdata.ns
    mxns = wdata.mxns

    # Print problem description 
    wdata.PrintIntro()

    # Loop over jpre and gstype (four cases) 
    for jpre in ('left','right',):
        for gstype in ('modified', 'classical'):
            # Initialize c and print heading 
            wdata.CInit(c)
            wdata.PrintHeader(jpre, gstype);

            # Call CVodeInit or CVodeReInit, then CVSpgmr to set up problem 
      
            firstrun = jpre == 'left' and gstype == 'modified'
            if firstrun:
                # Already done above in WebData cinit
                #cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
        
                wdata.initSolver(T0, c)
                wdata.setTolerances(reltol, abstol)
        
                wdata.setupIndirectLinearSolver('spgmr', prectype=jpre, maxl=MAXL, 
                                                gstype=gstype, epslin=DELT)
        

            else:
                
                wdata.reInitSolver(T0, c)
                wdata.spilsPrecType = jpre
                wdata.spilsGSType = gstype

        
      
            # Print initial values 
            if firstrun:
                wdata.PrintAllSpecies(c, T0)
          
            # Loop over output points, call CVode, print sample solution values. 
            tout = T1
          
            for iout in range(1, NOUT):
                flag, t = wdata.Solve(tout, c)
            
                wdata.PrintOutput(t);
                if firstrun and iout % 3 == 0:
                    wdata.PrintAllSpecies(c, t);
            
                if (tout > 0.9):
                    tout += DTOUT
                else:
                    tout *= TOUT_MULT; 
          
          
            # Print final statistics, and loop for next case 
            wdata.PrintFinalStats()
      


if __name__ == '__main__':
    main()

