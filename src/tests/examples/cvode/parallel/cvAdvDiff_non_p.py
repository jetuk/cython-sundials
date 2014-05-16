# -*- coding: utf-8 -*-
"""
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2007/10/25 20:03:28 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, George Byrne,
 *                and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the program for
 * its solution by CVODE. The problem is the semi-discrete
 * form of the advection-diffusion equation in 1-D:
 *   du/dt = d^2 u / dx^2 + .5 du/dx
 * on the interval 0 <= x <= 2, and the time interval 0 <= t <= 5.
 * Homogeneous Dirichlet boundary conditions are posed, and the
 * initial condition is the following:
 *   u(x,t=0) = x(2-x)exp(2x) .
 * The PDE is discretized on a uniform grid of size MX+2 with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX.
 * This program solves the problem with the option for nonstiff
 * systems: ADAMS method and functional iteration.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .5, 1.0, ..., 5.
 * Run statistics (optional outputs) are printed at the end.
 *
 * This version uses MPI for user routines.
 * Execute with Number of Processors = N,  with 1 <= N <= MX.
 * -----------------------------------------------------------------
"""

from pySundials.cvode import Cvode
from pySundials.sundials import NvectorMPINdarrayFloat64
import numpy as np
from mpi4py import MPI


# Problem Constants 

ZERO  =0.0

XMAX  =2.0    # domain boundary           
MX    =10             # mesh dimension            
NEQ   =MX             # number of equations       
ATOL  =1.0e-5 # scalar absolute tolerance 
T0    =ZERO           # initial time              
T1    =0.5    # first output time         
DTOUT =0.5    # output time increment     
NOUT  =10             # number of output times    

class AdvDiff(Cvode):
    """
    UserData 
        contains grid constants, parallel machine parameters, work array.
    
    """
    
    def __init__(self, comm, npes, my_pe, dx, hdcoef, hacoef, **kwds):
        Cvode.__init__(self, **kwds)
        
        self.dx = dx
        self.hdcoef = hdcoef
        self.hacoef = hacoef
        
        self.npes = npes
        self.my_pe = my_pe
        self.comm = comm
        self.z = np.empty( 100 )
        
    def SetIC(self, u, my_base):
        """
        Set initial conditions in u vector       
        """


        # Set pointer to data array and get local length of u. 
        udata = u.data
        my_length = udata.shape[0]

        # Load initial profile into u vector 
        for i in range(my_length):      
            iglobal = my_base + i;
            x = iglobal*self.dx;
            udata[i] = x*(XMAX - x)*np.exp(2.0*x);

    def RhsFn(self, t, u, udot):
        """
        f routine. Compute f(t,u). 
        """
        
        udata = u.data
        dudata = udot.data


        # Extract needed problem constants from data 
  
        hordc = self.hdcoef;
        horac = self.hacoef;

        # Extract parameters for parallel computation. 
        comm = self.comm;
        npes = self.npes;           # Number of processes.  
        my_pe = self.my_pe;         # Current process number. 
        my_length = udata.shape[0] # Number of local elements of u.  
        z = self.z;

        # Compute related parameters. 
        my_pe_m1 = my_pe - 1
        my_pe_p1 = my_pe + 1
        last_pe = npes - 1
        my_last = my_length - 1

        # Store local segment of u in the working array z. 
        for i in range(my_length):
            z[i+1] = udata[i];

        # Pass needed data to processes before and after current process. 
        if my_pe != 0:
            #MPI_Send(&z[1], 1, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm);
            comm.send(z[1], dest=my_pe_m1, tag=0)
        if my_pe != last_pe:
            #MPI_Send(&z[my_length], 1, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm);   
            comm.send(z[my_length], dest=my_pe_p1, tag=0)
            

        # Receive needed data from processes before and after current process. 
        if my_pe != 0:
            #MPI_Recv(&z[0], 1, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm, &status);
            z[0] = comm.recv(source=my_pe_m1, tag=0)
        else:
            z[0] = ZERO;
   
        if my_pe != last_pe:
            #MPI_Recv(&z[my_length+1], 1, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm,
            #  &status);   
            z[my_length+1] = comm.recv(source=my_pe_p1, tag=0)
        else:
            z[my_length+1] = ZERO;

        # Loop over all grid points in current process. 
        for i in range(my_length):
            # Extract u at x_i and two neighboring points 
            ui = z[i+1]
            ult = z[i]
            urt = z[i+2]

            # Set diffusion and advection terms and load into udot 
            hdiff = hordc*(ult - 2.0*ui + urt)
            hadv = horac*(urt - ult)
            dudata[i] = hdiff + hadv
        return 0

    def PrintFinalStats(self,):
        """
        Print some final statistics located in the iopt array 
        """
        print "\nFinal Statistics: "
        print "nst = %-6ld  nfe  = %-6ld " % (self.numSteps, self.numRhsEvals)
        print "nni = %-6ld  ncfn = %-6ld  netf = %ld\n" % (self.numNonlinSolvIters,
              self.numNonlinSolvConvFails, self.numErrTestFails)


def main():
    """
    **************************** Main Program *****************************
    """
    # Get processor number, total number of pe's, and my_pe. 
    #MPI_Init(&argc, &argv);
    comm = MPI.COMM_WORLD
    npes = comm.Get_size()
    my_pe = comm.Get_rank()

    # Set local vector length. 
    nperpe = NEQ/npes
    nrem = NEQ - npes*nperpe
    local_N = nperpe+1 if my_pe < nrem else nperpe
    my_base = my_pe*local_N if my_pe < nrem else my_pe*nperpe + nrem
    
    u = NvectorMPINdarrayFloat64(comm, local_N, NEQ)
    
  
    reltol = ZERO;  # Set the tolerances 
    abstol = ATOL;

    dx = XMAX/(float(MX+1))  # Set grid coefficients in data 
    hdcoef = 1.0/(dx*dx)
    hacoef = 0.5/(2.0*dx)

    # Call CVodeCreate to create the solver memory and specify the 
    # Adams-Moulton LMM and the use of a functional iteration   
    data = AdvDiff(comm, npes, my_pe, dx, hdcoef, hacoef, multistep='adams',
                   iteration='functional')

    data.SetIC(u, my_base);  # Initialize u vector 

    # Call CVodeInit to initialize the integrator memory and specify the
    # user's right hand side function in u'=f(t,u), the inital time T0, and
    # the initial dependent variable vector u. 
    data.initSolver( T0, u )

    # Call CVodeSStolerances to specify the scalar relative tolerance
    # and scalar absolute tolerances 
    data.setTolerances(reltol, abstol)
  
    if my_pe == 0: PrintIntro(npes)

    umax = u.MaxNorm()

    if my_pe == 0:
        t = T0
        PrintData(t, umax, 0)

    # In loop over output points, call CVode, print results, test for error 
    tout = T1
    for iout in range(NOUT):
        tout += DTOUT
        flag, t = data.Solve(tout, u)

        if flag == 'CV_SUCCESS':
            umax = u.MaxNorm()        
            nst = data.numSteps

            if my_pe == 0: PrintData(t, umax, nst);

    if my_pe == 0: data.PrintFinalStats()  # Print some final statistics 



#*********************** Private Helper Functions **********************


def PrintIntro(npes):
    """
    Print problem introduction     
    """
    print "\n 1-D advection-diffusion equation, mesh size =%3d" % MX
    print "\n Number of PEs = %3d \n" % npes

def PrintData(t, umax, nst):
    """
    Print data 
    """
    print "At t = %4.2f  max.norm(u) =%14.6e  nst =%4ld" % (t, umax, nst)



if __name__ == '__main__':
    main()
