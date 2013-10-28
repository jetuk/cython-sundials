cimport libsundials as sun
cimport libkinsol as kinsol

from sundials cimport N_Vector

from libc.stdlib cimport abort, malloc, free
from cpython cimport Py_INCREF, Py_DECREF

import cython

import sys

class KinsolError(Exception):
    pass


cdef class Kinsol:
    #cdef void *_kn
    #cdef int _ms
    #cdef int _it
    #cdef N_Vector tmpl
    
    def __cinit__(self, ): 
                     
        
        self._kn = kinsol.KINCreate()
        if not self._kn:
            raise MemoryError
        
        ret = kinsol.KINSetUserData(self._kn, <void *>self)  
        if ret != 0:
            raise KinsolError()
            
        print 'Initialised Kinsol'
        
        Py_INCREF(self)
        
        
    def Setup(self, N_Vector tmpl, 
              user_Jv=False,
              normtol=0.0, 
              scsteptol=0.0,
              mxnewtstep=None, 
              spils_maxl=0, 
              mxiter=200,
              dense=True,
              linear_solv='spgmr',
              user_Presolve=False,
              max_restarts=0,
              eta_choice='choice1'):
        print 'Setting up Kinsol...'
        cdef long N
        N = tmpl.Size()
        self.tmpl = tmpl
        print N
        sys.stdout.write('\t ..Init..')
        ret = kinsol.KINInit(self._kn, _KnRhsFn, tmpl._v)

        print ret
                
        if ret != 0:
            print 'Failed!'
            raise KinsolError()
        print 'Success'
        ret = kinsol.KINSetNumMaxIters(self._kn, mxiter)
        ret = kinsol.KINSetFuncNormTol(self._kn, normtol)        
        ret = kinsol.KINSetScaledStepTol(self._kn, scsteptol)
        print ret, N
        
        if not mxnewtstep is None:
            ret = kinsol.KINSetMaxNewtonStep(self._kn, mxnewtstep)
            print ret
        if dense:
            ret = kinsol.KINDense(self._kn, N)
        else:
            if linear_solv == 'spgmr':
                ret = kinsol.KINSpgmr(self._kn, spils_maxl)
                kinsol.KINSpilsSetMaxRestarts(self._kn, max_restarts)
            elif linear_solv == 'spbcg':
                ret = kinsol.KINSpbcg(self._kn, spils_maxl)
            elif linear_solv == 'sptfqmr':
                ret = kinsol.KINSptfqmr(self._kn, spils_maxl)
            else:
                raise ValueError
        
        print ret
        
        if eta_choice == 'choice1':
            ret = kinsol.KINSetEtaForm(self._kn, kinsol.KIN_ETACHOICE1)
        elif eta_choice == 'choice2':
            ret = kinsol.KINSetEtaForm(self._kn, kinsol.KIN_ETACHOICE2)
        elif eta_choice == 'constant':
            ret = kinsol.KINSetEtaForm(self._kn, kinsol.KIN_ETACONSTANT)
        else:
            raise ValueError
            
        if user_Presolve:
            ret = kinsol.KINSpilsSetPreconditioner(self._kn, _KnSpilsPrecSetupFn,
                                                   _KnSpilsPrecSolveFn)
        
        if user_Jv:
            ret = kinsol.KINSpilsSetJacTimesVecFn(self._kn, _KnJacTimesVecFn)
            
    def SetMaxSetupCalls(self, int msbset=0):
        ret = kinsol.KINSetMaxSetupCalls(self._kn, msbset)
        if ret == kinsol.KIN_SUCCESS:
            return
        if ret == kinsol.KIN_MEM_NULL:
            raise ValueError('Setup first must be called before SetMaxSetupCalls')
        if ret == kinsol.KIN_ILL_INPUT:
            raise ValueError('Illegal value')
        raise ValueError('Unknown error ({}))'.format(ret))
        
    def SetConstraints(self, N_Vector constraints):
        print 'Setting constraints...'
        
        ret = kinsol.KINSetConstraints(self._kn, constraints._v)
        if ret == kinsol.KIN_SUCCESS:
            return        
        if ret == kinsol.KIN_MEM_NULL:
            raise ValueError('Setup first must be called before SetConstraints')
        if ret == kinsol.KIN_ILL_INPUT:
            raise ValueError('Constraints contained illegal values')
        raise ValueError('Unknown error ({}))'.format(ret))
        
    def Solve(self, N_Vector uu, N_Vector u_scale, N_Vector f_scale,
              strategy = 'none', print_level=0):
                  
        cdef int strat
        if strategy == 'none':
            strat = kinsol.KIN_NONE
        elif strategy == 'linesearch':
            strat = kinsol.KIN_LINESEARCH
        else:
            raise ValueError                  
            
        flag = kinsol.KINSetPrintLevel(self._kn, print_level)
        
        flag = kinsol.KINSol(self._kn, uu._v, strat, u_scale._v, f_scale._v)
        
        print flag
        
    def RhsFn(self, uu, fval):
        raise NotImplementedError()        
        
    def JacTimesVec(self, uu, v, Jv, new_uu):
        raise NotImplementedError()
        
    def SpilsPrecSetup(self, uu, uscale, fval, fscale,
                                       vtemp1, vtemp2):
        """
    /*
     * -----------------------------------------------------------------
     * Type : KINSpilsPrecSetupFn
     * -----------------------------------------------------------------
     * The user-supplied preconditioner setup subroutine should
     * compute the right-preconditioner matrix P (stored in memory
     * block referenced by P_data pointer) used to form the
     * scaled preconditioned linear system:
     *
     *  (Df*J(uu)*(P^-1)*(Du^-1)) * (Du*P*x) = Df*(-F(uu))
     *
     * where Du and Df denote the diagonal scaling matrices whose
     * diagonal elements are stored in the vectors uscale and
     * fscale, repsectively.
     *
     * The preconditioner setup routine (referenced by iterative linear
     * solver modules via pset (type KINSpilsPrecSetupFn)) will not be
     * called prior to every call made to the psolve function, but will
     * instead be called only as often as necessary to achieve convergence
     * of the Newton iteration.
     *
     * Note: If the psolve routine requires no preparation, then a
     * preconditioner setup function need not be given.
     *
     *  uu  current iterate (unscaled) [input]
     *
     *  uscale  vector (type N_Vector) containing diagonal elements
     *          of scaling matrix for vector uu [input]
     *
     *  fval  vector (type N_Vector) containing result of nonliear
     *        system function evaluated at current iterate:
     *        fval = F(uu) [input]
     *
     *  fscale  vector (type N_Vector) containing diagonal elements
     *          of scaling matrix for fval [input]
     *
     *  user_data  pointer to user-allocated data memory block
     *
     *  vtemp1/vtemp2  available scratch vectors (temporary storage)
     *
     * If successful, the function should return 0 (zero). If an error
     * occurs, then the routine should return a non-zero integer value.
     * -----------------------------------------------------------------
     */               
        """                            
        raise NotImplementedError()
        
    def SpilsPrecSolve(uu, uscale, fval, fscale, vv, N_Vector vtemp):
        """
    /*
     * -----------------------------------------------------------------
     * Type : KINSpilsPrecSolveFn
     * -----------------------------------------------------------------
     * The user-supplied preconditioner solve subroutine (referenced
     * by iterative linear solver modules via psolve (type
     * KINSpilsPrecSolveFn)) should solve a (scaled) preconditioned
     * linear system of the generic form P*z = r, where P denotes the
     * right-preconditioner matrix computed by the pset routine.
     *
     *  uu  current iterate (unscaled) [input]
     *
     *  uscale  vector (type N_Vector) containing diagonal elements
     *          of scaling matrix for vector uu [input]
     *
     *  fval  vector (type N_Vector) containing result of nonliear
     *        system function evaluated at current iterate:
     *        fval = F(uu) [input]
     *
     *  fscale  vector (type N_Vector) containing diagonal elements
     *          of scaling matrix for fval [input]
     *
     *  vv  vector initially set to the right-hand side vector r, but
     *      which upon return contains a solution of the linear system
     *      P*z = r [input/output]
     *
     *  user_data  pointer to user-allocated data memory block
     *
     *  vtemp  available scratch vector (volatile storage)
     *
     * If successful, the function should return 0 (zero). If a
     * recoverable error occurs, then the subroutine should return
     * a positive integer value (in this case, KINSOL attempts to
     * correct by calling the preconditioner setup function if the 
     * preconditioner information is out of date). If an unrecoverable 
     * error occurs, then the preconditioner solve function should return 
     * a negative integer value.
     * -----------------------------------------------------------------
     */        
        """
        raise NotImplementedError()
        
        
    def PrintFinalStats(self, ):
        cdef long int nni
        cdef long int nfe
        cdef long int nje
        cdef long int nfeD
        
      
        flag = kinsol.KINGetNumNonlinSolvIters(self._kn, &nni);
        flag = kinsol.KINGetNumFuncEvals(self._kn, &nfe)            
        flag = kinsol.KINDlsGetNumJacEvals(self._kn, &nje)
        flag = kinsol.KINDlsGetNumFuncEvals(self._kn, &nfeD)

    
        print "Final Statistics:"
        print "  nni = {:5d}    nfe  = {:5d} ".format(nni, nfe)
        print "  nje = {:5d}    nfeD = {:5d} ".format(nje, nfeD)



cdef int _KnRhsFn(sun.N_Vector uu,
    		       sun.N_Vector fval, void *user_data):
    cdef object obj
    
    obj = <object>user_data
    
    pyuu = <object>uu.content
    pyfval = <object>fval.content
    return obj.RhsFn(pyuu, pyfval)
    
    
cdef int _KnJacTimesVecFn(sun.N_Vector v, sun.N_Vector Jv,
                                         sun.N_Vector uu, sun.booleantype *new_uu, 
                                         void *user_data):
    cdef object obj
    
    obj = <object>user_data
    
    pyuu = <object>uu.content
    pyv = <object>v.content
    pyJv = <object>Jv.content
    pynew_uu = <int>new_uu
    return obj.JacTimesVec(pyuu, pyv, pyJv, pynew_uu)
    

cdef int _KnSpilsPrecSetupFn(sun.N_Vector uu, sun.N_Vector uscale,
                                       sun.N_Vector fval, sun.N_Vector fscale,
                                       void *user_data, sun.N_Vector vtemp1,
    				   sun.N_Vector vtemp2):
               
    cdef object obj
    
    obj = <object>user_data
    
    pyuu = <object>uu.content
    pyuscale = <object>uscale.content
    pyfval = <object>fval.content
    pyfscale = <object>fscale.content
    pyvtemp1 = <object>vtemp1.content
    pyvtemp2 = <object>vtemp2.content
    
    return obj.SpilsPrecSetup(pyuu, pyuscale, pyfval, pyfscale, pyvtemp1, 
                              pyvtemp2)
    
cdef int _KnSpilsPrecSolveFn(sun.N_Vector uu, sun.N_Vector uscale, 
                                       sun.N_Vector fval, sun.N_Vector fscale, 
                                       sun.N_Vector vv, void *user_data,
                                       sun.N_Vector vtemp):
    
    cdef object obj    
    obj = <object>user_data
    
    pyuu = <object>uu.content
    pyuscale = <object>uscale.content
    pyfval = <object>fval.content
    pyfscale = <object>fscale.content
    pyvtemp = <object>vtemp.content
    
    
    return obj.SpilsPrecSolve(pyuu, pyuscale, pyfval, pyfscale, pyvtemp)    
        