cimport libsundials as sun
cimport libcvode as cvode

from sundials cimport N_Vector, pyDlsMat

from libc.stdlib cimport abort, malloc, free
from cpython cimport Py_INCREF, Py_DECREF

import cython

import numpy as np
cimport numpy as np
np.import_array() # initialize C API to call PyArray_SimpleNewFromData

include 'denseGET.pxi'

class CVodeError(Exception):
    pass

include 'cvode_properties.pxi'

PRE_CONDITION_CONSTS = {
    'none': sun.PREC_NONE, 'left': sun.PREC_LEFT, 
    'right': sun.PREC_RIGHT, 'both': sun.PREC_BOTH
}


cdef class Cvode(BaseCvode):
    
    def __init__(self, *args, **kwds):        
        self._nrtfn = kwds.pop('nrtfn', 0)
        
    @property
    def nrtfn(self, ):
        return self._nrtfn
        
        
    @property
    def multistep(self, ):
        if self._ms == cvode.CV_BDF:
            return 'bdf'.decode('UTF-8', 'strict')
        if self._ms == cvode.CV_ADAMS:
            return 'adams'.decode('UTF-8', 'strict')
        return None
        
    @property
    def iteration(self, ):
        if self._it == cvode.CV_FUNCTIONAL:
            return 'functional'.decode('UTF-8', 'strict')
        if self._it == cvode.CV_NEWTON:
            return 'newton'.decode('UTF-8', 'strict')
        return None        
        
    def initSolver(self, sun.realtype t0, N_Vector y0,  ):

        ret = cvode.CVodeInit(self._cv, _CvRhsFn, t0, y0._v)
        if ret != 0:
            raise CVodeError()
            
        if self._nrtfn > 0:
            ret = cvode.CVodeRootInit(self._cv, self._nrtfn, _CvRootFn, )
            if ret != 0:
                raise CVodeError()
   
    def setTolerances(self, sun.realtype reltol, abstol):
        
        if isinstance(abstol, N_Vector):
            ret = cvode.CVodeSVtolerances(self._cv, reltol, (<N_Vector>abstol)._v)
        else:
            ret = cvode.CVodeSStolerances(self._cv, reltol, abstol)
        if ret != 0:
            raise CVodeError()
            
    
    def errWeights(self, N_Vector nv):  

        ret = cvode.CVodeGetErrWeights(self._cv, nv._v)
        if ret == cvode.CV_SUCCESS:
            return
        if ret == cvode.CV_MEM_NULL:
            raise ValueError("CVODE memory pointer is NULL")
        if ret == cvode.CV_NO_SLDET:
            raise ValueError("Stability limit was not turned on")
        raise ValueError("Unknown error")
        

            
    def estLocalErrors(self, N_Vector nv):    

        ret = cvode.CVodeGetEstLocalErrors(self._cv, nv._v)
        if ret == cvode.CV_SUCCESS:
            return
        if ret == cvode.CV_MEM_NULL:
            raise ValueError("CVODE memory pointer is NULL")
        if ret == cvode.CV_NO_SLDET:
            raise ValueError("Stability limit was not turned on")
        raise ValueError("Unknown error")

    def rootInfo(self, ):
        
        
        if self._nrtfn < 1:
            raise CVodeError('Root finding not enabled (nrtfn < 1)')
            
        cdef int[::1] roots = np.empty( self._nrtfn, dtype=np.int32 )
        
        ret = cvode.CVodeGetRootInfo(self._cv, &roots[0])
        if ret == cvode.CV_SUCCESS:
            return roots
        if ret == cvode.CV_MEM_NULL:
            raise ValueError("CVODE memory pointer is NULL")
        if ret == cvode.CV_NO_SLDET:
            raise ValueError("Stability limit was not turned on")
        raise ValueError("Unknown error")

    ###################            
    # Linear solver setup routines
    #
    ###################
        
        
    def setupDenseLinearSolver(self, long int N, user_jac=False ):
        """
        Initialise Dense Linear Solver
        
        if user_jac is True then the subclass must reimplement Cvode.DlsDenseJacFn
        to perform the appropriate jacobian calculations.
        """
        ret = cvode.CVDense(self._cv, N)
        if user_jac:
            ret = cvode.CVDlsSetDenseJacFn(self._cv, _CvDlsDenseJacFn)
            
    def setupBandLinearSolver(self, long int N, long int mupper, long int mlower, user_jac=False ):
        """
        Initialise Dense Linear Solver
        
        if user_jac is True then the subclass must reimplement Cvode.DlsDenseJacFn
        to perform the appropriate jacobian calculations.
        """
        ret = cvode.CVBand(self._cv, N, mupper, mlower)
        if user_jac:
            ret = cvode.CVDlsSetBandJacFn(self._cv, _CvDlsBandJacFn)
            
    def setupIndirectLinearSolver(self, solver='spgmr', pretype='none', int maxl=0, user_pre=False, user_jac=False ):
        """
        Initialise Indirect Linear Solver
        
        
        """
        
        cdef int pre = PRE_CONDITION_CONSTS[pretype]
        
        
        if solver == 'spgmr':
            ret = cvode.CVSpgmr(self._cv, pre, maxl)
        elif solver == 'spbcg':
            ret = cvode.CVSpbcg(self._cv, pre, maxl)
        elif solver == 'sptfqmr':
            ret = cvode.CVSptfqmr(self._cv, pre, maxl)
        else:
            raise ValueError("Solver name not recognised")
            
        if ret == cvode.CVSPILS_MEM_NULL:
            raise ValueError("Null cvode memory pointer given")
        elif ret == cvode.CVSPILS_MEM_FAIL:
            raise ValueError("Allocating memory for linear solver failed")
        elif ret == cvode.CVSPILS_ILL_INPUT:
            raise ValueError("Illegal input")
        elif ret != cvode.CVSPILS_SUCCESS:
            raise ValueError("Unknown Error ({})".format(ret))
            
            
        if user_pre:
            ret = cvode.CVSpilsSetPreconditioner(self._cv, _CvSpilsPrecSetupFn,
                                                   _CvSpilsPrecSolveFn)
            self._handleSpilsSetReturn(ret, 'SpilsSetPreconditioner')
        if user_jac:
            ret = cvode.CVSpilsSetJacTimesVecFn(self._cv, _CvSpilsJacTimesVecFn)            
            self._handleSpilsSetReturn(ret, 'SpilsSetJacTimesVecFn')
        

        
    def Solve(self, sun.realtype tout, N_Vector yout):
        cdef sun.realtype tret     
                
        #tret = <sun.realtype *>malloc(cython.sizeof(sun.realtype))
        cdef int flag = cvode.CVode(self._cv, tout, yout._v, &tret, cvode.CV_NORMAL)   
        return cvode.CVodeGetReturnFlagName(flag), tret
        
    def RhsFn(self, t, y, ydot):
        raise NotImplementedError()
        
    def RootFn(self, t, y, np.float64_t[::1] gout):
        raise NotImplementedError()
        
    def DlsDenseJacFn(self, N, t, y, fy, J, tmp1, tmp2, tmp3):
        raise NotImplementedError()
        
    def DlsBandJacFn(self, N, mupper, mlower, t, y, fy, J, tmp1, tmp2, tmp3):        
        raise NotImplementedError()
        
    def SpilsJacTimesVec(self, uu, v, Jv, new_uu):
        raise NotImplementedError()
        
    def SpilsPrecSetup(self, uu, uscale, fval, fscale,
                                       vtemp1, vtemp2):    
        raise NotImplementedError()
        
        
cdef int _CvRhsFn(sun.realtype t, sun.N_Vector y,
    		       sun.N_Vector ydot, void *user_data):
    cdef object obj
    
    obj = <object>user_data
    
    pyy = <object>y.content
    pyydot = <object>ydot.content
    return obj.RhsFn(t, pyy, pyydot)
    
    
cdef int _CvRootFn(sun.realtype t, sun.N_Vector y, sun.realtype *gout,
                   void *user_data):
                       
    cdef object obj
    
    obj = <object>user_data
    
    pyy = <object>y.content
    cdef Py_ssize_t N = obj.nrtfn
    cdef np.float64_t[::1] pygout = <np.float64_t [:N]> gout
    
    try:
        return obj.RootFn(t, pyy, pygout)
    except Exception:
        return -1
    
cdef int _CvDlsDenseJacFn(long int N, sun.realtype t,
			       sun.N_Vector y, sun.N_Vector fy, 
			       sun.DlsMat Jac, void *user_data,
			       sun.N_Vector tmp1, sun.N_Vector tmp2, sun.N_Vector tmp3):
              
    cdef object obj = <object>user_data              
              
    pyy = <object>y.content
    pyfy = <object>fy.content
    pytmp1 = <object>tmp1.content
    pytmp2 = <object>tmp2.content
    pytmp3 = <object>tmp3.content
    
    pyJ = pyDlsMat()
    pyJ._m = Jac
              
    #try:
    return obj.DlsDenseJacFn(N, t, pyy, pyfy, pyJ, pytmp1, pytmp2, pytmp3)
        #raise NotImplementedError("Waiting on Cython wrapper of DlsMat")
    #except Exception:
    #    return -1
    
cdef int _CvDlsBandJacFn(long int N, long int mupper, long int mlower,
			      sun.realtype t, sun.N_Vector y, sun.N_Vector fy, 
			      sun.DlsMat Jac, void *user_data,
			      sun.N_Vector tmp1, sun.N_Vector tmp2, sun.N_Vector tmp3):
                  
    try:
        raise NotImplementedError("Waiting on Cython wrapper of DlsMat")       
    except Exception:
        return -1
    
    
cdef int _CvSpilsJacTimesVecFn(sun.N_Vector v, sun.N_Vector Jv, sun.realtype t,
                                        sun.N_Vector y, sun.N_Vector fy,
                                        void *user_data, sun.N_Vector tmp):
    cdef object obj
    
    obj = <object>user_data
    
    
    pyv = <object>v.content
    pyJv = <object>Jv.content
    pyy = <object>y.content
    pyfy = <object>fy.content
    pytmp = <object>tmp.content
    
    
    try:
        return obj.JacTimesVec(pyv, pyJv, t, pyy, pyfy, pytmp )
    except Exception:
        return -1
    
    

    

cdef int _CvSpilsPrecSetupFn(sun.realtype t, sun.N_Vector y, sun.N_Vector fy,
                                      sun.booleantype jok, sun.booleantype *jcurPtr,
                                      sun.realtype gamma, void *user_data,
                                      sun.N_Vector tmp1, sun.N_Vector tmp2,
                                      sun.N_Vector tmp3):
               
    cdef object obj
    
    obj = <object>user_data
    
    pyy = <object>y.content
    pyfy = <object>fy.content
    
    
    
    pyvtemp1 = <object>tmp1.content
    pyvtemp2 = <object>tmp2.content
    pyvtemp2 = <object>tmp3.content
    
    try:
        return obj.SpilsPrecSetup(t, pyy, pyfy)
    except Exception:
        return -1
    
cdef int _CvSpilsPrecSolveFn(sun.realtype t, sun.N_Vector y, sun.N_Vector fy,
                                      sun.N_Vector r, sun.N_Vector z,
                                      sun.realtype gamma, sun.realtype delta,
                                      int lr, void *user_data, sun.N_Vector tmp):
    
    cdef object obj    
    obj = <object>user_data
    
    pyy = <object>y.content
    pyfy = <object>fy.content
    pyr = <object>r.content
    pyz = <object>z.content
    pytmp = <object>tmp.content
    
    try:
        return obj.SpilsPrecSolve(t, pyy, pyfy, pyr, pyz, pytmp)     
    except Exception:
        return -1