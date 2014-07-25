cimport libsundials as sun
cimport libcvodes as cvode

from sundials cimport N_Vector, pyDlsMat
from sundials import PRE_CONDITION_CONSTS, GS_TYPES

from libc.stdlib cimport abort, malloc, free
from cpython cimport Py_INCREF, Py_DECREF

import cython

import numpy as np
cimport numpy as np
np.import_array() # initialize C API to call PyArray_SimpleNewFromData



class CVodeError(Exception):
    pass

include 'denseGET.pxi'
include 'cvodes_properties.pxi'


cdef class Cvodes(BaseCvodes):
    
    def initSens(self, int num_sensitivities, method, yS0, rhs_type='all'):
        cdef int ism
        
        if method == 'simultaneous':
            ism = cvode.CV_SIMULTANEOUS
        elif method == 'staggered':
            ism = cvode.CV_STAGGERED
        elif method == 'staggered1':
            ism = cvode.CV_STAGGERED1
        else:
            raise CVodeError()

        # Convert list of N_Vectors into pointers           
        cdef int n = len(yS0)
        cdef int i
        cdef sun.N_Vector *yS0_ptrs = <sun.N_Vector *>malloc(n*cython.sizeof(<sun.N_Vector*>NULL))
        
        for i in range(n):
            yS0_ptrs[i] = (<N_Vector>yS0[i])._v
        
        if rhs_type == 'all':
            ret = cvode.CVodeSensInit(self._cv, num_sensitivities, ism,
                                  _CvSensRhsFn, yS0_ptrs)
        else:
            ret = cvode.CVodeSensInit1(self._cv, num_sensitivities, ism,
                                  _CvSensRhsFn1, yS0_ptrs)
        if ret != 0:
            raise CVodeError()

    def setSensTolerances(self, reltol=None, abstol=None):
        
        if reltol is None and abstol is None:
            ret = cvode.CVodeSensEEtolerances(self._cv)
        else:
            raise NotImplementedError()
        if ret != 0:
            raise CVodeError()
            
    def sensErrCon(self, sun.booleantype errconS):
        
        ret = cvode.CVodeSetSensErrCon(self._cv, errconS)
        if ret != 0:
            raise CVodeError()

    def sensParams(self, sun.realtype[::1] p=None, sun.realtype[::1] pbar=None, int[::1] plist=None):
        
        cdef sun.realtype* p_ = NULL
        cdef sun.realtype* pbar_ = NULL
        cdef int* plist_ = NULL
        
        if not p is None: p_ = &p[0]
        if not pbar is None: pbar_ = &pbar[0]
        if not plist is None: plist_ = &plist[0]
            
        
        
        ret = cvode.CVodeSetSensParams(self._cv, p_, pbar_, plist_)
        if ret != 0:
            raise CVodeError()
            
            
    def SolveSens(self, ysout):
        cdef sun.realtype tret     
        
        # Convert list of N_Vectors into pointers           
        cdef int n = len(ysout)
        cdef int i
        cdef sun.N_Vector *yS_ptrs = <sun.N_Vector *>malloc(n*cython.sizeof(<sun.N_Vector*>NULL))
        
        for i in range(n):
            yS_ptrs[i] = (<N_Vector>ysout[i])._v        
        
        cdef int flag = cvode.CVodeGetSens(self._cv, &tret, yS_ptrs)
        return cvode.CVodeGetReturnFlagName(flag), tret
        
        
    def SensRhsFn(self, Ns, t, y, ydot, yS, ySdot, tmp1, tmp2):
        raise NotImplementedError()
        
    def SensRhsFn1(self, Ns, t, y, ydot, iS, yS, ySdot, tmp1, tmp2):
        raise NotImplementedError()        

cdef int _CvSensRhsFn(int Ns, sun.realtype t,
			   sun.N_Vector y, sun.N_Vector ydot,
			   sun.N_Vector *yS, sun.N_Vector *ySdot,
			   void *user_data,
			   sun.N_Vector tmp1, sun.N_Vector tmp2):
          
    cdef object obj
    cdef int i
    
    obj = <object>user_data
    
    pyy = <object>y.content
    pyydot = <object>ydot.content
    pytmp1 = <object>tmp1.content
    pytmp2 = <object>tmp2.content
    
    pyys = []
    pyysdot = []
    
    for i in range(Ns):
        pyys.append( <object>yS[i].content )
        pyysdot.append( <object>ySdot[i].content )

    return obj.SensRhsFn(Ns, t, pyy, pyydot, pyys, pyysdot, pytmp1, pytmp2)
    

           
cdef int _CvSensRhsFn1(int Ns, sun.realtype t,
			   sun.N_Vector y, sun.N_Vector ydot,
			   int iS, sun.N_Vector yS, sun.N_Vector ySdot,
			   void *user_data,
			   sun.N_Vector tmp1, sun.N_Vector tmp2):
          
    cdef object obj
    cdef int i
    
    obj = <object>user_data
    
    pyy = <object>y.content
    pyydot = <object>ydot.content
    pytmp1 = <object>tmp1.content
    pytmp2 = <object>tmp2.content
    
    pyys = <object>yS.content
    pyysdot = <object>ySdot.content

    return obj.SensRhsFn1(Ns, t, pyy, pyydot, iS, pyys, pyysdot, pytmp1, pytmp2)           