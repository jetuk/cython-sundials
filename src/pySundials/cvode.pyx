cimport libsundials as sun
cimport libcvode as cvode

from sundials cimport N_Vector

from libc.stdlib cimport abort, malloc, free
from cpython cimport Py_INCREF, Py_DECREF

import cython


cdef class Cvode:
    #cdef void *_cv
    #cdef int _ms
    #cdef int _it
    #cdef N_Vector y0
    
    def __cinit__(self, multistep='bdf', iteration='functional'): 
                     
        if multistep == 'bdf':
            self._ms = cvode.CV_BDF
        elif multistep == 'adams':
            self._ms = cvode.CV_ADAMS
        else:
            raise ValueError
            
        if iteration == 'functional':
            self._it = cvode.CV_FUNCTIONAL
        elif iteration == 'newton':
            self._it = cvode.CV_NEWTON
        else:
            raise ValueError
            
        print self._ms, self._it
        
        self._cv = cvode.CVodeCreate(self._ms, self._it)       
        if not self._cv:
            raise MemoryError
        
        ret = cvode.CVodeSetUserData(self._cv, <void *>self)  
        
        Py_INCREF(self)
    
    def Setup(self, N_Vector y0, sun.realtype t0, reltol=0.0,
              abstol=1E-5):
        #self.y0 = y0
        #self.t0 = t0  
        cdef long N
        N = 2
        self.y0 = y0
        ret = cvode.CVodeInit(self._cv, _CvRhsFn, t0, y0._v)
        
        ret = cvode.CVodeSStolerances(self._cv, reltol, abstol);
        print ret
        #ret = cvode.CVDense(self._cv, N)
        ret = cvode.CVSpgmr(self._cv, sun.PREC_NONE, 0)
        print ret
        
    def Solve(self, sun.realtype tout, N_Vector yout):
        cdef sun.realtype *tret     
        
        tret = <sun.realtype *>malloc(cython.sizeof(sun.realtype))
        flag = cvode.CVode(self._cv, tout, yout._v, tret, cvode.CV_NORMAL)

        
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
        
    def RhsFn(self, t, y, ydot):
        raise NotImplementedError()
        
        
cdef int _CvRhsFn(sun.realtype t, sun.N_Vector y,
    		       sun.N_Vector ydot, void *user_data):
    cdef object obj
    
    obj = <object>user_data
    
    pyy = <object>y.content
    pyydot = <object>ydot.content
    return obj.RhsFn(t, pyy, pyydot)