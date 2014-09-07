# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 20:14:11 2014

@author: james
"""

def propertyGenerator(name, variables, module, mem_ptr, 
                      getter=None, setter=None, indent=1,
                      getter_success=None, setter_success=None,
                      getter_errors={}, setter_errors={}, 
                      exception="ValueError", 
                      unknown_error_msg="Unknown error"):
    
    """ 
    Function to generate cython code for a series of similar properties

    Example:    
    
    property maxSetupCalls:
        def __get__(self, ):
            raise NotImplementedError()
        
        def __set__(self, int msbset):
            ret = kinsol.KINSetMaxSetupCalls(self._kn, msbset)
            if ret == kinsol.KIN_SUCCESS:
                return
            if ret == kinsol.KIN_MEM_NULL:
                raise ValueError('Setup first must be called before SetMaxSetupCalls')
            if ret == kinsol.KIN_ILL_INPUT:
                raise ValueError('Illegal value')
            raise ValueError('Unknown error ({}))'.format(ret))
    
    """
    fmt = lambda s,i:" "*i*4 + "{}\n".format(s)
    make_args = lambda f,vs,d:d.join([f.format(v) for v in vs])    
    
    vnames = [vn for vn,vt in variables]
    
    s = fmt("property {}:".format(name), indent)
    s += fmt("def __get__(self, ):",indent+1)
    
    if getter is None:
        s += fmt("raise NotImplementedError()", indent+2)
    else:
        for vname,vtype in variables:
            s += fmt("cdef {} {}".format(vtype,vname), indent+2)
            
        args = make_args("&{}", vnames, ", ")
        
        s += fmt("ret = {}.{}({},{})".format(module,getter,mem_ptr,args),indent+2 )
        
        s += fmt("if ret == {}.{}:".format(module,getter_success), indent+2)
        
        args = make_args("{}", vnames, ", ")        
        s += fmt("return {}".format(args), indent+3)

        # Write all errors
        for error,msg in getter_errors.iteritems():       
            s += fmt("if ret == {}.{}:".format(module, error), indent+2)
            s += fmt('raise {}("{}")'.format(exception, msg), indent+3)
        
        # Final catch all error
        s += fmt('raise {}("{}")'.format(exception, unknown_error_msg), indent+2)
        
    s += fmt("",indent+1)
    
    s += fmt("def __set__(self, value):", indent+1)
    if setter is None:    
        s += fmt("raise NotImplementedError()", indent+2)
    else:
        for i,(vname,vtype) in enumerate(variables):
            if len(variables) > 1: # Ensure the tuple is broken up if more than 1 variable
                s += fmt("cdef {} {} = value[{}]".format(vtype,vname,i), indent+2)
            else:
                s += fmt("cdef {} {} = value".format(vtype,vname), indent+2)
            
        args = make_args("{}", vnames, ", ")
        
        s += fmt("ret = {}.{}({},{})".format(module,setter,mem_ptr,args), indent+2)
        
        s += fmt("if ret == {}.{}:".format(module,setter_success), indent+2)
        
        s += fmt("return", indent+3)

        # Write all errors
        for error,msg in setter_errors.iteritems():       
            s += fmt("if ret == {}.{}:".format(module, error), indent+2)
            s += fmt('raise {}("{}")'.format(exception, msg), indent+3)
        
        # Final catch all error
        s += fmt('raise {}("{}")'.format(exception, unknown_error_msg), indent+2)
            
    return s
    
    
def createKinsolProperties(filename ):
    
    module = 'kinsol'
    mem_ptr = 'self._kn'


    
    setter_success = 'KIN_SUCCESS'
    setter_errors={
        'KIN_MEM_NULL':"KINSOL memory pointer is NULL",
        'KIN_ILL_INPUT':"Illegal value"}   
    
    properties = (
        ('printLevel',[('printfl','int'),],'KINSetPrintLevel'),
        ('numMaxIters',[('mxiter','long int'),], 'KINSetNumMaxIters'),
        ('noInitSeutp',[('noInitSeutp','sun.booleantype'),],'KINSetNoInitSetup'),
        ('noResMon',[('noNNIResMon','sun.booleantype'),],'KINSetNoResMon'),
        ('maxSetupCalls',[('msbset','long int'),], 'KINSetMaxSetupCalls'),
        ('maxSubSetupCalls',[('msbsetsub','long int'),], 'KINSetMaxSubSetupCalls'),
        ('etaForm',[('etachoice','int'),], 'KINSetEtaForm'),
        ('etaConstValue',[('eta','sun.realtype'),], 'KINSetEtaConstValue'),
        ('etaParams',[('egamma','sun.realtype'),('ealpha','sun.realtype')],'KINSetEtaParams'),
        ('resMonParams',[('omegamin','sun.realtype'),('omegamax','sun.realtype')],'KINSetResMonParams'),
        ('resMonConstValue',[('omegaconst','sun.realtype'),],'KINSetResMonConstValue'),
        ('noMinEps',[('noMinEps','sun.booleantype'),],'KINSetNoMinEps'),
        ('maxNewtonStep',[('mxnewtstep','sun.realtype'),],'KINSetMaxNewtonStep'),
        ('maxBetaFails',[('mxbcf','long int'),],'KINSetMaxBetaFails'),
        ('funcNormTol',[('fnormtol','sun.realtype'),],'KINSetFuncNormTol'),
        ('scaledStepTol',[('scsteptol','sun.realtype'),],'KINSetScaledStepTol'),
    )

    s = """

cdef class BaseKinsol:
    
    def __cinit__(self,*args,**kwds): 
                     
        self._kn = kinsol.KINCreate()
        if not self._kn:
            raise MemoryError
        
        ret = kinsol.KINSetUserData(self._kn, <void *>self)  
        if ret != 0:
            raise KinsolError()    


    """
    s += "\n"
    for name,variables,func in properties:
        s += propertyGenerator(name,variables,module,mem_ptr,setter=func,
                               setter_success=setter_success,setter_errors=setter_errors)
        s += "\n"
    
    
    
    fh = open(filename, 'wb')
    fh.write(s)
    fh.close()
    
def createCvodeProperties(filename, ):
    module = 'cvode'
    mem_ptr = 'self._cv'
    
    getter_success = 'CV_SUCCESS'
    getter_errors = {
        'CV_MEM_NULL':'CVODE memory pointer is NULL',
        'CV_NO_SLDET':'Stability limit was not turned on',
    }   
    
    get_properties = (    
        ('workSpace', [('lenrw', 'long int'),('leniw','long int')],'CVodeGetWorkSpace'),
        ('numSteps', [('nsteps', 'long int'),], 'CVodeGetNumSteps'),
        ('numRhsEvals', [('nfevals', 'long int'),], 'CVodeGetNumRhsEvals' ),
        ('numLinSolvSetups', [('nlinsetups', 'long int'),], 'CVodeGetNumLinSolvSetups'),
        ('numErrTestFails', [('netfails', 'long int'),], 'CVodeGetNumErrTestFails'),
        ('lastOrder', [('qlast', 'int'),], 'CVodeGetLastOrder'),
        ('currentOrder', [('qcur', 'int'),], 'CVodeGetCurrentOrder'),
        ('numStabLimOrderReds', [('nslred', 'long int'),], 'CVodeGetNumStabLimOrderReds'),
        ('actualInitStep', [('hinused', 'sun.realtype'),], 'CVodeGetActualInitStep'),
        ('lastStep', [('hinused', 'sun.realtype'),], 'CVodeGetLastStep'),
        ('currentStep', [('hcur', 'sun.realtype'),], 'CVodeGetCurrentStep'),
        ('currentTime', [('tcur', 'sun.realtype'),], 'CVodeGetCurrentTime'),
        ('tolScaleFactor', [('tolsfac', 'sun.realtype'),], 'CVodeGetTolScaleFactor'),
        ('numGEvals', [('ngevals', 'long int'),], 'CVodeGetNumGEvals'),
        ('numNonlinSolvIters', [('nniters', 'long int'),], 'CVodeGetNumNonlinSolvIters'),
        ('numNonlinSolvConvFails', [('nncfails', 'long int'),], 'CVodeGetNumNonlinSolvConvFails'),
        ('dlsWorkSpace', [('lenrwLS', 'long int'),('leniwLS', 'long int'),], 'CVDlsGetWorkSpace'),
        ('dlsNumJacEvals', [('njevals', 'long int'),], 'CVDlsGetNumJacEvals'),
        ('dlsNumRhsEvals', [('nfevalsLS', 'long int'),], 'CVDlsGetNumRhsEvals'),
        ('dlsLastFlag', [('flag', 'long int'),], 'CVDlsGetLastFlag'), 
        ('spilsWorkSpace', [('lenrwLS', 'long int'),('leniwLS', 'long int')], 'CVSpilsGetWorkSpace'),
        ('spilsNumPrecEvals', [('npevals', 'long int'),], 'CVSpilsGetNumPrecEvals'),
        ('spilsNumPrecSolves', [('npsolves', 'long int'),], 'CVSpilsGetNumPrecSolves'),
        ('spilsNumLinIters', [('nliters', 'long int'),], 'CVSpilsGetNumLinIters'),
        ('spilsNumConvFails', [('nlcfails', 'long int'),], 'CVSpilsGetNumConvFails'),
        ('spilsNumJtimesEvals', [('njvevals', 'long int'),], 'CVSpilsGetNumJtimesEvals'),
        ('spilsNumRhsEvals', [('nfevalsLS', 'long int'),], 'CVSpilsGetNumRhsEvals')
    
    )
        
    
    setter_success = 'CV_SUCCESS'
    setter_errors={
        'CV_MEM_NULL':"CVODE memory pointer is NULL",
        'CV_ILL_INPUT':"Illegal value"}   
    
    set_properties = (
        ('maxNumSteps', [('mxsteps','long int'),], 'CVodeSetMaxNumSteps'),
        ('maxOrd', [('maxord','int'),], 'CVodeSetMaxOrd'),
        ('maxHnilWarns', [('mxhnil','int'),], 'CVodeSetMaxHnilWarns'),
        ('stabLimDet', [('stldet', 'sun.booleantype'),], 'CVodeSetStabLimDet'),
        ('initStep', [('hin','sun.realtype'),], 'CVodeSetInitStep'),
        ('minStep', [('hmin', 'sun.realtype'),], 'CVodeSetMinStep'),
        ('maxStep', [('hmax', 'sun.realtype'),], 'CVodeSetMaxStep'),
        ('stopTime', [('tstop', 'sun.realtype'),], 'CVodeSetStopTime'),
        ('maxErrTestFails', [('maxnef', 'int'),], 'CVodeSetMaxErrTestFails'),
        ('maxNonlinIters', [('maxcor', 'int'),], 'CVodeSetMaxNonlinIters'),
        ('maxConvFails', [('maxncf', 'int'),], 'CVodeSetMaxConvFails'),
        ('nonlinConvCoef', [('nlscoef', 'sun.realtype'),], 'CVodeSetNonlinConvCoef'),
        ('iterType', [('iter', 'int'),], 'CVodeSetIterType'),

    )


    s = """

cdef class BaseCvode:
    
    def __cinit__(self, *args, **kwds): 
        multistep = kwds.pop('multistep', 'bdf')
        iteration = kwds.pop('iteration', 'functional')
                     
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
            
        
        self._cv = cvode.CVodeCreate(self._ms, self._it)       
        if not self._cv:
            raise MemoryError
        
        ret = cvode.CVodeSetUserData(self._cv, <void *>self)   


    """
    s += "\n"
    
    for name,variables,func in get_properties:
        s += propertyGenerator(name,variables,module,mem_ptr,getter=func,
                               getter_success=getter_success,getter_errors=getter_errors)
        s += "\n"    
    
    for name,variables,func in set_properties:
        s += propertyGenerator(name,variables,module,mem_ptr,setter=func,
                               setter_success=setter_success,setter_errors=setter_errors)
        s += "\n"
    
    
    
    fh = open(filename, 'wb')
    fh.write(s)
    fh.close()
    
if __name__ == '__main__':
    
    s = propertyGenerator('maxSetupCalls', [('msbset','int'),], 'kinsol', 'self._kn',
                        setter='KINSetMaxSetupCalls', setter_success='KIN_SUCCESS',
                        setter_errors={'KIN_MEM_NULL':"KINSOL memory pointer is NULL",
                                        'KIN_ILL_INPUT':"Illegal value"},)
                                        
                                        
    