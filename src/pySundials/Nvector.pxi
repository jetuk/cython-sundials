
"""
/*
 * -----------------------------------------------------------------
 * Functions exported by NVECTOR module
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * N_VClone
 *   Creates a new vector of the same type as an existing vector.
 *   It does not copy the vector, but rather allocates storage for
 *   the new vector.
 *
 * N_VCloneEmpty
 *   Creates a new vector of the same type as an existing vector,
 *   but does not allocate storage.
 *
 * N_VDestroy
 *   Destroys a vector created with N_VClone.
 *
 * N_VSpace
 *   Returns space requirements for one N_Vector (type 'realtype' in
 *   lrw and type 'long int' in liw).
 *
 * N_VGetArrayPointer
 *   Returns a pointer to the data component of the given N_Vector.
 *   NOTE: This function assumes that the internal data is stored
 *   as a contiguous 'realtype' array. This routine is only used in
 *   the solver-specific interfaces to the dense and banded linear
 *   solvers, as well as the interfaces to  the banded preconditioners
 *   distributed with SUNDIALS.
 *   
 * N_VSetArrayPointer
 *   Overwrites the data field in the given N_Vector with a user-supplied
 *   array of type 'realtype'.
 *   NOTE: This function assumes that the internal data is stored
 *   as a contiguous 'realtype' array. This routine is only used in
 *   the interfaces to the dense linear solver.
 *
 * N_VLinearSum
 *   Performs the operation z = a*x + b*y
 *
 * N_VConst
 *   Performs the operation z[i] = c for i = 0, 1, ..., N-1
 *
 * N_VProd
 *   Performs the operation z[i] = x[i]*y[i] for i = 0, 1, ..., N-1
 *
 * N_VDiv
 *   Performs the operation z[i] = x[i]/y[i] for i = 0, 1, ..., N-1
 *
 * N_VScale
 *   Performs the operation z = c*x
 *
 * N_VAbs
 *   Performs the operation z[i] = |x[i]| for i = 0, 1, ..., N-1
 *
 * N_VInv
 *   Performs the operation z[i] = 1/x[i] for i = 0, 1, ..., N-1
 *   This routine does not check for division by 0. It should be
 *   called only with an N_Vector x which is guaranteed to have
 *   all non-zero components.
 *
 * N_VAddConst
 *   Performs the operation z[i] = x[i] + b   for i = 0, 1, ..., N-1
 *
 * N_VDotProd
 *   Returns the dot product of two vectors:
 *         sum (i = 0 to N-1) {x[i]*y[i]}
 *
 * N_VMaxNorm
 *   Returns the maximum norm of x:
 *         max (i = 0 to N-1) ABS(x[i])
 *
 * N_VWrmsNorm
 *   Returns the weighted root mean square norm of x with weight 
 *   vector w:
 *         sqrt [(sum (i = 0 to N-1) {(x[i]*w[i])^2})/N]
 *
 * N_VWrmsNormMask
 *   Returns the weighted root mean square norm of x with weight
 *   vector w, masked by the elements of id:
 *         sqrt [(sum (i = 0 to N-1) {(x[i]*w[i]*msk[i])^2})/N]
 *   where msk[i] = 1.0 if id[i] > 0 and
 *         msk[i] = 0.0 if id[i] < 0
 *
 * N_VMin
 *   Returns the smallest element of x:
 *         min (i = 0 to N-1) x[i]
 *
 * N_VWL2Norm
 *   Returns the weighted Euclidean L2 norm of x with weight 
 *   vector w:
 *         sqrt [(sum (i = 0 to N-1) {(x[i]*w[i])^2})]
 *
 * N_VL1Norm
 *   Returns the L1 norm of x:
 *         sum (i = 0 to N-1) {ABS(x[i])}
 *
 * N_VCompare
 *   Performs the operation
 *          z[i] = 1.0 if ABS(x[i]) >= c   i = 0, 1, ..., N-1
 *                 0.0 otherwise
 *
 * N_VInvTest
 *   Performs the operation z[i] = 1/x[i] with a test for 
 *   x[i] == 0.0 before inverting x[i].
 *   This routine returns TRUE if all components of x are non-zero 
 *   (successful inversion) and returns FALSE otherwise.
 *
 * N_VConstrMask
 *   Performs the operation : 
 *       m[i] = 1.0 if constraint test fails for x[i]
 *       m[i] = 0.0 if constraint test passes for x[i]
 *   where the constraint tests are as follows:
 *      If c[i] = +2.0, then x[i] must be >  0.0.
 *      If c[i] = +1.0, then x[i] must be >= 0.0.
 *      If c[i] = -1.0, then x[i] must be <= 0.0.
 *      If c[i] = -2.0, then x[i] must be <  0.0.
 *   This routine returns a boolean FALSE if any element failed
 *   the constraint test, TRUE if all passed. It also sets a
 *   mask vector m, with elements equal to 1.0 where the
 *   corresponding constraint test failed, and equal to 0.0
 *   where the constraint test passed.
 *   This routine is specialized in that it is used only for
 *   constraint checking.
 *
 * N_VMinQuotient
 *   Performs the operation : 
 *       minq  = min ( num[i]/denom[i]) over all i such that   
 *       denom[i] != 0.
 *   This routine returns the minimum of the quotients obtained
 *   by term-wise dividing num[i] by denom[i]. A zero element
 *   in denom will be skipped. If no such quotients are found,
 *   then the large value BIG_REAL is returned.
 *
 * -----------------------------------------------------------------
 *
 * The following table lists the vector functions used by
 * different modules in SUNDIALS. The symbols in the table
 * have the following meaning:
 * S    -  called by the solver;
 * D    -  called by the dense linear solver module
 * B    -  called by the band linear solver module
 * Di   -  called by the diagonal linear solver module
 * I    -  called by the iterative linear solver module
 * BP   -  called by the band preconditioner module
 * BBDP -  called by the band-block diagonal preconditioner module
 * F    -  called by the Fortran-to-C interface
 *
 *                  ------------------------------------------------
 *                                         MODULES                  
 * NVECTOR          ------------------------------------------------
 * FUNCTIONS          CVODE/CVODES          IDA             KINSOL    
 * -----------------------------------------------------------------
 * N_VClone           S Di I                S I BBDP        S I BBDP
 * -----------------------------------------------------------------
 * N_VCloneEmpty      F                     F               F
 * -----------------------------------------------------------------
 * N_VDestroy         S Di I                S I BBDP        S I BBDP
 * -----------------------------------------------------------------
 * N_VSpace           S                     S               S         
 * -----------------------------------------------------------------
 * N_VGetArrayPointer D B BP BBDP F         D B BBDP        BBDP F     
 * -----------------------------------------------------------------
 * N_VSetArrayPointer D F                   D               F
 * -----------------------------------------------------------------
 * N_VLinearSum       S D Di I              S D I           S I       
 * -----------------------------------------------------------------
 * N_VConst           S I                   S I             I       
 * -----------------------------------------------------------------
 * N_VProd            S Di I                S I             S I       
 * -----------------------------------------------------------------
 * N_VDiv             S Di I                S I             S I
 * -----------------------------------------------------------------
 * N_VScale           S D B Di I BP BBDP    S D B I BBDP    S I BBDP  
 * -----------------------------------------------------------------
 * N_VAbs             S                     S               S         
 * -----------------------------------------------------------------
 * N_VInv             S Di                  S               S         
 * -----------------------------------------------------------------
 * N_VAddConst        S Di                  S                        
 * -----------------------------------------------------------------
 * N_VDotProd         I                     I               I         
 * -----------------------------------------------------------------
 * N_VMaxNorm         S                     S               S         
 * -----------------------------------------------------------------
 * N_VWrmsNorm        S D B I BP BBDP       S                         
 * -----------------------------------------------------------------
 * N_VWrmsNormMask                          S                         
 * -----------------------------------------------------------------
 * N_VMin             S                     S               S         
 * -----------------------------------------------------------------
 * N_VWL2Norm                                               S I       
 * -----------------------------------------------------------------
 * N_VL1Norm                                                I
 * -----------------------------------------------------------------
 * N_VCompare         Di                    S                         
 * -----------------------------------------------------------------
 * N_VInvTest         Di                                              
 * -----------------------------------------------------------------
 * N_VConstrMask                            S               S         
 * -----------------------------------------------------------------
 * N_VMinQuotient                           S               S         
 * -----------------------------------------------------------------
 */
"""

cdef class N_Vector:
    """

    """
    #cdef sun.N_Vector _v
    #cdef sun.N_Vector_Ops _ops
    #cdef long LengthRealType
    #cdef long LengthIntType
    
    def __cinit__(self, ):
        #/* Create vector */
        self._v = <sun.N_Vector>malloc(cython.sizeof(self._v))
        if not self._v:
            raise MemoryError
        #/* Create vector operation structure */
        self._ops = <sun.N_Vector_Ops>malloc(cython.sizeof(sun._generic_N_Vector_Ops))
        if not self._ops: 
            free(self._v)
            raise MemoryError

        self._v.ops = self._ops
        self._v.content = <void *>self
        
        self._ops.nvclone           = _NVClone
        self._ops.nvcloneempty      = _NVCloneEmpty
        self._ops.nvdestroy         = _NVDestroy
        self._ops.nvspace           = _NVSpace
        self._ops.nvgetarraypointer = _N_VGetArrayPointer
        self._ops.nvsetarraypointer = _N_VSetArrayPointer
        self._ops.nvlinearsum       = _NVLinearSum
        self._ops.nvconst           = _NVConstant
        self._ops.nvprod            = _NVProd
        self._ops.nvdiv             = _NVDiv
        self._ops.nvscale           = _NVScale
        self._ops.nvabs             = _NVAbs
        self._ops.nvinv             = _NVInv
        self._ops.nvaddconst        = _NVAddConst
        self._ops.nvdotprod         = _NVDotProd
        self._ops.nvmaxnorm         = _NVMaxNorm
        self._ops.nvwrmsnormmask    = _NVWrmsNormMask
        self._ops.nvwrmsnorm        = _NVWrmsNorm
        self._ops.nvmin             = _NVMin        
        self._ops.nvwl2norm         = _NVWL2Norm
        self._ops.nvl1norm          = _NVL1Norm
        self._ops.nvcompare         = _NVCompare
        self._ops.nvinvtest         = _NVInvTest
        self._ops.nvconstrmask      = _NVConstrMask
        self._ops.nvminquotient     = _NVMinQuotient
        
        Py_INCREF(self)
        
    def __dealloc__(self, ):    
        
        if not self._ops:
            free(self._ops)
        if not self._v:
            free(self._v)       
            
        
        
    def Clone(self, ):
        """
         * N_VClone
         *   Creates a new vector of the same type as an existing vector.
         *   It does not copy the vector, but rather allocates storage for
         *   the new vector.
        """
        raise NotImplementedError()
    
    def CloneEmpty(self, ):
        """
         * N_VCloneEmpty
         *   Creates a new vector of the same type as an existing vector,
         *   but does not allocate storage.        
        """
        raise NotImplementedError()
        
    def Destroy(self, ):
        """
         * N_VDestroy
         *   Destroys a vector created with N_VClone.
        """
        Py_DECREF(self)
        raise NotImplementedError()
        
    def GetArrayPointer(self):
        """
         * N_VGetArrayPointer
         *   Returns a pointer to the data component of the given N_Vector.
         *   NOTE: This function assumes that the internal data is stored
         *   as a contiguous 'realtype' array. This routine is only used in
         *   the solver-specific interfaces to the dense and banded linear
         *   solvers, as well as the interfaces to  the banded preconditioners
         *   distributed with SUNDIALS.
        """
        raise NotImplementedError()
         
    def SetArrayPointer(self, v_data):
        """
         * N_VSetArrayPointer
         *   Overwrites the data field in the given N_Vector with a user-supplied
         *   array of type 'realtype'.
         *   NOTE: This function assumes that the internal data is stored
         *   as a contiguous 'realtype' array. This routine is only used in
         *   the interfaces to the dense linear solver.    
        """
        raise NotImplementedError()
        
    def LinearSum(self, double a, double b, N_Vector y, N_Vector z):
        """
         * N_VLinearSum
         *   Performs the operation z = a*x + b*y        
        """
        raise NotImplementedError()
        
    def Constant(self, c):
        """
        * N_VConst
        *   Performs the operation z[i] = c for i = 0, 1, ..., N-1        
        """
        raise NotImplementedError()
        
    def Prod(self, y, z):
        """
         * N_VProd
         *   Performs the operation z[i] = x[i]*y[i] for i = 0, 1, ..., N-1
        """
        raise NotImplementedError()
        
    def Div(self, y, z):
        """
         * N_VDiv
         *   Performs the operation z[i] = x[i]/y[i] for i = 0, 1, ..., N-1
        """
        raise NotImplementedError()    
        
    def Scale(self, c, z):
        """
         * N_VScale
         *   Performs the operation z = c*x       
        """
        raise NotImplementedError()
        
    def Abs(self, z):
        """
         * N_VAbs
         *   Performs the operation z[i] = |x[i]| for i = 0, 1, ..., N-1
        """
        raise NotImplementedError()
        
        
    def Inv(self, z):
        """
         * N_VInv
         *   Performs the operation z[i] = 1/x[i] for i = 0, 1, ..., N-1
         *   This routine does not check for division by 0. It should be
         *   called only with an N_Vector x which is guaranteed to have
         *   all non-zero components.
        """
        raise NotImplementedError()
        
    def AddConst(self, b, z):
        """
         * N_VAddConst
         *   Performs the operation z[i] = x[i] + b   for i = 0, 1, ..., N-1
        """
        raise NotImplementedError()
        
    def DotProd(self, y):
        """
         * N_VDotProd
         *   Returns the dot product of two vectors:
         *         sum (i = 0 to N-1) {x[i]*y[i]}
        """
        raise NotImplementedError()
        
    def MaxNorm(self, ):
        """
         * N_VMaxNorm
         *   Returns the maximum norm of x:
         *         max (i = 0 to N-1) ABS(x[i])        
        """
        raise NotImplementedError()
        
    def WrmsNorm(self, N_Vector w):
        """
         * N_VWrmsNorm
         *   Returns the weighted root mean square norm of x with weight 
         *   vector w:
         *         sqrt [(sum (i = 0 to N-1) {(x[i]*w[i])^2})/N]
        """
        raise NotImplementedError()
        
    def WrmsNormMask(self, w, mask):
        """
         * N_VWrmsNormMask
         *   Returns the weighted root mean square norm of x with weight
         *   vector w, masked by the elements of id:
         *         sqrt [(sum (i = 0 to N-1) {(x[i]*w[i]*msk[i])^2})/N]
         *   where msk[i] = 1.0 if id[i] > 0 and
         *         msk[i] = 0.0 if id[i] < 0        
        """
        raise NotImplementedError()
        
    def Min(self, ):
        """
         * N_VMin
         *   Returns the smallest element of x:
         *         min (i = 0 to N-1) x[i]          
        """
        raise NotImplementedError()
        
    def WL2Norm(self, w):
        """
         * N_VWL2Norm
         *   Returns the weighted Euclidean L2 norm of x with weight 
         *   vector w:
         *         sqrt [(sum (i = 0 to N-1) {(x[i]*w[i])^2})]
        """
        raise NotImplementedError()
        
    def L1Norm(self, ):
        """
         * N_VL1Norm
         *   Returns the L1 norm of x:
         *         sum (i = 0 to N-1) {ABS(x[i])}                
        """
        raise NotImplementedError()
        
        
    def Compare(self, c, z):
        """
         * N_VCompare
         *   Performs the operation
         *          z[i] = 1.0 if ABS(x[i]) >= c   i = 0, 1, ..., N-1
         *                 0.0 otherwise
        """
        raise NotImplementedError()
        
    def Invtest(self, z):
        """
         * N_VInvTest
         *   Performs the operation z[i] = 1/x[i] with a test for 
         *   x[i] == 0.0 before inverting x[i].
         *   This routine returns TRUE if all components of x are non-zero 
         *   (successful inversion) and returns FALSE otherwise.
        """
        raise NotImplementedError()
        
    def ConstrMask(self, c, m):
        """
         * N_VConstrMask
         *   Performs the operation : 
         *       m[i] = 1.0 if constraint test fails for x[i]
         *       m[i] = 0.0 if constraint test passes for x[i]
         *   where the constraint tests are as follows:
         *      If c[i] = +2.0, then x[i] must be >  0.0.
         *      If c[i] = +1.0, then x[i] must be >= 0.0.
         *      If c[i] = -1.0, then x[i] must be <= 0.0.
         *      If c[i] = -2.0, then x[i] must be <  0.0.
         *   This routine returns a boolean FALSE if any element failed
         *   the constraint test, TRUE if all passed. It also sets a
         *   mask vector m, with elements equal to 1.0 where the
         *   corresponding constraint test failed, and equal to 0.0
         *   where the constraint test passed.
         *   This routine is specialized in that it is used only for
         *   constraint checking.
        """
        raise NotImplementedError()
        
    def MinQuotient(self, denom):
        """
         * N_VMinQuotient
         *   Performs the operation : 
         *       minq  = min ( num[i]/denom[i]) over all i such that   
         *       denom[i] != 0.
         *   This routine returns the minimum of the quotients obtained
         *   by term-wise dividing num[i] by denom[i]. A zero element
         *   in denom will be skipped. If no such quotients are found,
         *   then the large value BIG_REAL is returned.        
        """
        raise NotImplementedError()

        
        
#        self._ops.nvcloneempty      = N_VCloneEmpty_MyNvector
#        self._ops.nvdestroy         = N_VDestroy_MyNvector
#        self._ops.nvspace           = N_VSpace_MyNvector
#    #    ops->nvgetarraypointer = N_VGetArrayPointer_Serial;
#    #    ops->nvsetarraypointer = N_VSetArrayPointer_Serial;
#        self._ops.nvlinearsum       = N_VLinearSum_MyNvector
#        self._ops.nvconst           = N_VConst_MyNvector
#        self._ops.nvprod            = N_VProd_MyNvector
#        self._ops.nvdiv             = N_VDiv_MyNvector
#        self._ops.nvscale           = N_VScale_MyNvector
#        self._ops.nvabs             = N_VAbs_MyNvector
#        self._ops.nvinv             = N_VInv_MyNvector
#        self._ops.nvaddconst        = N_VAddConst_MyNvector
#        self._ops.nvdotprod         = N_VDotProd_MyNvector
#    #    ops->nvmaxnorm         = N_VMaxNorm_Serial;
#    #    ops->nvwrmsnormmask    = N_VWrmsNormMask_Serial;
#    #    ops->nvwrmsnorm        = N_VWrmsNorm_Serial;
#    #    ops->nvmin             = N_VMin_Serial;
#    #    ops->nvwl2norm         = N_VWL2Norm_Serial;
#    #    ops->nvl1norm          = N_VL1Norm_Serial;
#    #    ops->nvcompare         = N_VCompare_Serial;
#    #    ops->nvinvtest         = N_VInvTest_Serial;
#    #    ops->nvconstrmask      = N_VConstrMask_Serial;
#    #    ops->nvminquotient     = N_VMinQuotient_Serial;
#    
cdef sun.N_Vector _NVClone(sun.N_Vector w):
    cdef N_Vector pyw
    cdef N_Vector pyz
    
    pyw = <object>w.content
    pyz = pyw.Clone()
    return pyz._v
    
cdef sun.N_Vector _NVCloneEmpty(sun.N_Vector w):
    cdef N_Vector pyw
    cdef N_Vector pyz
    
    pyw = <object>w.content
    pyz = pyw.CloneEmpty()
    return pyz._v
    
cdef void _NVDestroy(sun.N_Vector w):
    cdef N_Vector pyw
    
    pyw = <object>w.content
    pyw.Destroy()
    
cdef void _NVSpace(sun.N_Vector w, long *lir, long *liw):
    cdef N_Vector pyw
    
    pyw = <object>w.content
    lir = <long *>pyw.LengthRealType
    liw = <long *>pyw.LengthIntType
    
    
cdef sun.realtype* _N_VGetArrayPointer(sun.N_Vector v):
    cdef N_Vector pyv    
    cdef uintptr_t res
    cdef sun.realtype *_res
    pyv = <object>v.content
    res = pyv.GetArrayPointer()
    _res = <sun.realtype *>res
    return _res


cdef void _N_VSetArrayPointer(sun.realtype *v_data, sun.N_Vector v ):
    cdef N_Vector pyv    
    cdef uintptr_t ptr
    pyv = <object>v.content
    ptr = <uintptr_t>v_data
    pyv.SetArrayPointer(ptr)


cdef void _NVLinearSum(sun.realtype a, sun.N_Vector x,
                       sun.realtype b, sun.N_Vector y, sun.N_Vector z):
    cdef N_Vector pyx
    cdef N_Vector pyy
    cdef N_Vector pyz
    
    pyx = <N_Vector>x.content
    pyy = <N_Vector>y.content
    pyz = <N_Vector>z.content
    
    pyx.LinearSum(a, b, pyy, pyz)
    
cdef void _NVConstant(sun.realtype c, sun.N_Vector z):
    cdef N_Vector pyz
    
    pyz = <object>z.content
    pyz.Constant(c)
    
cdef void _NVProd(sun.N_Vector x, sun.N_Vector y, sun.N_Vector z):
    cdef N_Vector pyx
    cdef N_Vector pyy
    cdef N_Vector pyz
    
    pyx = <object>x.content
    pyy = <object>y.content
    pyz = <object>z.content
    pyx.Prod(pyy,pyz)
        
cdef void _NVDiv(sun.N_Vector x, sun.N_Vector y, sun.N_Vector z):
    cdef N_Vector pyx
    cdef N_Vector pyy
    cdef N_Vector pyz
    
    pyx = <object>x.content
    pyy = <object>y.content
    pyz = <object>z.content
    pyx.Div(pyy,pyz)
    
cdef void _NVScale(sun.realtype c, sun.N_Vector x, sun.N_Vector z):
    cdef N_Vector pyx
    cdef N_Vector pyz
    
    pyx = <object>x.content
    pyz = <object>z.content
    pyx.Scale(c,pyz)
    
    
cdef void _NVAbs(sun.N_Vector x, sun.N_Vector z):
    cdef N_Vector pyx
    cdef N_Vector pyz
    
    pyx = <object>x.content
    pyz = <object>z.content
    pyx.Abs(pyz)
    
cdef void _NVInv(sun.N_Vector x, sun.N_Vector z):
    cdef N_Vector pyx
    cdef N_Vector pyz
    
    pyx = <object>x.content
    pyz = <object>z.content
    
    pyx.Inv(pyz)

cdef void _NVAddConst(sun.N_Vector x, sun.realtype b, sun.N_Vector z):
    cdef N_Vector pyx
    cdef N_Vector pyz
    
    pyx = <object>x.content
    pyz = <object>z.content
    
    pyx.AddConst(b, pyz)

cdef sun.realtype _NVDotProd(sun.N_Vector x, sun.N_Vector y):
    cdef N_Vector pyx
    cdef N_Vector pyy
    cdef sun.realtype res
    
    pyx = <object>x.content
    pyy = <object>y.content
    
    res = pyx.DotProd(pyy)
    return res
    
cdef sun.realtype _NVMaxNorm(sun.N_Vector x):
    cdef N_Vector pyx
    cdef sun.realtype res
    
    pyx = <object>x.content
    
    res = pyx.MaxNorm()
    return res    
    
cdef sun.realtype _NVWrmsNorm(sun.N_Vector x, sun.N_Vector w):
    cdef N_Vector pyx
    cdef N_Vector pyw
    cdef sun.realtype res
    
    pyx = <object>x.content
    pyw = <object>w.content
    
    res = pyx.WrmsNorm(pyw)
    return res
    
cdef sun.realtype _NVWrmsNormMask(sun.N_Vector x, sun.N_Vector w, sun.N_Vector id):
    cdef N_Vector pyx
    cdef N_Vector pyw
    cdef N_Vector pyid
    cdef sun.realtype res
    
    pyx = <object>x.content
    pyw = <object>w.content
    pyid = <object>id.content
    
    res = pyx.WrmsNormMask(pyw, pyid)
    return res
    
cdef sun.realtype _NVMin(sun.N_Vector x):
    cdef N_Vector pyx
    cdef sun.realtype res
    
    pyx = <object>x.content
    res = pyx.Min()
    return res

cdef sun.realtype _NVWL2Norm(sun.N_Vector x, sun.N_Vector w):
    cdef N_Vector pyx
    cdef N_Vector pyw
    cdef sun.realtype res
    
    pyx = <object>x.content
    pyw = <object>w.content
    
    res = pyx.WL2Norm(pyw)
    return res    
    
cdef sun.realtype _NVL1Norm(sun.N_Vector x):
    cdef N_Vector pyx
    cdef sun.realtype res
    
    pyx = <object>x.content
    res = pyx.L1Norm()
    return res
    

cdef void _NVCompare(sun.realtype c, sun.N_Vector x, sun.N_Vector z):
    cdef N_Vector pyx
    cdef N_Vector pyz
           
    pyx = <object>x.content
    pyz = <object>z.content
    
    pyx.Compare(c, pyz)
    
    
cdef sun.booleantype _NVInvTest(sun.N_Vector x, sun.N_Vector z):
    cdef N_Vector pyx
    cdef N_Vector pyz
    cdef sun.booleantype res
        
    pyx = <object>x.content
    pyz = <object>z.content
    
    res = pyx.InvTest(pyz)    
    return res
    
cdef sun.booleantype _NVConstrMask(sun.N_Vector c, sun.N_Vector x, sun.N_Vector m):
    cdef N_Vector pyx
    cdef N_Vector pyc
    cdef N_Vector pym
    cdef sun.booleantype res
        
    pyx = <object>x.content
    pyc = <object>c.content
    pym = <object>m.content
    
    res = pyx.ConstrMask(pyc, pym) 
    return res    

cdef sun.realtype _NVMinQuotient(sun.N_Vector num, sun.N_Vector denom):
    cdef N_Vector pynum
    cdef N_Vector pydenom
    cdef sun.realtype res
    
    pynum = <object>num.content
    pydenom = <object>denom.content
    
    res = pynum.MinQuotient(pydenom)
    return res   
    
    
include "NvectorImplementations/NvectorNumpy.pxi"
include "NvectorImplementations/NvectorMemoryView.pxi"
include "NvectorImplementations/NvectorVariables.pxi"