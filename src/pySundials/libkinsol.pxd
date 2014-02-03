from libsundials cimport *

cdef extern from "kinsol/kinsol.h":
    

#    /*
#     * =================================================================
#     *              K I N S O L     C O N S T A N T S
#     * =================================================================
#     */
#    
#    /*
#     * -----------------------------------------------------------------
#     * KINSOL return flags
#     * -----------------------------------------------------------------
#     */

    enum: KIN_SUCCESS
    enum: KIN_INITIAL_GUESS_OK
    enum: KIN_STEP_LT_STPTOL
    
    enum: KIN_WARNING
    
    enum: KIN_MEM_NULL
    enum: KIN_ILL_INPUT
    enum: KIN_NO_MALLOC
    enum: KIN_MEM_FAIL
    enum: KIN_LINESEARCH_NONCONV
    enum: KIN_MAXITER_REACHED
    enum: KIN_MXNEWT_5X_EXCEEDED
    enum: KIN_LINESEARCH_BCFAIL
    enum: KIN_LINSOLV_NO_RECOVERY
    enum: KIN_LINIT_FAIL
    enum: KIN_LSETUP_FAIL
    enum: KIN_LSOLVE_FAIL
    
    enum: KIN_SYSFUNC_FAIL
    enum: KIN_FIRST_SYSFUNC_ERR
    enum: KIN_REPTD_SYSFUNC_ERR
    
    
#    /*
#     * -----------------------------------------------------------------
#     * enum:eration for inputs to KINSetEtaForm (eta choice)
#     * -----------------------------------------------------------------
#     * KIN_ETACONSTANT : use constant value for eta (default value is
#     *                   0.1 but a different value can be specified via
#     *                   a call to KINSetEtaConstValue)
#     *
#     * KIN_ETACHOICE1 : use choice #1 as given in Eisenstat and Walker's
#     *                  paper of SIAM J.Sci.Comput.,17 (1996), pp 16-32,
#     *                  wherein eta is defined to be:
#     *
#     *              eta(k+1) = ABS(||F(u_k+1)||_L2-||F(u_k)+J(u_k)*p_k||_L2)
#     *                       ---------------------------------------------
#     *                                       ||F(u_k)||_L2
#     *
#     *                                                      1+sqrt(5)
#     *              eta_safe = eta(k)^ealpha where ealpha = ---------
#     *                                                          2
#     *
#     * KIN_ETACHOICE2 : use choice #2 as given in Eisenstat and Walker's
#     *                  paper wherein eta is defined to be:
#     *
#     *                                  [ ||F(u_k+1)||_L2 ]^ealpha
#     *              eta(k+1) = egamma * [ --------------- ]
#     *                                  [  ||F(u_k)||_L2  ]
#     *
#     *                  where egamma = [0,1] and ealpha = (1,2]
#     *
#     *              eta_safe = egamma*(eta(k)^ealpha)
#     *
#     *                  Note: The default values of the scalar
#     *                  coefficients egamma and ealpha (both required)
#     *                  are egamma = 0.9 and ealpha = 2.0, but the
#     *                  routine KINSetEtaParams can be used to specify
#     *                  different values.
#     *
#     * When using either KIN_ETACHOICE1 or KIN_ETACHOICE2, if
#     * eta_safe > 0.1 then the following safeguard is applied:
#     *
#     *  eta(k+1) = MAX {eta(k+1), eta_safe}
#     *
#     * The following safeguards are always applied when using either
#     * KIN_ETACHOICE1 or KIN_ETACHOICE2 so that eta_min <= eta <= eta_max:
#     *
#     *  eta(k+1) = MAX {eta(k+1), eta_min}
#     *  eta(k+1) = MIN {eta(k+1), eta_max}
#     *
#     * where eta_min = 1.0e-4 and eta_max = 0.9 (see KINForcingTerm).
#     * -----------------------------------------------------------------
#     */
    
    enum: KIN_ETACHOICE1
    enum: KIN_ETACHOICE2
    enum: KIN_ETACONSTANT
      
#    /*
#     * -----------------------------------------------------------------
#     * enum:eration for global strategy
#     * -----------------------------------------------------------------
#     * Choices are KIN_NONE and KIN_LINESEARCH.
#     * -----------------------------------------------------------------
#     */
      
    enum: KIN_NONE
    enum: KIN_LINESEARCH
    
#    /*
#     * =================================================================
#     *              F U N C T I O N   T Y P E S
#     * =================================================================
#     */
#    
#    /*
#     * -----------------------------------------------------------------
#     * Type : KINSysFn
#     * -----------------------------------------------------------------
#     * The user-supplied subroutine implementing the nonlinear system
#     * function (vector-valued function) F must take as input the
#     * dependent variable vector uu (type N_Vector), and set fval (type
#     * N_Vector) equal to F(uu) before returning. Additional workspace
#     * is allocated by the user and referenced by the user_data memory
#     * pointer.
#     * 
#     * Note: The user-defined routine (internally referenced by a
#     * a pointer (type KINSysFn) named func) should have an 'int' return
#     * value type. However, the return value is currently ignored.
#     * -----------------------------------------------------------------
#     */
    
    ctypedef int (*KINSysFn)(N_Vector uu, N_Vector fval, void *user_data )
    
    
#    /*
#     * -----------------------------------------------------------------
#     * Type : KINErrHandlerFn
#     * -----------------------------------------------------------------
#     * A function eh, which handles error messages, must have type
#     * KINErrHandlerFn.
#     * The function eh takes as input the error code, the name of the
#     * module reporting the error, the error message, and a pointer to
#     * user data, the same as that passed to KINSetUserData.
#     * 
#     * All error codes are negative, except KIN_WARNING which indicates 
#     * a warning (the solver continues).
#     *
#     * A KINErrHandlerFn has no return value.
#     * -----------------------------------------------------------------
#     */
    
    ctypedef void (*KINErrHandlerFn)(int error_code, 
    				const char *module, const char *function, 
    				char *msg, void *user_data)
    
    
#    /*
#     * -----------------------------------------------------------------
#     * Type : KINInfoHandlerFn
#     * -----------------------------------------------------------------
#     * A function ih, which handles info messages, must have type
#     * KINInfoHandlerFn.
#     * The function ih takes as input the name of the module and of the
#     * function reporting the info message and a pointer to
#     * user data, the same as that passed to KINSetfdata.
#     * 
#     * A KINInfoHandlerFn has no return value.
#     * -----------------------------------------------------------------
#     */
    
    ctypedef void (*KINInfoHandlerFn)(const char *module, const char *function, 
    				 char *msg, void *user_data)
    
#    /*
#     * ================================================================
#     *          U S E R - C A L L A B L E   R O U T I N E S           
#     * ================================================================
#     */
#    
#    /*
#     * -----------------------------------------------------------------
#     * Function : KINCreate
#     * -----------------------------------------------------------------
#     * KINCreate allocates and initializes an internal memory block for
#     * the KINSOL solver module.
#     *
#     * If successful, KINCreate returns a pointer to the initialized
#     * memory block which should be passed to KINInit. If an
#     * error occurs, then KINCreate returns a NULL pointer.
#     * -----------------------------------------------------------------
#     */
    
    void *KINCreate()
    
#    /*
#     * -----------------------------------------------------------------
#     * Optional Input Specification Functions (KINSOL)
#     * -----------------------------------------------------------------
#     * The following functions can be called to set optional inputs:
#     *
#     *     Function Name      |    Optional Input  [Default Value]
#     *                        |
#     * -----------------------------------------------------------------
#     *                        |
#     * KINSetErrHandlerFn     | user-provided ErrHandler function.
#     *                        | [internal]
#     *                        |
#     * KINSetErrFile          | pointer (type FILE) indicating where all
#     *                        | warning/error messages should be sent
#     *                        | if the default internal error handler 
#     *                        | is used
#     *                        | [stderr]
#     *                        |
#     * KINSetPrintLevel       | level of verbosity of output:
#     *                        |
#     *                        |  0  no statistical information is
#     *                        |     displayed (default level)
#     *                        |
#     *                        |  1  for each nonlinear iteration display
#     *                        |     the following information: the scaled
#     *                        |     norm (L2) of the system function
#     *                        |     evaluated at the current iterate, the
#     *                        |     scaled norm of the Newton step (only if
#     *                        |     using KIN_NONE), and the
#     *                        |     number of function evaluations performed
#     *                        |     thus far
#     *                        |
#     *                        |  2  display level 1 output and the
#     *                        |     following values for each iteration:
#     *                        |
#     *                        |       fnorm (L2) = ||fscale*func(u)||_L2
#     *                        |       (only for KIN_NONE)
#     *                        |
#     *                        |       scaled fnorm (for stopping) =
#     *                        |       ||fscale*ABS(func(u))||_L-infinity
#     *                        |       (for KIN_NONE and
#     *                        |       KIN_LINESEARCH)
#     *                        |
#     *                        |  3  display level 2 output plus additional
#     *                        |     values used by the global strategy
#     *                        |     (only if using KIN_LINESEARCH), and
#     *                        |     statistical information for the linear
#     *                        |     solver
#     *                        | [0]
#     *                        |
#     * KINSetInfoHandlerFn    | user-provided InfoHandler function.
#     *                        | [internal]
#     *                        |
#     * KINSetInfoFile         | pointer (type FILE) specifying where
#     *                        | informative (non-error) messages should
#     *                        | be sent if the default internal info
#     *                        | handler is used
#     *                        | [stdout]
#     *                        |
#     * KINSetUserData         | pointer to user-allocated memory that is
#     *                        | passed to the user-supplied subroutine
#     *                        | implementing the nonlinear system function
#     *                        | F(u)
#     *                        | [NULL]
#     *                        |
#     * KINSetNumMaxIters      | maximum number of nonlinear iterations
#     *                        | [MXITER_DEFAULT] (defined in kinsol_impl.h)
#     *                        |
#     * KINSetNoInitSetup      | flag controlling whether or not the
#     *                        | KINSol routine makes an initial call
#     *                        | to the linear solver setup routine (lsetup)
#     *                        | (possible values are TRUE and FALSE)
#     *                        | [FALSE]
#     *                        |
#     * KINSetNoResMon         | flag controlling whether or not the nonlinear
#     *                        | residual monitoring scheme is used to control
#     *                        | Jacobian updating (possible values are TRUE
#     *                        | and FALSE)
#     *                        | [FALSE if using direct linear solver]
#     *                        | [TRUE if using inexact linear solver]
#     *                        |
#     * KINSetMaxSetupCalls    | mbset, number of nonlinear iteraions, such 
#     *                        | that a call to the linear solver setup routine
#     *                        | (lsetup) is forced every mbset iterations.
#     *                        | If mbset=1, lsetup s called at every iteration.
#     *                        | [MSBSET_DEFAULT] (defined in kinsol_impl.h)
#     *                        |
#     * KINSetMaxSubSetupCalls | mbsetsub is the number of nonlinear iterations
#     *                        | between checks by the nonlinear residual
#     *                        | monitoring algorithm (specifies length of
#     *                        | subinterval)
#     *                        | NOTE: mbset should be a multiple of mbsetsub
#     *                        | [MSBSET_SUB_DEFAULT] (defined in kinsol_impl.h)
#     *                        |
#     * KINSetEtaForm          | flag indicating which method to use to
#     *                        | compute the value of the eta coefficient
#     *                        | used in the calculation of the linear
#     *                        | solver convergence tolerance:
#     *                        |
#     *                        |  eps = (eta+uround)*||fscale*func(u)||_L2
#     *                        |
#     *                        | the linear solver tests for convergence by
#     *                        | checking if the following inequality has
#     *                        | been satisfied:
#     *                        |
#     *                        |  ||fscale*(func(u)+J(u)*p)||_L2 <= eps
#     *                        |
#     *                        | where J(u) is the system Jacobian
#     *                        | evaluated at the current iterate, and p
#     *                        | denotes the Newton step
#     *                        |
#     *                        | choices for computing eta are as follows:
#     *                        |
#     *                        |  KIN_ETACHOICE1  (refer to KINForcingTerm)
#     *                        |
#     *                        |  eta = ABS(||F(u_k+1)||_L2-||F(u_k)+J(u_k)*p_k||_L2)
#     *                        |        ---------------------------------------------
#     *                        |                        ||F(u_k)||_L2
#     *                        | 
#     *                        |  KIN_ETACHOICE2  (refer to KINForcingTerm)
#     *                        |
#     *                        |                [ ||F(u_k+1)||_L2 ]^alpha
#     *                        |  eta = gamma * [ --------------- ]
#     *                        |                [  ||F(u_k)||_L2  ]
#     *                        |
#     *                        |  where gamma = [0,1] and alpha = (1,2]
#     *                        |
#     *                        |  KIN_ETACONSTANT  use a constant value for eta
#     *                        | [KIN_ETACHOICE1]
#     *                        |
#     * KINSetEtaConstValue    | constant value of eta - use with
#     *                        | KIN_ETACONSTANT option
#     *                        | [0.1]
#     *                        |
#     * KINSetEtaParams        | values of eta_gamma (egamma) and eta_alpha
#     *                        | (ealpha) coefficients - use with KIN_ETACHOICE2
#     *                        | option
#     *                        | [0.9 and 2.0]
#     *                        |
#     * KINSetResMonParams     | values of omega_min and omega_max scalars
#     *                        | used by nonlinear residual monitoring
#     *                        | algorithm (see KINStop)
#     *                        | [0.00001 and 0.9]
#     *                        |
#     * KINSetResMonConstValue | constant value used by residual monitoring
#     *                        | algorithm. If omega=0, then it is estimated
#     *                        | using omega_min and omega_max.
#     *                        | [0.0]
#     *                        |
#     * KINSetNoMinEps         | flag controlling whether or not the value
#     *                        | of eps is bounded below by 0.01*fnormtol
#     *                        | (see KINSetFuncNormTol)
#     *                        |
#     *                        |  FALSE  constrain value of eps by setting
#     *                        |         to the following:
#     *                        |
#     *                        |          eps = MAX{0.01*fnormtol, eps}
#     *                        |
#     *                        |  TRUE  do not constrain value of eps
#     *                        | [FALSE]
#     *                        |
#     * KINSetMaxNewtonStep    | maximum scaled length of Newton step
#     *                        | (reset to value of one if user-supplied
#     *                        | value is less than one)
#     *                        | [1000*||uscale*u_0||_L2]
#     *                        |
#     * KINSetMaxBetaFails     | maximum number of beta condition failures
#     *                        | in the line search algorithm.
#     *                        | [MXNBCF_DEFAULT] (defined in kinsol_impl.h)
#     *                        |
#     * KINSetRelErrFunc       | real scalar equal to realative error in
#     *                        | computing F(u) (used in difference-
#     *                        | quotient approximation of matrix-vector
#     *                        | product J(u)*v)
#     *                        | [(uround)^1/2]
#     *                        |
#     * KINSetFuncNormTol      | real scalar used as stopping tolerance on
#     *                        | ||fscale*ABS(func(u))||_L-infinity (see
#     *                        | KINStop and KINInitialStop)
#     *                        | [(uround)^1/3]
#     *                        |
#     * KINSetScaledStepTol    | real scalar used as stopping tolerance on
#     *                        | the maximum scaled step length:
#     *                        |
#     *                        |  ||    u_k+1 - u_k    ||
#     *                        |  || ----------------- ||_L-infinity
#     *                        |  || ABS(u_k+1)+uscale ||
#     *                        |
#     *                        | (see KINStop)
#     *                        | [(uround)^2/3]
#     *                        |
#     * KINSetConstraints      | pointer to an array (type N_Vector) of
#     *                        | constraints on the solution vector u
#     *                        | 
#     *                        | if constraints[i] =
#     *                        |
#     *                        |   0  u[i] not constrained
#     *                        |
#     *                        |  +1  u[i] constrained to be >= 0
#     *                        |  -1  u[i] constrained to be <= 0
#     *                        |
#     *                        |  +2  u[i] constrained to be > 0
#     *                        |  -2  u[i] constrained to be < 0
#     *                        |
#     *                        | if a NULL pointer is given, then no
#     *                        | constraints are applied to vector u
#     *                        | [NULL]
#     *                        |
#     * KINSetSysFunc          | set the user-provided routine which
#     *                        | defines the nonlinear problem to be 
#     *                        | solved
#     *                        | [none]
#     * -----------------------------------------------------------------
#     * The possible return values for the KINSet* subroutines are the
#     * following:
#     *
#     * KIN_SUCCESS : means the associated variable was successfully
#     *               set [0]
#     *
#     * KIN_MEM_NULL : means a NULL KINSOL memory block pointer was given
#     *                (must call the KINCreate and KINInit memory
#     *                allocation subroutines prior to calling KINSol) [-1]
#     *
#     * KIN_ILL_INPUT : means the supplied parameter was invalid (check
#     *                 error message) [-2]
#     * -----------------------------------------------------------------
#     * Note: If successful, these functions return KIN_SUCCESS. If an
#     * argument has an illegal value, then an error message is printed
#     * to the file specified by errfp and an error code is returned.
#     * -----------------------------------------------------------------
#     */
    
    int KINSetErrHandlerFn(void *kinmem, KINErrHandlerFn ehfun, void *eh_data)
    #int KINSetErrFile(void *kinmem, FILE *errfp)
    int KINSetInfoHandlerFn(void *kinmem, KINInfoHandlerFn ihfun, void *ih_data)
    #int KINSetInfoFile(void *kinmem, FILE *infofp)
    int KINSetUserData(void *kinmem, void *user_data)
    int KINSetPrintLevel(void *kinmemm, int printfl)
    int KINSetNumMaxIters(void *kinmem, long int mxiter)
    int KINSetNoInitSetup(void *kinmem, booleantype noInitSetup)
    int KINSetNoResMon(void *kinmem, booleantype noNNIResMon)
    int KINSetMaxSetupCalls(void *kinmem, long int msbset)
    int KINSetMaxSubSetupCalls(void *kinmem, long int msbsetsub)
    int KINSetEtaForm(void *kinmem, int etachoice)
    int KINSetEtaConstValue(void *kinmem, realtype eta)
    int KINSetEtaParams(void *kinmem, realtype egamma, realtype ealpha)
    int KINSetResMonParams(void *kinmem, realtype omegamin, realtype omegamax)
    int KINSetResMonConstValue(void *kinmem, realtype omegaconst)
    int KINSetNoMinEps(void *kinmem, booleantype noMinEps)
    int KINSetMaxNewtonStep(void *kinmem, realtype mxnewtstep)
    int KINSetMaxBetaFails(void *kinmem, long int mxnbcf)
    int KINSetRelErrFunc(void *kinmem, realtype relfunc)
    int KINSetFuncNormTol(void *kinmem, realtype fnormtol)
    int KINSetScaledStepTol(void *kinmem, realtype scsteptol)
    int KINSetConstraints(void *kinmem, N_Vector constraints)
    int KINSetSysFunc(void *kinmem, KINSysFn func)
    
#    /*
#     * -----------------------------------------------------------------
#     * Function : KINInit
#     * -----------------------------------------------------------------
#     * KINInit allocates additional memory for vector storage and
#     * sets a couple problem-specific KINSOL variables.
#     *
#     * Note: Additional vectors must be initialized by the user and
#     * passed to the KINSol routine.
#     *
#     *  kinmem  pointer to an internal memory block allocated during a
#     *          prior call to KINCreate
#     *
#     *  func  name of user-supplied subroutine implementing the
#     *        nonlinear function F(u)
#     *
#     *  tmpl  implementation-specific template vector (type N_Vector)
#     *        (created using either N_VNew_Serial or N_VNew_Parallel)
#     *
#     * KINInit return flags: KIN_SUCCESS, KIN_MEM_NULL, KIN_ILL_INPUT,
#     * and KIN_MEM_FAIL (see below). If an error occurs, then KINInit
#     * prints an error message.
#     *
#     * -----------------------------------------------------------------
#     * The possible return values for the KINInit subroutine are the
#     * following:
#     *
#     * KIN_SUCCESS : means the necessary system memory was successfully
#     *               allocated [0]
#     *
#     * KIN_MEM_NULL : means a NULL KINSOL memory block pointer was given
#     *                (must call the KINCreate routine before calling
#     *                KINInit) [-1]
#     *
#     * KIN_ILL_INPUT : means the name of a user-supplied subroutine
#     *                 implementing the nonlinear system function F(u)
#     *                 was not given [-2]
#     *
#     * KIN_MEM_FAIL : means an error occurred during memory allocation
#     *                (either insufficient system resources are available
#     *                or the vector kernel has not yet been initialized)
#     *                [-4]
#     * -----------------------------------------------------------------
#     */
    
    int KINInit(void *kinmem, KINSysFn func, N_Vector tmpl)
    
#    /*
#     * -----------------------------------------------------------------
#     * Function : KINSol
#     * -----------------------------------------------------------------
#     * KINSol (main KINSOL driver routine) manages the computational
#     * process of computing an approximate solution of the nonlinear
#     * system. If the initial guess (initial value assigned to vector u)
#     * doesn't violate any user-defined constraints, then the subroutine
#     * attempts to solve the system F(u) = 0 using a nonlinear Krylov
#     * subspace projection method. The Newton-Krylov iterations are
#     * stopped if either of the following conditions is satisfied:
#     *
#     *  ||F(u)||_L-infinity <= 0.01*fnormtol
#     *
#     *  ||u[i+1] - u[i]||_L-infinity <= scsteptol
#     *
#     * However, if the current iterate satisfies the second stopping
#     * criterion, it doesn't necessarily mean an approximate solution
#     * has been found since the algorithm may have stalled, or the
#     * user-specified step tolerance (scsteptol) may be too large.
#     *
#     *  kinmem  pointer to an internal memory block allocated during a
#     *          prior call to KINCreate
#     *
#     *  uu  vector set to initial guess by user before calling KINSol,
#     *      but which upon return contains an approximate solution of
#     *      the nonlinear system F(u) = 0
#     *
#     *  strategy  global strategy applied to Newton step if unsatisfactory
#     *            (KIN_NONE or KIN_LINESEARCH)
#     *
#     *  u_scale  vector containing diagonal elements of scaling matrix
#     *           for vector u chosen so that the components of
#     *           u_scale*u (as a matrix multiplication) all have
#     *           about the same magnitude when u is close to a root
#     *           of F(u)
#     *
#     *  f_scale  vector containing diagonal elements of scaling matrix
#     *           for F(u) chosen so that the components of
#     *           f_scale*F(u) (as a matrix multiplication) all have
#     *           roughly the same magnitude when u is not too near a
#     *           root of F(u)
#     *
#     * Note: The components of vectors u_scale and f_scale should be
#     * positive.
#     *
#     * If successful, KINSol returns a vector uu contains an approximate
#     * solution of the given nonlinear system. If an error occurs, then
#     * an error message is printed and an error code is returned.
#     *
#     * -----------------------------------------------------------------
#     * KINSol Return Values
#     * -----------------------------------------------------------------
#     *
#     * The possible return values for the KINSol subroutine are the
#     * following:
#     *
#     * KIN_SUCCESS : means ||fscale*ABS(func(u))||_L-infinity <= 0.01*fnormtol
#     *               and the current iterate uu is probably an approximate
#     *               solution of the nonlinear system F(u) = 0 [0]
#     *
#     * KIN_INITIAL_GUESS_OK : means the initial user-supplied guess
#     *                        already satisfies the stopping criterion
#     *                        given above [1]
#     *
#     * KIN_STEP_LT_STPTOL : means the following inequality has been
#     *                      satisfied (stopping tolerance on scaled
#     *                      step length):
#     *
#     *                    ||    u_k+1 - u_k    ||
#     *                    || ----------------- ||_L-infinity <= scsteptol
#     *                    || ABS(u_k+1)+uscale ||
#     *
#     *                      so the current iterate (denoted above by u_k+1)
#     *                      may be an approximate solution of the given
#     *                      nonlinear system, but it is also quite possible
#     *                      that the algorithm is "stalled" (making
#     *                      insufficient progress) near an invalid solution,
#     *                      or the real scalar scsteptol is too large [2]
#     *
#     * KIN_LINESEARCH_NONCONV : means the line search algorithm was unable
#     *                          to find an iterate sufficiently distinct
#     *                          from the current iterate
#     *
#     *                          failure to satisfy the sufficient decrease
#     *                          condition could mean the current iterate is
#     *                          "close" to an approximate solution of the given
#     *                          nonlinear system, the finite-difference
#     *                          approximation of the matrix-vector product
#     *                          J(u)*v is inaccurate, or the real scalar
#     *                          scsteptol is too large [-5]
#     *
#     * KIN_MAXITER_REACHED : means the maximum number of nonlinear iterations
#     *                       has been reached [-6]
#     *
#     * KIN_MXNEWT_5X_EXCEEDED : means five consecutive steps have been taken
#     *                          that satisfy the following inequality:
#     *
#     *                            ||uscale*p||_L2 > 0.99*mxnewtstep
#     *
#     *                          where p denotes the current step and
#     *                          mxnewtstep is a real scalar upper bound
#     *                          on the scaled step length
#     *
#     *                          such a failure may mean ||fscale*func(u)||_L2
#     *                          asymptotes from above to a finite value, or
#     *                          the real scalar mxnewtstep is too small [-7]
#     *
#     * KIN_LINESEARCH_BCFAIL : means the line search algorithm (implemented
#     *                         in KINLineSearch) was unable to satisfy the
#     *                         beta-condition for MXNBCF + 1 nonlinear
#     *                         iterations (not necessarily consecutive),
#     *                         which may indicate the algorithm is making
#     *                         poor progress [-8]
#     *
#     * KIN_LINSOLV_NO_RECOVERY : means the user-supplied routine psolve
#     *                           encountered a recoverable error, but
#     *                           the preconditioner is already current [-9]
#     *
#     * KIN_LINIT_FAIL : means the linear solver initialization routine (linit)
#     *                  encountered an error [-10]
#     *
#     * KIN_LSETUP_FAIL : means the user-supplied routine pset (used to compute
#     *                   the preconditioner) encountered an unrecoverable
#     *                   error [-11]
#     *
#     * KIN_LSOLVE_FAIL : means either the user-supplied routine psolve (used to
#     *                   to solve the preconditioned linear system) encountered
#     *                   an unrecoverable error, or the linear solver routine
#     *                   (lsolve) encountered an error condition [-12]
#     *
#     * KIN_MEM_NULL : means a NULL KINSOL memory block pointer was given
#     *                (must call the KINCreate and KINInit memory
#     *                allocation subroutines prior to calling KINSol) [-1]
#     *
#     * KIN_NO_MALLOC : means additional system memory has not yet been
#     *                 allocated for vector storage (forgot to call the
#     *                 KINInit routine) [-3]
#     *
#     * KIN_ILL_INPUT : means at least one input parameter was invalid
#     *                 (check error output) [-2]
#     * -----------------------------------------------------------------
#     */
    
    int KINSol(void *kinmem, N_Vector uu, int strategy,
    			   N_Vector u_scale, N_Vector f_scale)
    
#    /*
#     * -----------------------------------------------------------------
#     * Optional Output Extraction Functions (KINSOL)
#     * -----------------------------------------------------------------
#     * The following functions can be called to get optional outputs
#     * and statistical information related to the KINSOL solver:
#     *
#     *       Function Name       |      Returned Value
#     *                           |
#     * -----------------------------------------------------------------
#     *                           |
#     * KINGetWorkSpace           | returns both integer workspace size
#     *                           | (total number of long int-sized blocks
#     *                           | of memory allocated by KINSOL for
#     *                           | vector storage) and real workspace
#     *                           | size (total number of realtype-sized
#     *                           | blocks of memory allocated by KINSOL
#     *                           | for vector storage)
#     *                           |
#     * KINGetNumFuncEvals        | total number evaluations of the
#     *                           | nonlinear system function F(u)
#     *                           | (number of direct calls made to the
#     *                           | user-supplied subroutine by KINSOL
#     *                           | module member functions)
#     *                           |
#     * KINGetNumNonlinSolvIters  | total number of nonlinear iterations
#     *                           | performed
#     *                           |
#     * KINGetNumBetaCondFails    | total number of beta-condition
#     *                           | failures (see KINLineSearch)
#     *                           |
#     *                           | KINSOL halts if the number of such
#     *                           | failures exceeds the value of the
#     *                           | constant MXNBCF (defined in kinsol.c)
#     *                           |
#     * KINGetNumBacktrackOps     | total number of backtrack operations
#     *                           | (step length adjustments) performed
#     *                           | by the line search algorithm (see
#     *                           | KINLineSearch)
#     *                           |
#     * KINGetFuncNorm            | scaled norm of the nonlinear system
#     *                           | function F(u) evaluated at the
#     *                           | current iterate:
#     *                           |
#     *                           |  ||fscale*func(u)||_L2
#     *                           |
#     * KINGetStepLength          | scaled norm (or length) of the step
#     *                           | used during the previous iteration:
#     *                           |
#     *                           |  ||uscale*p||_L2
#     *                           |
#     * KINGetReturnFlagName      | returns the name of the constant
#     *                           | associated with a KINSOL return flag
#     *                           |
#     * -----------------------------------------------------------------
#     *
#     * The possible return values for the KINSet* subroutines are the
#     * following:
#     *
#     * KIN_SUCCESS : means the information was successfully retrieved [0]
#     * 
#     * KIN_MEM_NULL : means a NULL KINSOL memory block pointer was given
#     *                (must call the KINCreate and KINInit memory
#     *                allocation subroutines prior to calling KINSol) [-1]
#     * -----------------------------------------------------------------
#     */
    
    int KINGetWorkSpace(void *kinmem, long int *lenrw, long int *leniw)
    int KINGetNumNonlinSolvIters(void *kinmem, long int *nniters)
    int KINGetNumFuncEvals(void *kinmem, long int *nfevals)
    int KINGetNumBetaCondFails(void *kinmem, long int *nbcfails)
    int KINGetNumBacktrackOps(void *kinmem, long int *nbacktr)
    int KINGetFuncNorm(void *kinmem, realtype *fnorm)
    int KINGetStepLength(void *kinmem, realtype *steplength)
    char *KINGetReturnFlagName(long int flag)
    
#    /*
#     * -----------------------------------------------------------------
#     * Function : KINFree
#     * -----------------------------------------------------------------
#     * KINFree frees system memory resources reserved for the KINSOL
#     * solver module.
#     *
#     *  kinmem  pointer to an internal memory block allocated during
#     *          prior calls to KINCreate and KINInit
#     * -----------------------------------------------------------------
#     */
    
    void KINFree(void **kinmem)
    
cdef extern from "kinsol/kinsol_dense.h":
#    /*
#     * -----------------------------------------------------------------
#     * Function : KINDense
#     * -----------------------------------------------------------------
#     * A call to the KINDense function links the main solver with the
#     * KINDENSE linear solver. Its arguments are as follows:
#     *
#     * kinmem - pointer to an internal memory block allocated during a
#     *          prior call to KINCreate
#     *
#     * N      - problem size
#     *
#     * The return value of KINDense is one of:
#     *    KINDLS_SUCCESS   if successful
#     *    KINDLS_MEM_NULL  if the kinsol memory was NULL
#     *    KINDLS_MEM_FAIL  if there was a memory allocation failure
#     *    KINDLS_ILL_INPUT if a required vector operation is missing
#     * -----------------------------------------------------------------
#     */
#    
    int KINDense(void *kinmem, long int N)
   
cdef extern from "kinsol/kinsol_band.h":
#    /*
#     * -----------------------------------------------------------------
#     * Function : KINBand
#     * -----------------------------------------------------------------
#     * A call to the KINBand function links the main solver with the 
#     * KINBAND linear solver. Its arguments are as follows:
#     *
#     * kinmem - pointer to the integrator memory returned by KINCreate.
#     *
#     * N      - problem size
#     *
#     * mupper - upper bandwidth of the band Jacobian
#     *
#     * mlower - lower bandwidth of the band Jacobian
#     *
#     * The return value of KINBand is one of:
#     *    KINDLS_SUCCESS   if successful
#     *    KINDLS_MEM_NULL  if the kinsol memory was NULL
#     *    KINDLS_MEM_FAIL  if there was a memory allocation failure
#     *    KINDLS_ILL_INPUT if a required vector operation is missing
#     *                        or if a bandwidth has an illegal value.
#     * -----------------------------------------------------------------
#     */
    
    int KINBand(void *kinmem, long int N, long int mupper, long int mlower)

cdef extern from "kinsol/kinsol_direct.h":
#    /*
#     * =================================================================
#     *              K I N D I R E C T     C O N S T A N T S
#     * =================================================================
#     */
#    
#    /* 
#     * -----------------------------------------------------------------
#     * KINDLS return values 
#     * -----------------------------------------------------------------
#     */
#    
#    #define KINDLS_SUCCESS           0
#    #define KINDLS_MEM_NULL         -1
#    #define KINDLS_LMEM_NULL        -2
#    #define KINDLS_ILL_INPUT        -3
#    #define KINDLS_MEM_FAIL         -4
#    
#    /* Additional last_flag values */
#    
    enum: KINDLS_JACFUNC_UNRECVR
    enum: KINDLS_JACFUNC_RECVR
#    
#    /*
#     * =================================================================
#     *              F U N C T I O N   T Y P E S
#     * =================================================================
#     */
#    
#    /*
#     * -----------------------------------------------------------------
#     * Type: KINDlsDenseJacFn
#     * -----------------------------------------------------------------
#     *
#     * A dense Jacobian approximation function Jac must be of type 
#     * KINDlsDenseJacFn. Its parameters are:
#     *
#     * N        - problem size.
#     *
#     * u        - current iterate (unscaled) [input]
#     *
#     * fu       - vector (type N_Vector) containing result of nonlinear
#     *            system function evaluated at current iterate:
#     *            fu = F(u) [input]
#     *
#     * J        - dense matrix (of type DlsMat) that will be loaded
#     *            by a KINDlsDenseJacFn with an approximation to the
#     *            Jacobian matrix J = (dF_i/dy_j).
#     *
#     * user_data   - pointer to user data - the same as the user_data
#     *            parameter passed to KINSetFdata.
#     *
#     * tmp1, tmp2 - available scratch vectors (volatile storage)
#     *
#     * A KINDlsDenseJacFn should return 0 if successful, a positive 
#     * value if a recoverable error occurred, and a negative value if 
#     * an unrecoverable error occurred.
#     *
#     * -----------------------------------------------------------------
#     *
#     * NOTE: The following are two efficient ways to load a dense Jac:         
#     * (1) (with macros - no explicit data structure references)      
#     *     for (j=0; j < Neq; j++) {                                  
#     *       col_j = DENSE_COL(Jac,j);                                 
#     *       for (i=0; i < Neq; i++) {                                
#     *         generate J_ij = the (i,j)th Jacobian element           
#     *         col_j[i] = J_ij;                                       
#     *       }                                                        
#     *     }                                                          
#     * (2) (without macros - explicit data structure references)      
#     *     for (j=0; j < Neq; j++) {                                  
#     *       col_j = (Jac->data)[j];                                   
#     *       for (i=0; i < Neq; i++) {                                
#     *         generate J_ij = the (i,j)th Jacobian element           
#     *         col_j[i] = J_ij;                                       
#     *       }                                                        
#     *     }                                                          
#     * A third way, using the DENSE_ELEM(A,i,j) macro, is much less   
#     * efficient in general.  It is only appropriate for use in small 
#     * problems in which efficiency of access is NOT a major concern. 
#     *                                                                
#     * -----------------------------------------------------------------
#     */
#      
#      
    ctypedef int (*KINDlsDenseJacFn)(long int N,
    				N_Vector u, N_Vector fu, 
    				DlsMat J, void *user_data,
    				N_Vector tmp1, N_Vector tmp2)
#      
#    /*
#     * -----------------------------------------------------------------
#     * Type: KINDlsBandJacFn
#     * -----------------------------------------------------------------
#     *
#     * A band Jacobian approximation function Jac must have the
#     * prototype given below. Its parameters are:
#     *
#     * N is the problem size
#     *
#     * mupper is the upper half-bandwidth of the approximate banded
#     * Jacobian. This parameter is the same as the mupper parameter
#     * passed by the user to the linear solver initialization function.
#     *
#     * mlower is the lower half-bandwidth of the approximate banded
#     * Jacobian. This parameter is the same as the mlower parameter
#     * passed by the user to the linear solver initialization function.
#     *
#     * u        - current iterate (unscaled) [input]
#     *
#     * fu       - vector (type N_Vector) containing result of nonlinear
#     *            system function evaluated at current iterate:
#     *            fu = F(uu) [input]
#     *
#     * J        - band matrix (of type DlsMat) that will be loaded by a
#     *            KINDlsBandJacFn with an approximation to the Jacobian
#     *            matrix Jac = (dF_i/dy_j).
#     *
#     * user_data   - pointer to user data - the same as the user_data
#     *            parameter passed to KINSetFdata.
#     *
#     * tmp1, tmp2 - available scratch vectors (volatile storage)
#     *
#     * A KINDlsBandJacFn should return 0 if successful, a positive value
#     * if a recoverable error occurred, and a negative value if an 
#     * unrecoverable error occurred.
#     *
#     * -----------------------------------------------------------------
#     *
#     * NOTE. Three efficient ways to load J are:
#     *
#     * (1) (with macros - no explicit data structure references)
#     *    for (j=0; j < n; j++) {
#     *       col_j = BAND_COL(Jac,j);
#     *       for (i=j-mupper; i <= j+mlower; i++) {
#     *         generate J_ij = the (i,j)th Jacobian element
#     *         BAND_COL_ELEM(col_j,i,j) = J_ij;
#     *       }
#     *     }
#     *
#     * (2) (with BAND_COL macro, but without BAND_COL_ELEM macro)
#     *    for (j=0; j < n; j++) {
#     *       col_j = BAND_COL(Jac,j);
#     *       for (k=-mupper; k <= mlower; k++) {
#     *         generate J_ij = the (i,j)th Jacobian element, i=j+k
#     *         col_j[k] = J_ij;
#     *       }
#     *     }
#     *
#     * (3) (without macros - explicit data structure references)
#     *     offset = Jac->smu;
#     *     for (j=0; j < n; j++) {
#     *       col_j = ((Jac->data)[j])+offset;
#     *       for (k=-mupper; k <= mlower; k++) {
#     *         generate J_ij = the (i,j)th Jacobian element, i=j+k
#     *         col_j[k] = J_ij;
#     *       }
#     *     }
#     * Caution: Jac->smu is generally NOT the same as mupper.
#     *
#     * The BAND_ELEM(A,i,j) macro is appropriate for use in small
#     * problems in which efficiency of access is NOT a major concern.
#     *
#     * -----------------------------------------------------------------
#     */
#    
    ctypedef int (*KINDlsBandJacFn)(long int N, long int mupper, long int mlower,
    			       N_Vector u, N_Vector fu, 
    			       DlsMat J, void *user_data,
    			       N_Vector tmp1, N_Vector tmp2)
#    
#    /*
#     * =================================================================
#     *            E X P O R T E D    F U N C T I O N S 
#     * =================================================================
#     */
#    
#    /*
#     * -----------------------------------------------------------------
#     * Optional inputs to the KINDLS linear solver
#     * -----------------------------------------------------------------
#     *
#     * KINDlsSetDenseJacFn specifies the dense Jacobian approximation
#     * routine to be used for a direct dense linear solver.
#     *
#     * KINDlsSetBandJacFn specifies the band Jacobian approximation
#     * routine to be used for a direct band linear solver.
#     *
#     * By default, a difference quotient approximation, supplied with
#     * the solver is used.
#     *
#     * The return value is one of:
#     *    KINDLS_SUCCESS   if successful
#     *    KINDLS_MEM_NULL  if the KINSOL memory was NULL
#     *    KINDLS_LMEM_NULL if the linear solver memory was NULL
#     * -----------------------------------------------------------------
#     */
#    
    int KINDlsSetDenseJacFn(void *kinmem, KINDlsDenseJacFn jac)
    int KINDlsSetBandJacFn(void *kinmem, KINDlsBandJacFn jac)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Optional outputs from a KINDLS linear solver
#     * -----------------------------------------------------------------
#     *
#     * KINDlsGetWorkSpace    returns the real and integer workspace used
#     *                       by the KINDLS linear solver.
#     * KINDlsGetNumJacEvals  returns the number of calls made to the
#     *                       Jacobian evaluation routine.
#     * KINDlsGetNumFuncEvals returns the number of calls to the user's F
#     *                       routine due to finite difference Jacobian
#     *                       evaluation.
#     * KINDlsGetLastFlag     returns the last error flag set by any of
#     *                       the KINDLS interface functions.
#     * KINDlsGetReturnFlagName returns the name of the constant
#     *                       associated with a KINDLS return flag
#     *
#     * The return value of KINDlsGet* is one of:
#     *    KINDLS_SUCCESS   if successful
#     *    KINDLS_MEM_NULL  if the KINSOL memory was NULL
#     *    KINDLS_LMEM_NULL if the linear solver memory was NULL
#     * -----------------------------------------------------------------
#     */
#    
    int KINDlsGetWorkSpace(void *kinmem, long int *lenrwB, long int *leniwB)
    int KINDlsGetNumJacEvals(void *kinmem, long int *njevalsB)
    int KINDlsGetNumFuncEvals(void *kinmem, long int *nfevalsB)
    int KINDlsGetLastFlag(void *kinmem, long int *flag)
    char *KINDlsGetReturnFlagName(long int flag)

    
cdef extern from "kinsol/kinsol_spils.h":
#    /*
#     * -----------------------------------------------------------------
#     * KINSPILS return values
#     * -----------------------------------------------------------------
#     */
#    
    enum: KINSPILS_SUCCESS   # 0
#    
    enum: KINSPILS_MEM_NULL  #-1
    enum: KINSPILS_LMEM_NULL #-2
    enum: KINSPILS_ILL_INPUT #-3
    enum: KINSPILS_MEM_FAIL  #-4
    enum: KINSPILS_PMEM_NULL #-5
#    
#    /*
#     * -----------------------------------------------------------------
#     * KINSPILS solver constant
#     * -----------------------------------------------------------------
#     * KINSPILS_MAXL : default maximum dimension of Krylov subspace
#     * -----------------------------------------------------------------
#     */
#    
    enum: KINSPILS_MAXL #10
#    
#    /*
#     * -----------------------------------------------------------------
#     * Type : KINSpilsPrecSetupFn
#     * -----------------------------------------------------------------
#     * The user-supplied preconditioner setup subroutine should
#     * compute the right-preconditioner matrix P (stored in memory
#     * block referenced by P_data pointer) used to form the
#     * scaled preconditioned linear system:
#     *
#     *  (Df*J(uu)*(P^-1)*(Du^-1)) * (Du*P*x) = Df*(-F(uu))
#     *
#     * where Du and Df denote the diagonal scaling matrices whose
#     * diagonal elements are stored in the vectors uscale and
#     * fscale, repsectively.
#     *
#     * The preconditioner setup routine (referenced by iterative linear
#     * solver modules via pset (type KINSpilsPrecSetupFn)) will not be
#     * called prior to every call made to the psolve function, but will
#     * instead be called only as often as necessary to achieve convergence
#     * of the Newton iteration.
#     *
#     * Note: If the psolve routine requires no preparation, then a
#     * preconditioner setup function need not be given.
#     *
#     *  uu  current iterate (unscaled) [input]
#     *
#     *  uscale  vector (type N_Vector) containing diagonal elements
#     *          of scaling matrix for vector uu [input]
#     *
#     *  fval  vector (type N_Vector) containing result of nonliear
#     *        system function evaluated at current iterate:
#     *        fval = F(uu) [input]
#     *
#     *  fscale  vector (type N_Vector) containing diagonal elements
#     *          of scaling matrix for fval [input]
#     *
#     *  user_data  pointer to user-allocated data memory block
#     *
#     *  vtemp1/vtemp2  available scratch vectors (temporary storage)
#     *
#     * If successful, the function should return 0 (zero). If an error
#     * occurs, then the routine should return a non-zero integer value.
#     * -----------------------------------------------------------------
#     */
#    
    ctypedef int (*KINSpilsPrecSetupFn)(N_Vector uu, N_Vector uscale,
                                       N_Vector fval, N_Vector fscale,
                                       void *user_data, N_Vector vtemp1,
    				   N_Vector vtemp2)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Type : KINSpilsPrecSolveFn
#     * -----------------------------------------------------------------
#     * The user-supplied preconditioner solve subroutine (referenced
#     * by iterative linear solver modules via psolve (type
#     * KINSpilsPrecSolveFn)) should solve a (scaled) preconditioned
#     * linear system of the generic form P*z = r, where P denotes the
#     * right-preconditioner matrix computed by the pset routine.
#     *
#     *  uu  current iterate (unscaled) [input]
#     *
#     *  uscale  vector (type N_Vector) containing diagonal elements
#     *          of scaling matrix for vector uu [input]
#     *
#     *  fval  vector (type N_Vector) containing result of nonliear
#     *        system function evaluated at current iterate:
#     *        fval = F(uu) [input]
#     *
#     *  fscale  vector (type N_Vector) containing diagonal elements
#     *          of scaling matrix for fval [input]
#     *
#     *  vv  vector initially set to the right-hand side vector r, but
#     *      which upon return contains a solution of the linear system
#     *      P*z = r [input/output]
#     *
#     *  user_data  pointer to user-allocated data memory block
#     *
#     *  vtemp  available scratch vector (volatile storage)
#     *
#     * If successful, the function should return 0 (zero). If a
#     * recoverable error occurs, then the subroutine should return
#     * a positive integer value (in this case, KINSOL attempts to
#     * correct by calling the preconditioner setup function if the 
#     * preconditioner information is out of date). If an unrecoverable 
#     * error occurs, then the preconditioner solve function should return 
#     * a negative integer value.
#     * -----------------------------------------------------------------
#     */
#    
    ctypedef int (*KINSpilsPrecSolveFn)(N_Vector uu, N_Vector uscale, 
                                       N_Vector fval, N_Vector fscale, 
                                       N_Vector vv, void *user_data,
                                       N_Vector vtemp)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Type : KINSpilsJacTimesVecFn
#     * -----------------------------------------------------------------
#     * The (optional) user-supplied matrix-vector product subroutine
#     * (referenced internally via jtimes (type KINSpilsJacTimesVecFn))
#     * is used to compute Jv = J(uu)*v (system Jacobian applied to a
#     * given vector). If a user-defined routine is not given, then the
#     * private routine is used.
#     *
#     *  v  unscaled variant of vector to be multiplied by J(uu) [input]
#     *
#     *  Jv  vector containing result of matrix-vector product J(uu)*v
#     *      [output]
#     *
#     *  uu  current iterate (unscaled) [input]
#     *
#     *  new_uu  flag (reset by user) indicating if the iterate uu
#     *          has been updated in the interim - Jacobian needs
#     *          to be updated/reevaluated, if appropriate, unless
#     *          new_uu = FALSE [input/output]
#     *
#     *  user_data  pointer to user data, the same as the user_data
#     *          parameter passed to the KINSetUserData function.
#     *
#     * If successful, the function should return 0 (zero). If an error
#     * occurs, then the routine should return a non-zero integer value.
#     * -----------------------------------------------------------------
#     */
#    
    ctypedef int (*KINSpilsJacTimesVecFn)(N_Vector v, N_Vector Jv,
                                         N_Vector uu, booleantype *new_uu, 
                                         void *J_data)
    
    
    
    
    
#    /*
#     * -----------------------------------------------------------------
#     * Optional Input Specification Functions
#     * -----------------------------------------------------------------
#     * The following functions can be called to set optional inputs:
#     *
#     *       Function Name       |   Optional Input  [Default Value]
#     *                           |
#     * -----------------------------------------------------------------
#     *                           |
#     * KINSpilsSetMaxRestarts    | maximum number of times the SPGMR
#     *                           | (scaled preconditioned GMRES) linear
#     *                           | solver can be restarted
#     *                           | [0]
#     *                           |
#     * KINSpilsSetPreconditioner | used to set the following:
#     *                           |   (a) name of user-supplied routine
#     *                           |       used to compute a preconditioner
#     *                           |       matrix for the given linear
#     *                           |       system (pset)
#     *                           |       [NULL]
#     *                           |   (b) name of user-supplied routine
#     *                           |       used to apply preconditioner to
#     *                           |       linear system (psolve)
#     *                           |       [NULL]
#     *                           |
#     * KINSpilsSetJacTimesVecFn  | used to set the following the name
#     *                           | of user-supplied subroutine used to 
#     *                           | compute the matrix-vector product J(u)*v,
#     *                           | where J denotes the system Jacobian.
#     *                           | [KINSpilsDQJtimes]
#     * -----------------------------------------------------------------
#     */
#    
    int KINSpilsSetMaxRestarts(void *kinmem, int maxrs)
    int KINSpilsSetPreconditioner(void *kinmem,
    					      KINSpilsPrecSetupFn pset,
    					      KINSpilsPrecSolveFn psolve)
    int KINSpilsSetJacTimesVecFn(void *kinmem,
                                                 KINSpilsJacTimesVecFn jtv)
#    
#    /*
#     * -----------------------------------------------------------------
#     * KINSpilsSet* Return Values
#     * -----------------------------------------------------------------
#     * The possible return values for the KINSpilsSet* subroutines
#     * are the following:
#     *
#     * KINSPILS_SUCCESS : means the associated parameter was successfully
#     *                    set [0]
#     *
#     * KINSPILS_ILL_INPUT : means the supplied parameter was invalid
#     *                      (check error message) [-3]
#     *
#     * KINSPILS_MEM_NULL : means a NULL KINSOL memory block pointer
#     *                     was given [-1]
#     *
#     * KINSPILS_LMEM_NULL : means system memory has not yet been
#     *                      allocated for the linear solver 
#     *                      (lmem == NULL) [-2]
#     * -----------------------------------------------------------------
#     */
#    
#    /*
#     * -----------------------------------------------------------------
#     * Optional Output Extraction Functions
#     * -----------------------------------------------------------------
#     * The following functions can be called to get optional outputs
#     * and statistical information related to the KINSPILS linear
#     * solvers:
#     *
#     *        Function Name       |      Returned Value
#     *                            |
#     * -----------------------------------------------------------------
#     *                            |
#     * KINSpilsGetWorkSpace       | returns both integer workspace size
#     *                            | (total number of long int-sized blocks
#     *                            | of memory allocated  for
#     *                            | vector storage), and real workspace
#     *                            | size (total number of realtype-sized
#     *                            | blocks of memory allocated
#     *                            | for vector storage)
#     *                            |
#     * KINSpilsGetNumPrecEvals    | total number of preconditioner
#     *                            | evaluations (number of calls made
#     *                            | to the user-defined pset routine)
#     *                            |
#     * KINSpilsGetNumPrecSolves   | total number of times preconditioner
#     *                            | was applied to linear system (number
#     *                            | of calls made to the user-supplied
#     *                            | psolve function)
#     *                            |
#     * KINSpilsGetNumLinIters     | total number of linear iterations
#     *                            | performed
#     *                            |
#     * KINSpilsGetNumConvFails    | total number of linear convergence
#     *                            | failures
#     *                            |
#     * KINSpilsGetNumJtimesEvals  | total number of times the matrix-
#     *                            | vector product J(u)*v was computed
#     *                            | (number of calls made to the jtimes
#     *                            | subroutine)
#     *                            |
#     * KINSpilsGetNumFuncEvals    | total number of evaluations of the
#     *                            | system function F(u) (number of
#     *                            | calls made to the user-supplied
#     *                            | func routine by the linear solver
#     *                            | module member subroutines)
#     *                            |
#     * KINSpilsGetLastFlag        | returns the last flag returned by
#     *                            | the linear solver
#     *                            |
#     * KINSpilsGetReturnFlagName  | returns the name of the constant
#     *                            | associated with a KINSPILS return flag
#     * -----------------------------------------------------------------
#     */
#    
    int KINSpilsGetWorkSpace(void *kinmem, long int *lenrwSG, long int *leniwSG)
    int KINSpilsGetNumPrecEvals(void *kinmem, long int *npevals)
    int KINSpilsGetNumPrecSolves(void *kinmem, long int *npsolves)
    int KINSpilsGetNumLinIters(void *kinmem, long int *nliters)
    int KINSpilsGetNumConvFails(void *kinmem, long int *nlcfails)
    int KINSpilsGetNumJtimesEvals(void *kinmem, long int *njvevals)
    int KINSpilsGetNumFuncEvals(void *kinmem, long int *nfevalsS)
    int KINSpilsGetLastFlag(void *kinmem, long int *flag)
    char *KINSpilsGetReturnFlagName(long int flag)
    
cdef extern from "kinsol/kinsol_spgmr.h":
#    /*
#     * -----------------------------------------------------------------
#     * Function : KINSpgmr
#     * -----------------------------------------------------------------
#     * KINSpgmr links the main KINSOL solver module with the SPGMR
#     * linear solver module. The routine establishes the inter-module
#     * interface by setting the generic KINSOL pointers linit, lsetup,
#     * lsolve, and lfree to KINSpgmrInit, KINSpgmrSetup, KINSpgmrSolve,
#     * and KINSpgmrFree, respectively.
#     *
#     *  kinmem  pointer to an internal memory block allocated during a
#     *          prior call to KINCreate
#     *
#     *  maxl  maximum allowable dimension of Krylov subspace (passing
#     *        a value of 0 (zero) will cause the default value
#     *        KINSPILS_MAXL (predefined constant) to be used)
#     *
#     * -----------------------------------------------------------------
#     * KINSpgmr Return Values
#     * -----------------------------------------------------------------
#     *
#     * The possible return values for the KINSpgmr subroutine are the
#     * following:
#     *
#     * KINSPILS_SUCCESS : means the KINSPGMR linear solver module
#     *                    (implementation of the GMRES method) was
#     *                    successfully initialized - allocated system
#     *                    memory and set shared variables to default
#     *                    values [0]
#     *
#     * KINSPILS_MEM_NULL : means a NULL KINSOL memory block pointer was
#     *                     given (must call the KINCreate and KINMalloc
#     *                     memory allocation subroutines prior to
#     *                     calling KINSpgmr) [-1]
#     *
#     * KINSPILS_MEM_FAIL : means either insufficient system resources
#     *                     were available to allocate memory for the main
#     *                     KINSPGMR data structure (type KINSpgmrMemRec),
#     *                     or the SpgmrMalloc subroutine failed (unable
#     *                     to allocate enough system memory for vector
#     *                     storage and/or the main SPGMR data structure
#     *                     (type SpgmrMemRec)) [-4]
#     *
#     * KINSPILS_ILL_INPUT : means a supplied parameter was invalid
#     *                      (check error message) [-3]
#     *
#     * The above constants are defined in kinsol_spils.h
#     * -----------------------------------------------------------------
#     */
    
    int KINSpgmr(void *kinmem, int maxl)
    
cdef extern from "kinsol/kinsol_spbcgs.h":
#    /*
#     * -----------------------------------------------------------------
#     * Function : KINSpbcg
#     * -----------------------------------------------------------------
#     * KINSpbcg links the main KINSOL solver module with the SPBCG
#     * linear solver module. The routine establishes the inter-module
#     * interface by setting the generic KINSOL pointers linit, lsetup,
#     * lsolve, and lfree to KINSpbcgInit, KINSpbcgSetup, KINSpbcgSolve,
#     * and KINSpbcgFree, respectively.
#     *
#     *  kinmem  pointer to an internal memory block allocated during a
#     *          prior call to KINCreate
#     *
#     *  maxl  maximum allowable dimension of Krylov subspace (passing
#     *        a value of 0 (zero) will cause the default value
#     *        KINSPILS_MAXL (predefined constant) to be used)
#     *
#     * If successful, KINSpbcg returns KINSPILS_SUCCESS. If an error
#     * occurs, then KINSpbcg returns an error code (negative integer
#     * value).
#     *
#     * -----------------------------------------------------------------
#     * KINSpbcg Return Values
#     * -----------------------------------------------------------------
#     * The possible return values for the KINSpbcg subroutine are the
#     * following:
#     *
#     * KINSPILS_SUCCESS : means the KINSPBCG linear solver module
#     *                    (implementation of the Bi-CGSTAB method) was
#     *                    successfully initialized - allocated system
#     *                    memory and set shared variables to default
#     *                    values [0]
#     *
#     * KINSPILS_MEM_NULL : means a NULL KINSOL memory block pointer
#     *                     was given (must call the KINCreate and
#     *                     KINMalloc memory allocation subroutines
#     *                     prior to calling KINSpbcg) [-1]
#     *
#     * KINSPILS_MEM_FAIL : means either insufficient system resources
#     *                     were available to allocate memory for the
#     *                     main KINSPBCG data structure (type
#     *                     KINSpbcgMemRec), or the SpbcgMalloc subroutine
#     *                     failed (unable to allocate enough system
#     *                     memory for vector storate and/or the main
#     *                     SPBCG data structure (type SpbcgMemRec)) [-4]
#     *
#     * KINSPILS_ILL_INPUT : means either a supplied parameter was invalid,
#     *                      or the NVECTOR implementation is NOT
#     *                      compatible [-3]
#     *
#     * The above constants are defined in kinsol_spils.h
#     * -----------------------------------------------------------------
#     */
    
    int KINSpbcg(void *kinmem, int maxl)
    
cdef extern from "kinsol/kinsol_sptfqmr.h":
#    
#    /*
#     * -----------------------------------------------------------------
#     * Function : KINSptfqmr
#     * -----------------------------------------------------------------
#     * KINSptfqmr links the main KINSOL solver module with the SPTFQMR
#     * linear solver module. The routine establishes the inter-module
#     * interface by setting the generic KINSOL pointers linit, lsetup,
#     * lsolve, and lfree to KINSptfqmrInit, KINSptfqmrSetup, KINSptfqmrSolve,
#     * and KINSptfqmrFree, respectively.
#     *
#     *  kinmem  pointer to an internal memory block allocated during a
#     *          prior call to KINCreate
#     *
#     *  maxl  maximum allowable dimension of Krylov subspace (passing
#     *        a value of 0 (zero) will cause the default value
#     *        KINSPTFQMR_MAXL (predefined constant) to be used)
#     *
#     * If successful, KINSptfqmr returns KINSPTFQMR_SUCCESS. If an error
#     * occurs, then KINSptfqmr returns an error code (negative integer
#     * value).
#     *
#     * -----------------------------------------------------------------
#     * KINSptfqmr Return Values
#     * -----------------------------------------------------------------
#     * The possible return values for the KINSptfqmr subroutine are the
#     * following:
#     *
#     * KINSPTFQMR_SUCCESS : means the KINSPTFQMR linear solver module
#     *                      (implementation of the TFQMR method) was
#     *                      successfully initialized - allocated system
#     *                      memory and set shared variables to default
#     *                      values [0]
#     *
#     * KINSPTFQMR_MEM_NULL : means a NULL KINSOL memory block pointer
#     *                       was given (must call the KINCreate and
#     *                       KINMalloc memory allocation subroutines
#     *                       prior to calling KINSptfqmr) [-1]
#     *
#     * KINSPTFQMR_MEM_FAIL : means either insufficient system resources
#     *                       were available to allocate memory for the
#     *                       main KINSPTFQMR data structure (type
#     *                       KINSptfqmrMemRec), or the SptfqmrMalloc
#     *                       subroutine failed (unable to allocate enough
#     *                       system memory for vector storate and/or the
#     *                       main SPTFQMR data structure
#     *                       (type SptfqmrMemRec)) [-4]
#     *
#     * KINSPTFQMR_ILL_INPUT : means either a supplied parameter was invalid,
#     *                        or the NVECTOR implementation is NOT
#     *                        compatible [-3]
#     *
#     * The above constants are defined in kinsol_spils.h
#     * -----------------------------------------------------------------
#     */

    int KINSptfqmr(void *kinmem, int maxl)
    
    
    
cdef extern from *: # actually defined in private header kinsol_impl.h
    cdef struct KINMemRec:

        realtype kin_uround        #/* machine epsilon (or unit roundoff error) 
				 #(defined in sundials_types.h)                */

        #/* problem specification data */

        #KINSysFn kin_func;           /* nonlinear system function implementation     */
        #void *kin_user_data;         /* work space available to func routine         */
        realtype kin_fnormtol #      /* stopping tolerance on L2-norm of function
                    		#		  value                                        */
        realtype kin_scsteptol  #    /* scaled step length tolerance                 */
        int kin_globalstrategy  #    /* choices are KIN_NONE and KIN_LINESEARCH      */
        int kin_printfl         #    /* level of verbosity of output                 */
        long int kin_mxiter     #    /* maximum number of nonlinear iterations       */
    