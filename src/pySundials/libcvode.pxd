from libsundials cimport *

cdef extern from "cvode/cvode.h":
    
#    /*
#     * =================================================================
#     *              C V O D E     C O N S T A N T S
#     * =================================================================
#     */
#    
#    /*
#     * -----------------------------------------------------------------
#     * Enumerations for inputs to CVodeCreate and CVode.
#     * -----------------------------------------------------------------
#     * Symbolic constants for the lmm and iter parameters to CVodeCreate
#     * and the input parameter itask to CVode, are given below.
#     *
#     * lmm:   The user of the CVODE package specifies whether to use the
#     *        CV_ADAMS (Adams-Moulton) or CV_BDF (Backward Differentiation
#     *        Formula) linear multistep method. The BDF method is
#     *        recommended for stiff problems, and the CV_ADAMS method is
#     *        recommended for nonstiff problems.
#     *
#     * iter:  At each internal time step, a nonlinear equation must
#     *        be solved. The user can specify either CV_FUNCTIONAL
#     *        iteration, which does not require linear algebra, or a
#     *        CV_NEWTON iteration, which requires the solution of linear
#     *        systems. In the CV_NEWTON case, the user also specifies a
#     *        CVODE linear solver. CV_NEWTON is recommended in case of
#     *        stiff problems.
#     *
#     * itask: The itask input parameter to CVode indicates the job
#     *        of the solver for the next user step. The CV_NORMAL
#     *        itask is to have the solver take internal steps until
#     *        it has reached or just passed the user specified tout
#     *        parameter. The solver then interpolates in order to
#     *        return an approximate value of y(tout). The CV_ONE_STEP
#     *        option tells the solver to just take one internal step
#     *        and return the solution at the point reached by that step.
#     * -----------------------------------------------------------------
#     */
    
    #/* lmm */
    enum: CV_ADAMS
    enum: CV_BDF
    
    #/* iter */
    enum: CV_FUNCTIONAL
    enum: CV_NEWTON
    
    #/* itask */
    enum: CV_NORMAL
    enum: CV_ONE_STEP

    

#    /*
#     * ----------------------------------------
#     * CVODE return flags
#     * ----------------------------------------
#     */

    enum: CV_SUCCESS
    enum: CV_TSTOP_RETURN
    enum: CV_ROOT_RETURN
    
    enum: CV_WARNING
    
    enum: CV_TOO_MUCH_WORK
    enum: CV_TOO_MUCH_ACC
    enum: CV_ERR_FAILURE
    enum: CV_CONV_FAILURE
    
    enum: CV_LINIT_FAIL
    enum: CV_LSETUP_FAIL
    enum: CV_LSOLVE_FAIL
    enum: CV_RHSFUNC_FAIL
    enum: CV_FIRST_RHSFUNC_ERR
    enum: CV_REPTD_RHSFUNC_ERR    
    enum: CV_UNREC_RHSFUNC_ERR
    enum: CV_RTFUNC_FAIL
    
    enum: CV_MEM_FAIL
    enum: CV_MEM_NULL
    enum: CV_ILL_INPUT
    enum: CV_NO_MALLOC
    enum: CV_BAD_K
    enum: CV_BAD_T
    enum: CV_BAD_DKY
    enum: CV_TOO_CLOSE


#    /*
#     * =================================================================
#     *              F U N C T I O N   T Y P E S
#     * =================================================================
#     */
#    
#    /*
#     * -----------------------------------------------------------------
#     * Type : CVRhsFn
#     * -----------------------------------------------------------------
#     * The f function which defines the right hand side of the ODE
#     * system y' = f(t,y) must have type CVRhsFn.
#     * f takes as input the independent variable value t, and the
#     * dependent variable vector y.  It stores the result of f(t,y)
#     * in the vector ydot.  The y and ydot arguments are of type
#     * N_Vector.
#     * (Allocation of memory for ydot is handled within CVODE)
#     * The user_data parameter is the same as the user_data
#     * parameter set by the user through the CVodeSetUserData routine.
#     * This user-supplied pointer is passed to the user's f function
#     * every time it is called.
#     *
#     * A CVRhsFn should return 0 if successful, a negative value if
#     * an unrecoverable error occured, and a positive value if a 
#     * recoverable error (e.g. invalid y values) occured. 
#     * If an unrecoverable occured, the integration is halted. 
#     * If a recoverable error occured, then (in most cases) CVODE
#     * will try to correct and retry.
#     * -----------------------------------------------------------------
#     */

    ctypedef int (*CVRhsFn)(realtype t, N_Vector y,
    		       N_Vector ydot, void *user_data)

#    /*
#     * -----------------------------------------------------------------
#     * Type : CVRootFn
#     * -----------------------------------------------------------------
#     * A function g, which defines a set of functions g_i(t,y) whose
#     * roots are sought during the integration, must have type CVRootFn.
#     * The function g takes as input the independent variable value
#     * t, and the dependent variable vector y.  It stores the nrtfn
#     * values g_i(t,y) in the realtype array gout.
#     * (Allocation of memory for gout is handled within CVODE.)
#     * The user_data parameter is the same as that passed by the user
#     * to the CVodeSetUserData routine.  This user-supplied pointer is
#     * passed to the user's g function every time it is called.
#     *
#     * A CVRootFn should return 0 if successful or a non-zero value
#     * if an error occured (in which case the integration will be halted).
#     * -----------------------------------------------------------------
#     */
#    
    ctypedef int (*CVRootFn)(realtype t, N_Vector y, realtype *gout, void *user_data)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Type : CVEwtFn
#     * -----------------------------------------------------------------
#     * A function e, which sets the error weight vector ewt, must have
#     * type CVEwtFn.
#     * The function e takes as input the current dependent variable y.
#     * It must set the vector of error weights used in the WRMS norm:
#     * 
#     *   ||y||_WRMS = sqrt [ 1/N * sum ( ewt_i * y_i)^2 ]
#     *
#     * Typically, the vector ewt has components:
#     * 
#     *   ewt_i = 1 / (reltol * |y_i| + abstol_i)
#     *
#     * The user_data parameter is the same as that passed by the user
#     * to the CVodeSetUserData routine.  This user-supplied pointer is
#     * passed to the user's e function every time it is called.
#     * A CVEwtFn e must return 0 if the error weight vector has been
#     * successfuly set and a non-zero value otherwise.
#     * -----------------------------------------------------------------
#     */
#    
    ctypedef int (*CVEwtFn)(N_Vector y, N_Vector ewt, void *user_data)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Type : CVErrHandlerFn
#     * -----------------------------------------------------------------
#     * A function eh, which handles error messages, must have type
#     * CVErrHandlerFn.
#     * The function eh takes as input the error code, the name of the
#     * module reporting the error, the error message, and a pointer to
#     * user data, the same as that passed to CVodeSetUserData.
#     * 
#     * All error codes are negative, except CV_WARNING which indicates 
#     * a warning (the solver continues).
#     *
#     * A CVErrHandlerFn has no return value.
#     * -----------------------------------------------------------------
#     */
#    
    ctypedef void (*CVErrHandlerFn)(int error_code,
    			       char *module, char *function,
    			       char *msg, void *user_data)
#    
#    /*
#     * =================================================================
#     *          U S E R - C A L L A B L E   R O U T I N E S
#     * =================================================================
#     */
#    
#    /*
#     * -----------------------------------------------------------------
#     * Function : CVodeCreate
#     * -----------------------------------------------------------------
#     * CVodeCreate creates an internal memory block for a problem to
#     * be solved by CVODE.
#     *
#     * lmm   is the type of linear multistep method to be used.
#     *       The legal values are CV_ADAMS and CV_BDF (see previous
#     *       description).
#     *
#     * iter  is the type of iteration used to solve the nonlinear
#     *       system that arises during each internal time step.
#     *       The legal values are CV_FUNCTIONAL and CV_NEWTON.
#     *
#     * If successful, CVodeCreate returns a pointer to initialized
#     * problem memory. This pointer should be passed to CVodeInit.
#     * If an initialization error occurs, CVodeCreate prints an error
#     * message to standard err and returns NULL.
#     * -----------------------------------------------------------------
#     */
#    
    void *CVodeCreate(int lmm, int iter)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Integrator optional input specification functions
#     * -----------------------------------------------------------------
#     * The following functions can be called to set optional inputs
#     * to values other than the defaults given below:
#     *
#     * Function                |  Optional input / [ default value ]
#     * -----------------------------------------------------------------
#     *                         |
#     * CVodeSetErrHandlerFn    | user-provided ErrHandler function.
#     *                         | [internal]
#     *                         |
#     * CVodeSetErrFile         | the file pointer for an error file
#     *                         | where all CVODE warning and error
#     *                         | messages will be written if the default
#     *                         | internal error handling function is used. 
#     *                         | This parameter can be stdout (standard 
#     *                         | output), stderr (standard error), or a 
#     *                         | file pointer (corresponding to a user 
#     *                         | error file opened for writing) returned 
#     *                         | by fopen.
#     *                         | If not called, then all messages will
#     *                         | be written to the standard error stream.
#     *                         | [stderr]
#     *                         |
#     * CVodeSetUserData        | a pointer to user data that will be
#     *                         | passed to the user's f function every
#     *                         | time f is called.
#     *                         | [NULL]
#     *                         |
#     * CVodeSetMaxOrd          | maximum lmm order to be used by the
#     *                         | solver.
#     *                         | [12 for Adams , 5 for BDF]
#     *                         |
#     * CVodeSetMaxNumSteps     | maximum number of internal steps to be
#     *                         | taken by the solver in its attempt to
#     *                         | reach tout.
#     *                         | [500]
#     *                         |
#     * CVodeSetMaxHnilWarns    | maximum number of warning messages
#     *                         | issued by the solver that t+h==t on the
#     *                         | next internal step. A value of -1 means
#     *                         | no such messages are issued.
#     *                         | [10]
#     *                         |
#     * CVodeSetStabLimDet      | flag to turn on/off stability limit
#     *                         | detection (TRUE = on, FALSE = off).
#     *                         | When BDF is used and order is 3 or
#     *                         | greater, CVsldet is called to detect
#     *                         | stability limit.  If limit is detected,
#     *                         | the order is reduced.
#     *                         | [FALSE]
#     *                         |
#     * CVodeSetInitStep        | initial step size.
#     *                         | [estimated by CVODE]
#     *                         |
#     * CVodeSetMinStep         | minimum absolute value of step size
#     *                         | allowed.
#     *                         | [0.0]
#     *                         |
#     * CVodeSetMaxStep         | maximum absolute value of step size
#     *                         | allowed.
#     *                         | [infinity]
#     *                         |
#     * CVodeSetStopTime        | the independent variable value past
#     *                         | which the solution is not to proceed.
#     *                         | [infinity]
#     *                         |
#     * CVodeSetMaxErrTestFails | Maximum number of error test failures
#     *                         | in attempting one step.
#     *                         | [7]
#     *                         |
#     * CVodeSetMaxNonlinIters  | Maximum number of nonlinear solver
#     *                         | iterations at one solution.
#     *                         | [3]
#     *                         |
#     * CVodeSetMaxConvFails    | Maximum number of convergence failures
#     *                         | allowed in attempting one step.
#     *                         | [10]
#     *                         |
#     * CVodeSetNonlinConvCoef  | Coefficient in the nonlinear
#     *                         | convergence test.
#     *                         | [0.1]
#     *                         |
#     * -----------------------------------------------------------------
#     *                         |
#     * CVodeSetIterType        | Changes the current nonlinear iteration
#     *                         | type.
#     *                         | [set by CVodecreate]
#     *                         |
#     * -----------------------------------------------------------------
#     *                            |
#     * CVodeSetRootDirection      | Specifies the direction of zero
#     *                            | crossings to be monitored
#     *                            | [both directions]
#     *                            |
#     * CVodeSetNoInactiveRootWarn | disable warning about possible
#     *                            | g==0 at beginning of integration
#     *                            | 
#     * -----------------------------------------------------------------
#    
#     * -----------------------------------------------------------------
#     * Return flag:
#     *   CV_SUCCESS   if successful
#     *   CV_MEM_NULL  if the cvode memory is NULL
#     *   CV_ILL_INPUT if an argument has an illegal value
#     * -----------------------------------------------------------------
#     */
#    
    int CVodeSetErrHandlerFn(void *cvode_mem, CVErrHandlerFn ehfun, void *eh_data)
    #int CVodeSetErrFile(void *cvode_mem, FILE *errfp)
    int CVodeSetUserData(void *cvode_mem, void *user_data)
    int CVodeSetMaxOrd(void *cvode_mem, int maxord)
    int CVodeSetMaxNumSteps(void *cvode_mem, long int mxsteps)
    int CVodeSetMaxHnilWarns(void *cvode_mem, int mxhnil)
    int CVodeSetStabLimDet(void *cvode_mem, booleantype stldet)
    int CVodeSetInitStep(void *cvode_mem, realtype hin)
    int CVodeSetMinStep(void *cvode_mem, realtype hmin)
    int CVodeSetMaxStep(void *cvode_mem, realtype hmax)
    int CVodeSetStopTime(void *cvode_mem, realtype tstop)
    int CVodeSetMaxErrTestFails(void *cvode_mem, int maxnef)
    int CVodeSetMaxNonlinIters(void *cvode_mem, int maxcor)
    int CVodeSetMaxConvFails(void *cvode_mem, int maxncf)
    int CVodeSetNonlinConvCoef(void *cvode_mem, realtype nlscoef)
    
    int CVodeSetIterType(void *cvode_mem, int iter)
    
    int CVodeSetRootDirection(void *cvode_mem, int *rootdir)
    int CVodeSetNoInactiveRootWarn(void *cvode_mem)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Function : CVodeInit
#     * -----------------------------------------------------------------
#     * CVodeInit allocates and initializes memory for a problem to
#     * to be solved by CVODE.
#     *
#     * cvode_mem is pointer to CVODE memory returned by CVodeCreate.
#     *
#     * f       is the name of the C function defining the right-hand
#     *         side function in y' = f(t,y).
#     *
#     * t0      is the initial value of t.
#     *
#     * y0      is the initial condition vector y(t0).
#     *
#     * Return flag:
#     *  CV_SUCCESS if successful
#     *  CV_MEM_NULL if the cvode memory was NULL
#     *  CV_MEM_FAIL if a memory allocation failed
#     *  CV_ILL_INPUT f an argument has an illegal value.
#     * -----------------------------------------------------------------
#     */
#    
    int CVodeInit(void *cvode_mem, CVRhsFn f, realtype t0, N_Vector y0)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Function : CVodeReInit
#     * -----------------------------------------------------------------
#     * CVodeReInit re-initializes CVode for the solution of a problem,
#     * where a prior call to CVodeInit has been made with the same
#     * problem size N. CVodeReInit performs the same input checking
#     * and initializations that CVodeInit does.
#     * But it does no memory allocation, assuming that the existing
#     * internal memory is sufficient for the new problem.
#     *
#     * The use of CVodeReInit requires that the maximum method order,
#     * maxord, is no larger for the new problem than for the problem
#     * specified in the last call to CVodeInit.  This condition is
#     * automatically fulfilled if the multistep method parameter lmm
#     * is unchanged (or changed from CV_ADAMS to CV_BDF) and the default
#     * value for maxord is specified.
#     *
#     * All of the arguments to CVodeReInit have names and meanings
#     * identical to those of CVodeInit.
#     *
#     * The return value of CVodeReInit is equal to CV_SUCCESS = 0 if
#     * there were no errors; otherwise it is a negative int equal to:
#     *   CV_MEM_NULL      indicating cvode_mem was NULL (i.e.,
#     *                    CVodeCreate has not been called).
#     *   CV_NO_MALLOC     indicating that cvode_mem has not been
#     *                    allocated (i.e., CVodeInit has not been
#     *                    called).
#     *   CV_ILL_INPUT     indicating an input argument was illegal
#     *                    (including an attempt to increase maxord).
#     * In case of an error return, an error message is also printed.
#     * -----------------------------------------------------------------
#     */
#    
    int CVodeReInit(void *cvode_mem, realtype t0, N_Vector y0)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Functions : CVodeSStolerances
#     *             CVodeSVtolerances
#     *             CVodeWFtolerances
#     * -----------------------------------------------------------------
#     *
#     * These functions specify the integration tolerances. One of them
#     * MUST be called before the first call to CVode.
#     *
#     * CVodeSStolerances specifies scalar relative and absolute tolerances.
#     * CVodeSVtolerances specifies scalar relative tolerance and a vector
#     *   absolute tolerance (a potentially different absolute tolerance 
#     *   for each vector component).
#     * CVodeWFtolerances specifies a user-provides function (of type CVEwtFn)
#     *   which will be called to set the error weight vector.
#     *
#     * The tolerances reltol and abstol define a vector of error weights,
#     * ewt, with components
#     *   ewt[i] = 1/(reltol*abs(y[i]) + abstol)      (in the SS case), or
#     *   ewt[i] = 1/(reltol*abs(y[i]) + abstol[i])   (in the SV case).
#     * This vector is used in all error and convergence tests, which
#     * use a weighted RMS norm on all error-like vectors v:
#     *    WRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v[i]*ewt[i])^2 ),
#     * where N is the problem dimension.
#     *
#     * The return value of these functions is equal to CV_SUCCESS = 0 if
#     * there were no errors; otherwise it is a negative int equal to:
#     *   CV_MEM_NULL      indicating cvode_mem was NULL (i.e.,
#     *                    CVodeCreate has not been called).
#     *   CV_NO_MALLOC     indicating that cvode_mem has not been
#     *                    allocated (i.e., CVodeInit has not been
#     *                    called).
#     *   CV_ILL_INPUT     indicating an input argument was illegal
#     *                    (e.g. a negative tolerance)
#     * In case of an error return, an error message is also printed.
#     * -----------------------------------------------------------------
#     */
#    
    int CVodeSStolerances(void *cvode_mem, realtype reltol, realtype abstol)
    int CVodeSVtolerances(void *cvode_mem, realtype reltol, N_Vector abstol)
    int CVodeWFtolerances(void *cvode_mem, CVEwtFn efun)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Function : CVodeRootInit
#     * -----------------------------------------------------------------
#     * CVodeRootInit initializes a rootfinding problem to be solved
#     * during the integration of the ODE system.  It must be called
#     * after CVodeCreate, and before CVode.  The arguments are:
#     *
#     * cvode_mem = pointer to CVODE memory returned by CVodeCreate.
#     *
#     * nrtfn     = number of functions g_i, an int >= 0.
#     *
#     * g         = name of user-supplied function, of type CVRootFn,
#     *             defining the functions g_i whose roots are sought.
#     *
#     * If a new problem is to be solved with a call to CVodeReInit,
#     * where the new problem has no root functions but the prior one
#     * did, then call CVodeRootInit with nrtfn = 0.
#     *
#     * The return value of CVodeRootInit is CV_SUCCESS = 0 if there were
#     * no errors; otherwise it is a negative int equal to:
#     *   CV_MEM_NULL     indicating cvode_mem was NULL, or
#     *   CV_MEM_FAIL     indicating a memory allocation failed.
#     *                    (including an attempt to increase maxord).
#     *   CV_ILL_INPUT   indicating nrtfn > 0 but g = NULL.
#     * In case of an error return, an error message is also printed.
#     * -----------------------------------------------------------------
#     */
#    
    int CVodeRootInit(void *cvode_mem, int nrtfn, CVRootFn g)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Function : CVode
#     * -----------------------------------------------------------------
#     * CVode integrates the ODE over an interval in t.
#     * If itask is CV_NORMAL, then the solver integrates from its
#     * current internal t value to a point at or beyond tout, then
#     * interpolates to t = tout and returns y(tout) in the user-
#     * allocated vector yout. If itask is CV_ONE_STEP, then the solver
#     * takes one internal time step and returns in yout the value of
#     * y at the new internal time. In this case, tout is used only
#     * during the first call to CVode to determine the direction of
#     * integration and the rough scale of the t variable. If tstop is
#     * enabled (through a call to CVodeSetStopTime), then CVode returns
#     * the solution at tstop. Once the integrator returns at a tstop
#     * time, any future testing for tstop is disabled (and can be 
#     * reenabled only though a new call to CVodeSetStopTime).
#     * The time reached by the solver is placed in (*tret). The
#     * user is responsible for allocating the memory for this value.
#     *
#     * cvode_mem is the pointer to CVODE memory returned by
#     *           CVodeCreate.
#     *
#     * tout  is the next time at which a computed solution is desired.
#     *
#     * yout  is the computed solution vector. In CV_NORMAL mode with no
#     *       errors and no roots found, yout=y(tout).
#     *
#     * tret  is a pointer to a real location. CVode sets (*tret) to
#     *       the time reached by the solver and returns
#     *       yout=y(*tret).
#     *
#     * itask is CV_NORMAL or CV_ONE_STEP. These two modes are described above.
#     *
#     * Here is a brief description of each return value:
#     *
#     * CV_SUCCESS:      CVode succeeded and no roots were found.
#     *
#     * CV_ROOT_RETURN:  CVode succeeded, and found one or more roots.
#     *                  If nrtfn > 1, call CVodeGetRootInfo to see
#     *                  which g_i were found to have a root at (*tret).
#     *
#     * CV_TSTOP_RETURN: CVode succeeded and returned at tstop.
#     *
#     * CV_MEM_NULL:     The cvode_mem argument was NULL.
#     *
#     * CV_NO_MALLOC:    cvode_mem was not allocated.
#     *
#     * CV_ILL_INPUT:    One of the inputs to CVode is illegal. This
#     *                  includes the situation when a component of the
#     *                  error weight vectors becomes < 0 during
#     *                  internal time-stepping.  It also includes the
#     *                  situation where a root of one of the root
#     *                  functions was found both at t0 and very near t0.
#     *                  The ILL_INPUT flag will also be returned if the
#     *                  linear solver routine CV--- (called by the user
#     *                  after calling CVodeCreate) failed to set one of
#     *                  the linear solver-related fields in cvode_mem or
#     *                  if the linear solver's init routine failed. In
#     *                  any case, the user should see the printed
#     *                  error message for more details.
#     *
#     * CV_TOO_MUCH_WORK: The solver took mxstep internal steps but
#     *                  could not reach tout. The default value for
#     *                  mxstep is MXSTEP_DEFAULT = 500.
#     *
#     * CV_TOO_MUCH_ACC: The solver could not satisfy the accuracy
#     *                  demanded by the user for some internal step.
#     *
#     * CV_ERR_FAILURE:  Error test failures occurred too many times
#     *                  (= MXNEF = 7) during one internal time step or
#     *                  occurred with |h| = hmin.
#     *
#     * CV_CONV_FAILURE: Convergence test failures occurred too many
#     *                  times (= MXNCF = 10) during one internal time
#     *                  step or occurred with |h| = hmin.
#     *
#     * CV_LINIT_FAIL:   The linear solver's initialization function 
#     *                  failed.
#     *
#     * CV_LSETUP_FAIL:  The linear solver's setup routine failed in an
#     *                  unrecoverable manner.
#     *
#     * CV_LSOLVE_FAIL:  The linear solver's solve routine failed in an
#     *                  unrecoverable manner.
#     * -----------------------------------------------------------------
#     */
#    
    int CVode(void *cvode_mem, realtype tout, N_Vector yout,
    			  realtype *tret, int itask)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Function : CVodeGetDky
#     * -----------------------------------------------------------------
#     * CVodeGetDky computes the kth derivative of the y function at
#     * time t, where tn-hu <= t <= tn, tn denotes the current
#     * internal time reached, and hu is the last internal step size
#     * successfully used by the solver. The user may request
#     * k=0, 1, ..., qu, where qu is the order last used. The
#     * derivative vector is returned in dky. This vector must be
#     * allocated by the caller. It is only legal to call this
#     * function after a successful return from CVode.
#     *
#     * cvode_mem is the pointer to CVODE memory returned by
#     *           CVodeCreate.
#     *
#     * t   is the time at which the kth derivative of y is evaluated.
#     *     The legal range for t is [tn-hu,tn] as described above.
#     *
#     * k   is the order of the derivative of y to be computed. The
#     *     legal range for k is [0,qu] as described above.
#     *
#     * dky is the output derivative vector [((d/dy)^k)y](t).
#     *
#     * The return value for CVodeGetDky is one of:
#     *
#     *   CV_SUCCESS:  CVodeGetDky succeeded.
#     *
#     *   CV_BAD_K:    k is not in the range 0, 1, ..., qu.
#     *
#     *   CV_BAD_T:    t is not in the interval [tn-hu,tn].
#     *
#     *   CV_BAD_DKY:  The dky argument was NULL.
#     *
#     *   CV_MEM_NULL: The cvode_mem argument was NULL.
#     * -----------------------------------------------------------------
#     */
#    
    int CVodeGetDky(void *cvode_mem, realtype t, int k, N_Vector dky)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Integrator optional output extraction functions
#     * -----------------------------------------------------------------
#     * The following functions can be called to get optional outputs
#     * and statistics related to the main integrator.
#     * -----------------------------------------------------------------
#     * CVodeGetWorkSpace returns the CVODE real and integer workspaces
#     * CVodeGetNumSteps returns the cumulative number of internal
#     *                  steps taken by the solver
#     * CVodeGetNumRhsEvals returns the number of calls to the user's
#     *                     f function
#     * CVodeGetNumLinSolvSetups returns the number of calls made to
#     *                          the linear solver's setup routine
#     * CVodeGetNumErrTestFails returns the number of local error test
#     *                         failures that have occured
#     * CVodeGetLastOrder returns the order used during the last
#     *                   internal step
#     * CVodeGetCurrentOrder returns the order to be used on the next
#     *                      internal step
#     * CVodeGetNumStabLimOrderReds returns the number of order
#     *                             reductions due to stability limit
#     *                             detection
#     * CVodeGetActualInitStep returns the actual initial step size
#     *                        used by CVODE
#     * CVodeGetLastStep returns the step size for the last internal
#     *                  step
#     * CVodeGetCurrentStep returns the step size to be attempted on
#     *                     the next internal step
#     * CVodeGetCurrentTime returns the current internal time reached
#     *                     by the solver
#     * CVodeGetTolScaleFactor returns a suggested factor by which the
#     *                        user's tolerances should be scaled when
#     *                        too much accuracy has been requested for
#     *                        some internal step
#     * CVodeGetErrWeights returns the current error weight vector.
#     *                    The user must allocate space for eweight.
#     * CVodeGetEstLocalErrors returns the vector of estimated local
#     *                        errors. The user must allocate space
#     *                        for ele.
#     * CVodeGetNumGEvals returns the number of calls to the user's
#     *                   g function (for rootfinding)
#     * CVodeGetRootInfo returns the indices for which g_i was found to 
#     *                  have a root. The user must allocate space for 
#     *                  rootsfound. For i = 0 ... nrtfn-1, 
#     *                  rootsfound[i] = 1 if g_i has a root, and = 0 if not.
#     *
#     * CVodeGet* return values:
#     *   CV_SUCCESS   if succesful
#     *   CV_MEM_NULL  if the cvode memory was NULL
#     *   CV_NO_SLDET  if stability limit was not turned on
#     * -----------------------------------------------------------------
#     */
#    
    int CVodeGetWorkSpace(void *cvode_mem, long int *lenrw, long int *leniw)
    int CVodeGetNumSteps(void *cvode_mem, long int *nsteps)
    int CVodeGetNumRhsEvals(void *cvode_mem, long int *nfevals)
    int CVodeGetNumLinSolvSetups(void *cvode_mem, long int *nlinsetups)
    int CVodeGetNumErrTestFails(void *cvode_mem, long int *netfails)
    int CVodeGetLastOrder(void *cvode_mem, int *qlast)
    int CVodeGetCurrentOrder(void *cvode_mem, int *qcur)
    int CVodeGetNumStabLimOrderReds(void *cvode_mem, long int *nslred)
    int CVodeGetActualInitStep(void *cvode_mem, realtype *hinused)
    int CVodeGetLastStep(void *cvode_mem, realtype *hlast)
    int CVodeGetCurrentStep(void *cvode_mem, realtype *hcur)
    int CVodeGetCurrentTime(void *cvode_mem, realtype *tcur)
    int CVodeGetTolScaleFactor(void *cvode_mem, realtype *tolsfac)
    int CVodeGetErrWeights(void *cvode_mem, N_Vector eweight)
    int CVodeGetEstLocalErrors(void *cvode_mem, N_Vector ele)
    int CVodeGetNumGEvals(void *cvode_mem, long int *ngevals)
    int CVodeGetRootInfo(void *cvode_mem, int *rootsfound)
#    
#    /*
#     * -----------------------------------------------------------------
#     * As a convenience, the following functions provides the
#     * optional outputs in one group.
#     * -----------------------------------------------------------------
#     */
#    
    int CVodeGetIntegratorStats(void *cvode_mem, long int *nsteps,
    					    long int *nfevals, long int *nlinsetups,
    					    long int *netfails, int *qlast,
    					    int *qcur, realtype *hinused, realtype *hlast,
    					    realtype *hcur, realtype *tcur)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Nonlinear solver optional output extraction functions
#     * -----------------------------------------------------------------
#     * The following functions can be called to get optional outputs
#     * and statistics related to the nonlinear solver.
#     * -----------------------------------------------------------------
#     * CVodeGetNumNonlinSolvIters returns the number of nonlinear
#     *                            solver iterations performed.
#     * CVodeGetNumNonlinSolvConvFails returns the number of nonlinear
#     *                                convergence failures.
#     * -----------------------------------------------------------------
#     */
#    
    int CVodeGetNumNonlinSolvIters(void *cvode_mem, long int *nniters)
    int CVodeGetNumNonlinSolvConvFails(void *cvode_mem, long int *nncfails)
#    
#    /*
#     * -----------------------------------------------------------------
#     * As a convenience, the following function provides the
#     * nonlinear solver optional outputs in a group.
#     * -----------------------------------------------------------------
#     */
#    
    int CVodeGetNonlinSolvStats(void *cvode_mem, long int *nniters,
    					    long int *nncfails)
#    
#    /*
#     * -----------------------------------------------------------------
#     * The following function returns the name of the constant 
#     * associated with a CVODE return flag
#     * -----------------------------------------------------------------
#     */
#    
    char *CVodeGetReturnFlagName(long int flag)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Function : CVodeFree
#     * -----------------------------------------------------------------
#     * CVodeFree frees the problem memory cvode_mem allocated by
#     * CVodeCreate and CVodeInit. Its only argument is the pointer
#     * cvode_mem returned by CVodeCreate.
#     * -----------------------------------------------------------------
#     */
#    
    void CVodeFree(void **cvode_mem)
#    
#    #ifdef __cplusplus
#    }
#    #endif
#    
#    #endif

cdef extern from "cvode/cvode_dense.h":
    int CVDense(void *cvode_mem, long int N)
    
cdef extern from "cvode/cvode_band.h":
#    /*
#     * -----------------------------------------------------------------
#     * Function : CVBand
#     * -----------------------------------------------------------------
#     * A call to the CVBand function links the main CVODE integrator
#     * with the CVBAND linear solver.
#     *
#     * cvode_mem is the pointer to the integrator memory returned by
#     *           CVodeCreate.
#     *
#     * N is the size of the ODE system.
#     *
#     * mupper is the upper bandwidth of the band Jacobian
#     *        approximation.
#     *
#     * mlower is the lower bandwidth of the band Jacobian
#     *        approximation.
#     *
#     * The return value of CVBand is one of:
#     *    CVDLS_SUCCESS   if successful
#     *    CVDLS_MEM_NULL  if the cvode memory was NULL
#     *    CVDLS_MEM_FAIL  if there was a memory allocation failure
#     *    CVDLS_ILL_INPUT if a required vector operation is missing or
#     *                       if a bandwidth has an illegal value.
#     * -----------------------------------------------------------------
#     */

    int CVBand(void *cvode_mem, long int N, long int mupper, long int mlower)
    
cdef extern from "cvode/cvode_direct.h":
#/*
# * =================================================================
# *              C V D I R E C T     C O N S T A N T S
# * =================================================================
# */
#
#/* 
# * -----------------------------------------------------------------
# * CVDLS return values 
# * -----------------------------------------------------------------
# */
#
    enum: CVDLS_SUCCES
    enum: CVDLS_MEM_NULL
    enum: CVDLS_LMEM_NULL
    enum: CVDLS_ILL_INPUT
    enum: CVDLS_MEM_FAIL
#
#/* Additional last_flag values */
#
    enum: CVDLS_JACFUNC_UNRECVR
    enum: CVDLS_JACFUNC_RECVR
#
#/*
# * =================================================================
# *              F U N C T I O N   T Y P E S
# * =================================================================
# */
#
#/*
# * -----------------------------------------------------------------
# * Type: CVDlsDenseJacFn
# * -----------------------------------------------------------------
# *
# * A dense Jacobian approximation function Jac must be of type 
# * CVDlsDenseJacFn. Its parameters are:
# *
# * N   is the problem size.
# *
# * Jac is the dense matrix (of type DlsMat) that will be loaded
# *     by a CVDlsDenseJacFn with an approximation to the Jacobian 
# *     matrix J = (df_i/dy_j) at the point (t,y). 
# *
# * t   is the current value of the independent variable.
# *
# * y   is the current value of the dependent variable vector,
# *     namely the predicted value of y(t).
# *
# * fy  is the vector f(t,y).
# *
# * user_data is a pointer to user data - the same as the user_data
# *     parameter passed to CVodeSetFdata.
# *
# * tmp1, tmp2, and tmp3 are pointers to memory allocated for
# * vectors of length N which can be used by a CVDlsDenseJacFn
# * as temporary storage or work space.
# *
# * A CVDlsDenseJacFn should return 0 if successful, a positive 
# * value if a recoverable error occurred, and a negative value if 
# * an unrecoverable error occurred.
# *
# * -----------------------------------------------------------------
# *
# * NOTE: The following are two efficient ways to load a dense Jac:         
# * (1) (with macros - no explicit data structure references)      
# *     for (j=0; j < Neq; j++) {                                  
# *       col_j = DENSE_COL(Jac,j);                                 
# *       for (i=0; i < Neq; i++) {                                
# *         generate J_ij = the (i,j)th Jacobian element           
# *         col_j[i] = J_ij;                                       
# *       }                                                        
# *     }                                                          
# * (2) (without macros - explicit data structure references)      
# *     for (j=0; j < Neq; j++) {                                  
# *       col_j = (Jac->data)[j];                                   
# *       for (i=0; i < Neq; i++) {                                
# *         generate J_ij = the (i,j)th Jacobian element           
# *         col_j[i] = J_ij;                                       
# *       }                                                        
# *     }                                                          
# * A third way, using the DENSE_ELEM(A,i,j) macro, is much less   
# * efficient in general.  It is only appropriate for use in small 
# * problems in which efficiency of access is NOT a major concern. 
# *                                                                
# * NOTE: If the user's Jacobian routine needs other quantities,   
# *     they are accessible as follows: hcur (the current stepsize)
# *     and ewt (the error weight vector) are accessible through   
# *     CVodeGetCurrentStep and CVodeGetErrWeights, respectively 
# *     (see cvode.h). The unit roundoff is available as 
# *     UNIT_ROUNDOFF defined in sundials_types.h.
# *
# * -----------------------------------------------------------------
# */
#  
#  
    ctypedef int (*CVDlsDenseJacFn)(long int N, realtype t,
			       N_Vector y, N_Vector fy, 
			       DlsMat Jac, void *user_data,
			       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
#  
#/*
# * -----------------------------------------------------------------
# * Type: CVDlsBandJacFn
# * -----------------------------------------------------------------
# *
# * A band Jacobian approximation function Jac must have the
# * prototype given below. Its parameters are:
# *
# * N is the length of all vector arguments.
# *
# * mupper is the upper half-bandwidth of the approximate banded
# * Jacobian. This parameter is the same as the mupper parameter
# * passed by the user to the linear solver initialization function.
# *
# * mlower is the lower half-bandwidth of the approximate banded
# * Jacobian. This parameter is the same as the mlower parameter
# * passed by the user to the linear solver initialization function.
# *
# * t is the current value of the independent variable.
# *
# * y is the current value of the dependent variable vector,
# *      namely the predicted value of y(t).
# *
# * fy is the vector f(t,y).
# *
# * Jac is the band matrix (of type DlsMat) that will be loaded
# * by a CVDlsBandJacFn with an approximation to the Jacobian matrix
# * Jac = (df_i/dy_j) at the point (t,y).
# * Three efficient ways to load J are:
# *
# * (1) (with macros - no explicit data structure references)
# *    for (j=0; j < n; j++) {
# *       col_j = BAND_COL(Jac,j);
# *       for (i=j-mupper; i <= j+mlower; i++) {
# *         generate J_ij = the (i,j)th Jacobian element
# *         BAND_COL_ELEM(col_j,i,j) = J_ij;
# *       }
# *     }
# *
# * (2) (with BAND_COL macro, but without BAND_COL_ELEM macro)
# *    for (j=0; j < n; j++) {
# *       col_j = BAND_COL(Jac,j);
# *       for (k=-mupper; k <= mlower; k++) {
# *         generate J_ij = the (i,j)th Jacobian element, i=j+k
# *         col_j[k] = J_ij;
# *       }
# *     }
# *
# * (3) (without macros - explicit data structure references)
# *     offset = Jac->smu;
# *     for (j=0; j < n; j++) {
# *       col_j = ((Jac->data)[j])+offset;
# *       for (k=-mupper; k <= mlower; k++) {
# *         generate J_ij = the (i,j)th Jacobian element, i=j+k
# *         col_j[k] = J_ij;
# *       }
# *     }
# * Caution: Jac->smu is generally NOT the same as mupper.
# *
# * The BAND_ELEM(A,i,j) macro is appropriate for use in small
# * problems in which efficiency of access is NOT a major concern.
# *
# * user_data is a pointer to user data - the same as the user_data
# *          parameter passed to CVodeSetFdata.
# *
# * NOTE: If the user's Jacobian routine needs other quantities,
# *     they are accessible as follows: hcur (the current stepsize)
# *     and ewt (the error weight vector) are accessible through
# *     CVodeGetCurrentStep and CVodeGetErrWeights, respectively
# *     (see cvode.h). The unit roundoff is available as
# *     UNIT_ROUNDOFF defined in sundials_types.h
# *
# * tmp1, tmp2, and tmp3 are pointers to memory allocated for
# * vectors of length N which can be used by a CVDlsBandJacFn
# * as temporary storage or work space.
# *
# * A CVDlsBandJacFn should return 0 if successful, a positive value
# * if a recoverable error occurred, and a negative value if an 
# * unrecoverable error occurred.
# * -----------------------------------------------------------------
# */
#
    ctypedef int (*CVDlsBandJacFn)(long int N, long int mupper, long int mlower,
			      realtype t, N_Vector y, N_Vector fy, 
			      DlsMat Jac, void *user_data,
			      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
#
#/*
# * =================================================================
# *            E X P O R T E D    F U N C T I O N S 
# * =================================================================
# */
#
#/*
# * -----------------------------------------------------------------
# * Optional inputs to the CVDLS linear solver
# * -----------------------------------------------------------------
# *
# * CVDlsSetDenseJacFn specifies the dense Jacobian approximation
# * routine to be used for a direct dense linear solver.
# *
# * CVDlsSetBandJacFn specifies the band Jacobian approximation
# * routine to be used for a direct band linear solver.
# *
# * By default, a difference quotient approximation, supplied with
# * the solver is used.
# *
# * The return value is one of:
# *    CVDLS_SUCCESS   if successful
# *    CVDLS_MEM_NULL  if the CVODE memory was NULL
# *    CVDLS_LMEM_NULL if the linear solver memory was NULL
# * -----------------------------------------------------------------
# */
#
    int CVDlsSetDenseJacFn(void *cvode_mem, CVDlsDenseJacFn jac)
    int CVDlsSetBandJacFn(void *cvode_mem, CVDlsBandJacFn jac)
#
#/*
# * -----------------------------------------------------------------
# * Optional outputs from the CVDLS linear solver
# * -----------------------------------------------------------------
# *
# * CVDlsGetWorkSpace   returns the real and integer workspace used
# *                     by the direct linear solver.
# * CVDlsGetNumJacEvals returns the number of calls made to the
# *                     Jacobian evaluation routine jac.
# * CVDlsGetNumRhsEvals returns the number of calls to the user
# *                     f routine due to finite difference Jacobian
# *                     evaluation.
# * CVDlsGetLastFlag    returns the last error flag set by any of
# *                     the CVDLS interface functions.
# *
# * The return value of CVDlsGet* is one of:
# *    CVDLS_SUCCESS   if successful
# *    CVDLS_MEM_NULL  if the CVODE memory was NULL
# *    CVDLS_LMEM_NULL if the linear solver memory was NULL
# * -----------------------------------------------------------------
# */
#
    int CVDlsGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS)
    int CVDlsGetNumJacEvals(void *cvode_mem, long int *njevals)
    int CVDlsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS)
    int CVDlsGetLastFlag(void *cvode_mem, long int *flag)
#
#/*
# * -----------------------------------------------------------------
# * The following function returns the name of the constant 
# * associated with a CVDLS return flag
# * -----------------------------------------------------------------
# */
#
#SUNDIALS_EXPORT char *CVDlsGetReturnFlagName(long int flag);    
    
cdef extern from "cvode/cvode_spils.h":
#    /*
#     * -----------------------------------------------------------------
#     * CVSPILS return values 
#     * -----------------------------------------------------------------
#     */
    
    enum: CVSPILS_SUCCESS
    enum: CVSPILS_MEM_NULL
    enum: CVSPILS_LMEM_NULL
    enum: CVSPILS_ILL_INPUT
    enum: CVSPILS_MEM_FAIL
    enum: CVSPILS_PMEM_NULL
    
#    /*
#     * -----------------------------------------------------------------
#     * CVSPILS solver constants
#     * -----------------------------------------------------------------
#     * CVSPILS_MAXL   : default value for the maximum Krylov
#     *                  dimension
#     *
#     * CVSPILS_MSBPRE : maximum number of steps between
#     *                  preconditioner evaluations
#     *
#     * CVSPILS_DGMAX  : maximum change in gamma between
#     *                  preconditioner evaluations
#     *
#     * CVSPILS_EPLIN  : default value for factor by which the
#     *                  tolerance on the nonlinear iteration is
#     *                  multiplied to get a tolerance on the linear
#     *                  iteration
#     * -----------------------------------------------------------------
#     */
#    
    enum: CVSPILS_MAXL
    enum: CVSPILS_MSBPRE
    enum: CVSPILS_DGMAX
    enum: CVSPILS_EPLIN
        
        
#    /*
#     * -----------------------------------------------------------------
#     * Type : CVSpilsPrecSetupFn
#     * -----------------------------------------------------------------
#     * The user-supplied preconditioner setup function PrecSetup and
#     * the user-supplied preconditioner solve function PrecSolve
#     * together must define left and right preconditoner matrices
#     * P1 and P2 (either of which may be trivial), such that the
#     * product P1*P2 is an approximation to the Newton matrix
#     * M = I - gamma*J.  Here J is the system Jacobian J = df/dy,
#     * and gamma is a scalar proportional to the integration step
#     * size h.  The solution of systems P z = r, with P = P1 or P2,
#     * is to be carried out by the PrecSolve function, and PrecSetup
#     * is to do any necessary setup operations.
#     *
#     * The user-supplied preconditioner setup function PrecSetup
#     * is to evaluate and preprocess any Jacobian-related data
#     * needed by the preconditioner solve function PrecSolve.
#     * This might include forming a crude approximate Jacobian,
#     * and performing an LU factorization on the resulting
#     * approximation to M.  This function will not be called in
#     * advance of every call to PrecSolve, but instead will be called
#     * only as often as necessary to achieve convergence within the
#     * Newton iteration.  If the PrecSolve function needs no
#     * preparation, the PrecSetup function can be NULL.
#     *
#     * For greater efficiency, the PrecSetup function may save
#     * Jacobian-related data and reuse it, rather than generating it
#     * from scratch.  In this case, it should use the input flag jok
#     * to decide whether to recompute the data, and set the output
#     * flag *jcurPtr accordingly.
#     *
#     * Each call to the PrecSetup function is preceded by a call to
#     * the RhsFn f with the same (t,y) arguments.  Thus the PrecSetup
#     * function can use any auxiliary data that is computed and
#     * saved by the f function and made accessible to PrecSetup.
#     *
#     * A function PrecSetup must have the prototype given below.
#     * Its parameters are as follows:
#     *
#     * t       is the current value of the independent variable.
#     *
#     * y       is the current value of the dependent variable vector,
#     *          namely the predicted value of y(t).
#     *
#     * fy      is the vector f(t,y).
#     *
#     * jok     is an input flag indicating whether Jacobian-related
#     *         data needs to be recomputed, as follows:
#     *           jok == FALSE means recompute Jacobian-related data
#     *                  from scratch.
#     *           jok == TRUE  means that Jacobian data, if saved from
#     *                  the previous PrecSetup call, can be reused
#     *                  (with the current value of gamma).
#     *         A Precset call with jok == TRUE can only occur after
#     *         a call with jok == FALSE.
#     *
#     * jcurPtr is a pointer to an output integer flag which is
#     *         to be set by PrecSetup as follows:
#     *         Set *jcurPtr = TRUE if Jacobian data was recomputed.
#     *         Set *jcurPtr = FALSE if Jacobian data was not recomputed,
#     *                        but saved data was reused.
#     *
#     * gamma   is the scalar appearing in the Newton matrix.
#     *
#     * user_data  is a pointer to user data - the same as the user_data
#     *         parameter passed to the CVodeSetUserData function.
#     *
#     * tmp1, tmp2, and tmp3 are pointers to memory allocated
#     *                      for N_Vectors which can be used by
#     *                      CVSpilsPrecSetupFn as temporary storage or
#     *                      work space.
#     *
#     * NOTE: If the user's preconditioner needs other quantities,
#     *       they are accessible as follows: hcur (the current stepsize)
#     *       and ewt (the error weight vector) are accessible through
#     *       CVodeGetCurrentStep and CVodeGetErrWeights, respectively).
#     *       The unit roundoff is available as UNIT_ROUNDOFF defined in
#     *       sundials_types.h.
#     *
#     * Returned value:
#     * The value to be returned by the PrecSetup function is a flag
#     * indicating whether it was successful.  This value should be
#     *   0   if successful,
#     *   > 0 for a recoverable error (step will be retried),
#     *   < 0 for an unrecoverable error (integration is halted).
#     * -----------------------------------------------------------------
#     */
#    
    ctypedef int (*CVSpilsPrecSetupFn)(realtype t, N_Vector y, N_Vector fy,
                                      booleantype jok, booleantype *jcurPtr,
                                      realtype gamma, void *user_data,
                                      N_Vector tmp1, N_Vector tmp2,
                                      N_Vector tmp3)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Type : CVSpilsPrecSolveFn
#     * -----------------------------------------------------------------
#     * The user-supplied preconditioner solve function PrecSolve
#     * is to solve a linear system P z = r in which the matrix P is
#     * one of the preconditioner matrices P1 or P2, depending on the
#     * type of preconditioning chosen.
#     *
#     * A function PrecSolve must have the prototype given below.
#     * Its parameters are as follows:
#     *
#     * t      is the current value of the independent variable.
#     *
#     * y      is the current value of the dependent variable vector.
#     *
#     * fy     is the vector f(t,y).
#     *
#     * r      is the right-hand side vector of the linear system.
#     *
#     * z      is the output vector computed by PrecSolve.
#     *
#     * gamma  is the scalar appearing in the Newton matrix.
#     *
#     * delta  is an input tolerance for use by PSolve if it uses
#     *        an iterative method in its solution.  In that case,
#     *        the residual vector Res = r - P z of the system
#     *        should be made less than delta in weighted L2 norm,
#     *        i.e., sqrt [ Sum (Res[i]*ewt[i])^2 ] < delta.
#     *        Note: the error weight vector ewt can be obtained
#     *        through a call to the routine CVodeGetErrWeights.
#     *
#     * lr     is an input flag indicating whether PrecSolve is to use
#     *        the left preconditioner P1 or right preconditioner
#     *        P2: lr = 1 means use P1, and lr = 2 means use P2.
#     *
#     * user_data  is a pointer to user data - the same as the user_data
#     *         parameter passed to the CVodeSetUserData function.
#     *
#     * tmp    is a pointer to memory allocated for an N_Vector
#     *        which can be used by PSolve for work space.
#     *
#     * Returned value:
#     * The value to be returned by the PrecSolve function is a flag
#     * indicating whether it was successful.  This value should be
#     *   0 if successful,
#     *   positive for a recoverable error (step will be retried),
#     *   negative for an unrecoverable error (integration is halted).
#     * -----------------------------------------------------------------
#     */
#    
    ctypedef int (*CVSpilsPrecSolveFn)(realtype t, N_Vector y, N_Vector fy,
                                      N_Vector r, N_Vector z,
                                      realtype gamma, realtype delta,
                                      int lr, void *user_data, N_Vector tmp)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Type : CVSpilsJacTimesVecFn
#     * -----------------------------------------------------------------
#     * The user-supplied function jtimes is to generate the product
#     * J*v for given v, where J is the Jacobian df/dy, or an
#     * approximation to it, and v is a given vector. It should return
#     * 0 if successful a positive value for a recoverable error or 
#     * a negative value for an unrecoverable failure.
#     *
#     * A function jtimes must have the prototype given below. Its
#     * parameters are as follows:
#     *
#     *   v        is the N_Vector to be multiplied by J.
#     *
#     *   Jv       is the output N_Vector containing J*v.
#     *
#     *   t        is the current value of the independent variable.
#     *
#     *   y        is the current value of the dependent variable
#     *            vector.
#     *
#     *   fy       is the vector f(t,y).
#     *
#     *   user_data   is a pointer to user data, the same as the user_data
#     *            parameter passed to the CVodeSetUserData function.
#     *
#     *   tmp      is a pointer to memory allocated for an N_Vector
#     *            which can be used by Jtimes for work space.
#     * -----------------------------------------------------------------
#     */
#    
    ctypedef int (*CVSpilsJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t,
                                        N_Vector y, N_Vector fy,
                                        void *user_data, N_Vector tmp)
#    
#    
#    
#    /*
#     * -----------------------------------------------------------------
#     * Optional inputs to the CVSPILS linear solver
#     * -----------------------------------------------------------------
#     *
#     * CVSpilsSetPrecType resets the type of preconditioner, pretype,
#     *                from the value previously set.
#     *                This must be one of PREC_NONE, PREC_LEFT, 
#     *                PREC_RIGHT, or PREC_BOTH.
#     *
#     * CVSpilsSetGSType specifies the type of Gram-Schmidt
#     *                orthogonalization to be used. This must be one of
#     *                the two enumeration constants MODIFIED_GS or
#     *                CLASSICAL_GS defined in iterative.h. These correspond
#     *                to using modified Gram-Schmidt and classical
#     *                Gram-Schmidt, respectively.
#     *                Default value is MODIFIED_GS.
#     *
#     * CVSpilsSetMaxl resets the maximum Krylov subspace size, maxl,
#     *                from the value previously set.
#     *                An input value <= 0, gives the default value.
#     *
#     * CVSpilsSetEpsLin specifies the factor by which the tolerance on
#     *                the nonlinear iteration is multiplied to get a
#     *                tolerance on the linear iteration.
#     *                Default value is 0.05.
#     *
#     * CVSpilsSetPreconditioner specifies the PrecSetup and PrecSolve functions.
#     *                Default is NULL for both arguments (no preconditioning)
#     *
#     * CVSpilsSetJacTimesVecFn specifies the jtimes function. Default is to 
#     *                use an internal finite difference approximation routine.
#     *
#     * The return value of CVSpilsSet* is one of:
#     *    CVSPILS_SUCCESS   if successful
#     *    CVSPILS_MEM_NULL  if the cvode memory was NULL
#     *    CVSPILS_LMEM_NULL if the linear solver memory was NULL
#     *    CVSPILS_ILL_INPUT if an input has an illegal value
#     * -----------------------------------------------------------------
#     */
#    
    int CVSpilsSetPrecType(void *cvode_mem, int pretype)
    int CVSpilsSetGSType(void *cvode_mem, int gstype)
    int CVSpilsSetMaxl(void *cvode_mem, int maxl)
    int CVSpilsSetEpsLin(void *cvode_mem, realtype eplifac)
    int CVSpilsSetPreconditioner(void *cvode_mem, 
                                                 CVSpilsPrecSetupFn pset,
                                                 CVSpilsPrecSolveFn psolve)
    int CVSpilsSetJacTimesVecFn(void *cvode_mem, 
                                                CVSpilsJacTimesVecFn jtv)
#    
#    /*
#     * -----------------------------------------------------------------
#     * Optional outputs from the CVSPILS linear solver
#     * -----------------------------------------------------------------
#     * CVSpilsGetWorkSpace returns the real and integer workspace used
#     *                by the SPILS module.
#     *
#     * CVSpilsGetNumPrecEvals returns the number of preconditioner
#     *                 evaluations, i.e. the number of calls made
#     *                 to PrecSetup with jok==FALSE.
#     *
#     * CVSpilsGetNumPrecSolves returns the number of calls made to
#     *                 PrecSolve.
#     *
#     * CVSpilsGetNumLinIters returns the number of linear iterations.
#     *
#     * CVSpilsGetNumConvFails returns the number of linear
#     *                 convergence failures.
#     *
#     * CVSpilsGetNumJtimesEvals returns the number of calls to jtimes.
#     *
#     * CVSpilsGetNumRhsEvals returns the number of calls to the user
#     *                 f routine due to finite difference Jacobian
#     *                 times vector evaluation.
#     *
#     * CVSpilsGetLastFlag returns the last error flag set by any of
#     *                 the CVSPILS interface functions.
#     *
#     * The return value of CVSpilsGet* is one of:
#     *    CVSPILS_SUCCESS   if successful
#     *    CVSPILS_MEM_NULL  if the cvode memory was NULL
#     *    CVSPILS_LMEM_NULL if the linear solver memory was NULL
#     * -----------------------------------------------------------------
#     */
#    
    int CVSpilsGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS)
    int CVSpilsGetNumPrecEvals(void *cvode_mem, long int *npevals)
    int CVSpilsGetNumPrecSolves(void *cvode_mem, long int *npsolves)
    int CVSpilsGetNumLinIters(void *cvode_mem, long int *nliters)
    int CVSpilsGetNumConvFails(void *cvode_mem, long int *nlcfails)
    int CVSpilsGetNumJtimesEvals(void *cvode_mem, long int *njvevals)
    int CVSpilsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS)
    int CVSpilsGetLastFlag(void *cvode_mem, long int *flag)
#    
#    /*
#     * -----------------------------------------------------------------
#     * The following function returns the name of the constant 
#     * associated with a CVSPILS return flag
#     * -----------------------------------------------------------------
#     */

#SUNDIALS_EXPORT char *CVSpilsGetReturnFlagName(long int flag);
        
cdef extern from "cvode/cvode_spgmr.h":
#    /*
#     * -----------------------------------------------------------------
#     * Function : CVSpgmr
#     * -----------------------------------------------------------------
#     * A call to the CVSpgmr function links the main CVODE integrator
#     * with the CVSPGMR linear solver.
#     *
#     * cvode_mem is the pointer to the integrator memory returned by
#     *           CVodeCreate.
#     *
#     * pretype   is the type of user preconditioning to be done.
#     *           This must be one of the four enumeration constants
#     *           PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined 
#     *           in sundials_iterative.h.
#     *           These correspond to no preconditioning,
#     *           left preconditioning only, right preconditioning
#     *           only, and both left and right preconditioning,
#     *           respectively.
#     *
#     * maxl      is the maximum Krylov dimension. This is an
#     *           optional input to the CVSPGMR solver. Pass 0 to
#     *           use the default value CVSPGMR_MAXL=5.
#     *
#     * The return value of CVSpgmr is one of:
#     *    CVSPILS_SUCCESS   if successful
#     *    CVSPILS_MEM_NULL  if the cvode memory was NULL
#     *    CVSPILS_MEM_FAIL  if there was a memory allocation failure
#     *    CVSPILS_ILL_INPUT if a required vector operation is missing
#     * The above constants are defined in cvode_spils.h
#     *
#     * -----------------------------------------------------------------
#     */
#    
    int CVSpgmr(void *cvode_mem, int pretype, int maxl)
    
    
cdef extern from "cvode/cvode_spbcgs.h":
#    /*
#     * -----------------------------------------------------------------
#     * Function : CVSpbcg
#     * -----------------------------------------------------------------
#     * A call to the CVSpbcg function links the main CVODE integrator
#     * with the CVSPBCG linear solver.
#     *
#     * cvode_mem is the pointer to the integrator memory returned by
#     *           CVodeCreate.
#     *
#     * pretype   is the type of user preconditioning to be done.
#     *           This must be one of the four enumeration constants
#     *           PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined
#     *           in iterative.h. These correspond to no preconditioning,
#     *           left preconditioning only, right preconditioning
#     *           only, and both left and right preconditioning,
#     *           respectively.
#     *
#     * maxl      is the maximum Krylov dimension. This is an
#     *           optional input to the CVSPBCG solver. Pass 0 to
#     *           use the default value CVSPBCG_MAXL=5.
#     *
#     * The return value of CVSpbcg is one of:
#     *    CVSPILS_SUCCESS   if successful
#     *    CVSPILS_MEM_NULL  if the cvode memory was NULL
#     *    CVSPILS_MEM_FAIL  if there was a memory allocation failure
#     *    CVSPILS_ILL_INPUT if a required vector operation is missing
#     * The above constants are defined in cvode_spils.h
#     *
#     * -----------------------------------------------------------------
#     */
#    
    int CVSpbcg(void *cvode_mem, int pretype, int maxl)
    
cdef extern from "cvode/cvode_sptfqmr.h":    
#    /*
#     * -----------------------------------------------------------------
#     * Function : CVSptfqmr
#     * -----------------------------------------------------------------
#     * A call to the CVSptfqmr function links the main CVODE integrator
#     * with the CVSPTFQMR linear solver.
#     *
#     * cvode_mem is the pointer to the integrator memory returned by
#     *           CVodeCreate.
#     *
#     * pretype   is the type of user preconditioning to be done.
#     *           This must be one of the four enumeration constants
#     *           PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined
#     *           in iterative.h. These correspond to no preconditioning,
#     *           left preconditioning only, right preconditioning
#     *           only, and both left and right preconditioning,
#     *           respectively.
#     *
#     * maxl      is the maximum Krylov dimension. This is an
#     *           optional input to the CVSPTFQMR solver. Pass 0 to
#     *           use the default value CVSPILS_MAXL=5.
#     *
#     * The return value of CVSptfqmr is one of:
#     *    CVSPILS_SUCCESS   if successful
#     *    CVSPILS_MEM_NULL  if the cvode memory was NULL
#     *    CVSPILS_MEM_FAIL  if there was a memory allocation failure
#     *    CVSPILS_ILL_INPUT if a required vector operation is missing
#     * The above constants are defined in cvode_spils.h
#     *
#     * -----------------------------------------------------------------
#     */
#    
    int CVSptfqmr(void *cvode_mem, int pretype, int maxl)