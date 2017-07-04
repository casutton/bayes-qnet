# Pyrex binding to the GLPK library
#  Designed to be simple and efficient

cdef extern from "glpk.h":

    ctypedef struct glp_prob

    glp_prob *lpx_read_model(char *model, char *data, char *output)
    int lpx_simplex(glp_prob *lp)

    int lpx_get_status(glp_prob *lp)

    int lpx_get_num_cols(glp_prob *lp)
    char *lpx_get_col_name(glp_prob *lp, int j)
    double lpx_get_col_prim(glp_prob *lp, int j)

    enum:
      LPX_E_OK       = 200  # /* success */
      LPX_E_EMPTY    = 201  # /* empty problem */
      LPX_E_BADB     = 202  # /* invalid initial basis */
      LPX_E_INFEAS   = 203  # /* infeasible initial solution */
      LPX_E_FAULT    = 204  # /* unable to start the search */
      LPX_E_OBJLL    = 205  # /* objective lower limit reached */
      LPX_E_OBJUL    = 206  # /* objective upper limit reached */
      LPX_E_ITLIM    = 207  # /* iterations limit exhausted */
      LPX_E_TMLIM    = 208  # /* time limit exhausted */
      LPX_E_NOFEAS   = 209  # /* no feasible solution */
      LPX_E_INSTAB   = 210  # /* numerical instability */
      LPX_E_SING     = 211  # /* problems with basis matrix */
      LPX_E_NOCONV   = 212  # /* no convergence (interior) */
      LPX_E_NOPFS    = 213  # /* no primal feas. sol. (LP presolver) */
      LPX_E_NODFS    = 214  # /* no dual feas. sol. (LP presolver) */

    enum:
       LPX_OPT       =  180  # /* optimal */
       LPX_FEAS      =  181  # /* feasible */
       LPX_INFEAS    =  182  # /* infeasible */
       LPX_NOFEAS    =  183  # /* no feasible */
       LPX_UNBND     =  184  # /* unbounded */
       LPX_UNDEF     =  185  # /* undefined */



def solve_lp (fname):
    """Returns a solution to a problem in the GNU MathProb language.
        Uses the simplex method
        FNAME is the name of the problem file.
        Returns a hash where the key is the name of the variable in the file, 
        and the value is the value in the solution.
       If an optimal solution was not found, returns None."""
    print "Reading model from %s ... " % fname

    cdef glp_prob *lp
    lp = lpx_read_model (fname, NULL, NULL)
    result = lpx_simplex (lp)
    if result != LPX_E_OK:
        print "Error in solver: %d" % result
        return None

    status = lpx_get_status (lp)
    if status != LPX_OPT and status != LPX_FEAS:
        raise Exception ("GLPK error: Solution not good %d" % result)

    print "GLPK solution: result %d status %d" % (result, status)

    sol = dict()

    nvars = lpx_get_num_cols(lp)
    for i from 1 <= i <= nvars:
        name = lpx_get_col_name (lp, i)
        val = lpx_get_col_prim (lp, i)
        sol[name] = val

    return sol
    
