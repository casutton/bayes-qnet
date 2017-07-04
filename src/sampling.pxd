from pwfun cimport Pwfun

cdef extern from "randomkit.h":
    ctypedef struct rk_state

cdef double _slice (double x0, g, double lower, double upper) except *
cdef double _slice_pwfun (Pwfun fn, double x_initial) except *

cdef int _roll_die_unnorm (p)

cdef rk_state *_state()

cdef double rand ()
cdef double uniform (double L, double U)
cdef double exponential (double scale)

