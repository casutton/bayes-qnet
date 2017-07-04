cdef class Pwfun:
    cdef readonly object xs
    cdef readonly object derivs
    cdef readonly object fns
    cdef int N
    
    cdef object OOR

    cdef double value (self, double x)
    cdef int _find_knot (self, double x)
    cdef double L (self)
    cdef double U (self)
    
    

cdef class Pwlin:
    cdef object xs      # Points at which slope changes xs[0] is lower bound xs[-1] upper
    cdef object x_tan   # Points at which the tangent is taken.
    cdef object heights # f(x_tan)
    cdef object derivs  # df(x_tan)

    cdef object f
    cdef object fprime

    cdef double value (self, double x)
    cdef c_argmin (self, double x0)
