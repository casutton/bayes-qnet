from numpy import random
import numpy

from scipy import integrate

DEBUG = 0
cdef int counter = 0

def mcint (fn, a, b, iter=100):
    return _mcint (fn, a, b, iter)

def integrate_by_points (fn, points):
    total = 0
    for l,u in zip(points, points[1:]):
        total += integrate.quad (fn, l, u, full_output=DEBUG)[0]
    return total

# Wrapper over scipy quadratue method    
cdef double _integrate (fn, double a, double b):
    qtup = integrate.quad (fn, a, b, limit=1000, full_output=DEBUG)
    value = qtup[0]

    if DEBUG:
        if len(qtup) > 3:
            print "Integration error:"
            print qtup[2:]
            debug_write_function (fn, a, b)
        
    if not numpy.isnan(value):
        print qtup
        return value
    else:
        # last resort
        print "Warning: Using last-resort Monte Carlo integral"
        return _mcint (fn, a, b, 100)
    
cdef double _mcint (fn, double a, double b, int max_n):
    cdef double the_sum = 0
    cdef int N
    
    for N from 1 <= N < max_n:
        x = random.uniform(a,b)
        the_sum += fn(x)
        
    return (b-a)*the_sum/N
    
cdef void debug_write_function (fn, double a, double b):
    global counter
    
    if numpy.isinf(a) or numpy.isinf(b):
        print "Warning: can't write debugging info for range %f ... %f" % (a,b)
        return
    cdef int N = 1000
    x = a
    eps = (b-a)/N
    
    print counter
    f = open ("misc_debug_%d.txt" % counter, "w")
    for i from 0 <= i < N:
        f.write ("%.17f %.17f\n" % (x, fn(x)))
        x += eps
    f.write ("%.17f %.17f\n" % (b, fn(b)))
    f.close ()
    
    counter += 1
