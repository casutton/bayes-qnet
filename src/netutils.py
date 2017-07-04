from math import exp, sqrt
from numpy import random
from scipy import special
import numpy
import sampling
import pwfun

# This file is really annoying.  Various function that aren't in Cython b/c
#  they're easier to do with closures

def expify (fn, C=0): return lambda x: exp(fn(x)+C)

def all_true (*args): return True

def not_in_set2 (inc_set):
    return lambda a,e: e.tid not in inc_set

# uses Miller(1980) conjugate prior, with p=1.0, q=r=s=0.1
def sample_gamma_posterior (data, shape0, scale0):
    S = numpy.sum(data)
    Slog = numpy.sum(map(numpy.log, data))
    N = len(data)

    logP = 1.0+Slog
    q = 1.0 + S
    r = 1.0 + N
    s = 1.0 + N
#    lnf_shape = lambda al: (al-1)*logP + special.gammaln(s*al+1) - (1-al*s)*numpy.log(q) -r*special.gammaln(al)

    lnf_shape = lambda al: (al-1)*logP - q/scale0 - r*special.gammaln(al) - (al*s)*numpy.log(scale0)
    shape_new = sampling.slice (lnf_shape, shape0, thin=5, N=1, lower=0, upper=numpy.inf)[0]
#    print "Sfoo", shape0, lnf_shape(shape0)
#    print "....", shape_new, lnf_shape(shape_new)

#    lnf_rate = lambda beta: -beta*q + s*shape_new*numpy.log(beta)
#    rate_new = sampling.slice (lnf_rate, 1.0/scale0, thin=5, N=1, lower=0, upper=numpy.inf)[0]
    lnf_scale = lambda beta: (shape_new-1)*logP - q/beta - r*special.gammaln(shape_new) - (shape_new*s)*numpy.log(beta)
    scale_new = sampling.slice (lnf_scale, scale0, thin=5, N=1, lower=0, upper=numpy.inf)[0]
#    print "Rfoo", 1.0/scale0, lnf_scale(scale0)
#    print "...", scale_new, lnf_scale(scale_new)
    
    return [shape_new, scale_new]

def gamma_transition_kernel (data, shape0, scale0, shape1, scale1):
    S = numpy.sum(data)
    Slog = numpy.sum(map(numpy.log, data))
    N = len(data)

    logP = 1.0+Slog
    q = 1.0 + S
    r = 1.0 + N
    s = 1.0 + N

    lnf_shape = lambda al: (al-1)*logP - q/scale0 - r*special.gammaln(al) - (al*s)*numpy.log(scale0)
    p_shape = lnf_shape(shape1)

    lnf_scale = lambda beta: (shape1-1)*logP - q/beta - r*special.gammaln(shape1) - (shape1*s)*numpy.log(beta)
    p_scale = lnf_scale (scale1)

    return p_shape + p_scale


# proposal

def uniform_proposal (L, U):
    if U < numpy.inf:
        return pwfun.Pwfun ([L, U], [lambda x: 1.0], [lambda x: 0])
    else:
        return pwfun.Pwfun ([L, U], [lambda x: -x], [lambda x: -1])

# bar proposal

NBARS=10
def bar_pair_proposal (arrv, e0, e1):
    dfn = e0.queue().pyDepartureLik (arrv, e0)
    afn = e1.queue().pyArrivalLik (arrv, e1)

    L = e0.a
    U = e1.d
    eps = (U-L)/NBARS
    x = L

    knots = [ L + i*eps for i in xrange(NBARS) ]
    vals = [ dfn(x) + afn(x) for x in knots ]
    fns = [ lambda d: v for v in vals ]
    derivs = [lambda d: 0] * NBARS
    knots.append(U)

    return pwfun.Pwfun (knots, fns, derivs)
        
def bar_final_proposal (arrv, e0):
    return uniform_proposal (e0.a, numpy.inf)
    
# trapezoid proposal
def zoid_pair_proposal (arrv, e0, e1):
    dfn = e0.queue().pyDepartureLik (arrv, e0)
    afn = e1.queue().pyArrivalLik (arrv, e1)
    
#     print "ZOID" 
#     print e0
#     print e1
#     print afn.dump_table()
#     print dfn.dump_table()

    L0, U0 = afn.range()
    L1, U1 = dfn.range()
    L = max(L0,L1)
    U = min(U0,U1)

    eps = (U-L)/NBARS
    x = L

    knots = [ L + i*eps for i in xrange(1,NBARS) ]
    vals = [ dfn(x) + afn(x) for x in knots ]
    for i in xrange(len(vals)):
        if numpy.isnan(vals[i]):
            vals[i] = -numpy.inf
#    print "V", vals

    dxs = [ dfn.fprime(x) + afn.fprime(x) for x in knots ]
    derivs = [lambda d: dx for dx in dxs ]

    fns = []
    for v,x,dx in zip(vals,knots,dxs):
        def genf (v, x, dx):
            return lambda d: v + dx*(d-x)
        f = genf(v,x,dx)
        fns.append (f)
    knots.append(U)

    result = pwfun.Pwfun (knots, fns, derivs)

    #GGG
#     print "----------------\nX    VAL   A+D   RESULT(x)   FNS[i](x)    DXS[i]"
#     for i in xrange(len(knots)-1):
#         x = knots[i]
#         print fns[i], x, fns[i](x)
#         print "%.5f %.5f %.5f %.5f %.5f %.5f" % (x, vals[i], afn(x)+dfn(x), result(x), (fns[i])(x), dxs[i])
        
    return result

def zoid_final_proposal (arrv, e0):
    dfn = e0.queue().pyDepartureLik (arrv, e0)
    L, U = dfn.range()

    eps = 0.1
    x = L+eps

    knots = []
    vals = []
    dxs = []
    last_dx = 0

    while last_dx >= 0:
        knots.append(x)
        vals.append(dfn(x))
        dxs.append(dfn.fprime(x))
        last_dx = dxs[-1]
        x += eps
        eps *= 2

    derivs = [lambda d: dx for dx in dxs ]
    fns = []
    for v,x,dx in zip(vals,knots,dxs):
        def genf (v, x, dx):
            return lambda d: v + dx*(d-x)
        f = genf(v,x,dx)
        fns.append (f)
    knots.append(U)

    result = pwfun.Pwfun (knots, fns, derivs)

    #GGG
#     print "----------------\nX    VAL   A+D   RESULT(x)   FNS[i](x)    DXS[i]"
#     for i in xrange(len(knots)-1):
#         x = knots[i]
#         print fns[i], x, fns[i](x)
#         print "%.5f %.5f %.5f %.5f %.5f %.5f" % (x, vals[i], dfn(x), result(x), (fns[i])(x), dxs[i])
        
    return result

# for unit tests
def check_quantiles (unittest, s1, s2, N):
    # check each 10% quantile
    s1.sort()
    s2.sort()

    bad = None
    
    for quantile in xrange(3,8):
        i1 = (len(s1)/10)*quantile
        i2 = (len(s2)/10)*quantile        
        diff_max = min (0.5, 15 / sqrt(N))

        print "Decile %d ... %.4f %.4f [delta %.4f max %.4f]" % (quantile, s1[i1], s2[i2], abs(s1[i1]-s2[i2]), diff_max)

        if abs(s1[i1] - s2[i2]) > diff_max:
           bad = "Mismatch %s quantile: %.4f (rejection) vs %.4f (ARS) [diff_max %.4f]" % (quantile, s1[i1], s2[i2], diff_max) 
        
    if bad:
        raise AssertionError (bad)


# for unit tests
def compute_hist (lst):
    l2 = lst[:]
    l2.sort()
    idxs = [ i*len(l2)/10 for i in range(10) ]
    idxs.append (len(l2)-1)
    return [ l2[i] for i in idxs ]

def print_hist (lst):
    print "\n".join (str(x) for x in compute_hist(lst))
