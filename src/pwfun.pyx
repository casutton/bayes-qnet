
from scipy.integrate import quad
from numpy.random import rand
import netutils
import numpy 

from misc cimport _integrate

cdef extern from "math.h":
    double exp(double x)
    double log(double x)

# Piecewise  function
cdef class Pwfun:

    def __init__ (self, knots, fns, derivs=[]):
        """Knots: double[] of abscissae.  knots[0] : lower bound of fn, knots[-1] : upper bound.
           Fns : List of functions on each interval.  length(fns) must == length(knots) - 1
           knots: List of abscissae of function.  If knots = [ x0, x1, x2], then fns[0] is used in [x0,x1), fns[1] is [x1,x2), etc.
           """
        self.xs = knots
        self.fns = fns
        self.derivs = derivs
        self.N = len(knots)-1
        self.OOR = -numpy.inf
        assert (len(fns) == len(knots)-1)
        
    # Returns the list of knots
    def knots (self): return self.xs
    
    def derivatives (self): return self.derivs
    
    cdef int _find_knot (self, double x):    
        return find_gt (self.xs, x)-1
        
    def __call__ (self, x): 
        return self.value (x)
        i = self._find_knot (x)

    cdef double value (self, double x):
        i = self._find_knot (x)
        # x is in interval (knots[i] .. knots[i+1])
        if 0 <= i and i < self.N:
            return self.fns[i](x)
        else: return self.OOR
    
    def __repr__ (self):
        ret = "Pwfun L: %.4f .. U: %.4f\n" % (self.xs[0], self.xs[-1])
        for x_fx_dx in zip(self.xs, self.fns, self.derivs):
            ret = ret + "   %.4f  %s  %s\n" % x_fx_dx            
        return ret
        
    def dump_table (self):
        ret = "Pwfun L: %.17g .. U: %.17g\n" % (self.xs[0], self.xs[-1])
        for x,f,d in zip(self.xs, self.fns, self.derivs):
            ret = ret + ("%-20.17g%-20.17g%-20.17g\n" % (x, f(x), d(x)))
        return ret
        
    cdef double L (self): return self.xs[0]
    cdef double U (self): return self.xs[-1]
    def range (self): return (self.L(), self.U())

    def fprime (self, x): 
        i = self._find_knot (x)
        if 0 <= i and i < len(self.derivs):
            # x is in interval (knots[i-1] .. knots[i])
            return self.derivs[i](x)
        else:
            raise Exception ("Invalid x value: was %.5f range %.5f .. %.5f" % (x, self.L(), self.U()))
                
    def integrate_one_section (self, i):
        """Returns the integral of section i of the function"""
        return _integrate (self.fns[i], self.xs[i], self.xs[i+1])
        
    def integrate_sections (self):
        """Returns list of integrals for all sectinons"""
        return map (self.integrate_one_section, xrange(self.N))
        
    def integrate (self):
        """Returns total integral"""
        return sum(self.integrate_sections())
            
    def integrate_exp_section (self, i):
        fn = netutils.expify(self.fns[i])
        return _integrate(fn, self.xs[i], self.xs[i+1])
        
    def integrate_exp_sections (self):
        return map(self.integrate_exp_section, xrange(self.N))
        
    def validate (self):
        for i in xrange(len(self.xs)-1):
            assert self.xs[i] <= self.xs[i+1], "Knots invalid at %d\n %s" % (i, self.xs)
        return True
        
        
cdef class ProdFn:
    cdef object sub_f1
    cdef object sub_f2
    cdef object sub_df1
    cdef object sub_df2
    
    def __init__ (self, f1, f2, df1, df2):
        self.sub_f1 = f1
        self.sub_f2 = f2
        self.sub_df1 = df1
        self.sub_df2 = df2
        
    def f1 (self, x): 
#        print "PROD_FN %.4f %.4f %.4f" % (x, self.sub_f1(x), self.sub_f2(x))
        return self.sub_f1(x) * self.sub_f2(x)
        
    def df1 (self, x):
        return self.sub_f1(x) * self.sub_df2(x) + self.sub_df1(x) * self.sub_f2(x)
        

cdef class SumFn:
    cdef object sub_f1
    cdef object sub_f2
    cdef object sub_df1
    cdef object sub_df2

    def __init__ (self, f1, f2, df1, df2):
        self.sub_f1 = f1
        self.sub_f2 = f2
        self.sub_df1 = df1
        self.sub_df2 = df2

    def f1 (self, x): 
#        print "PROD_FN %.4f %.4f %.4f" % (x, self.sub_f1(x), self.sub_f2(x))
        return self.sub_f1(x) + self.sub_f2(x)

    def df1 (self, x):
        return self.sub_df1(x) + self.sub_df2(x)


cdef object compute_combined_knots (Pwfun f1, Pwfun f2):
    L = max(f1.xs[0], f2.xs[0])
    U = min(f1.xs[-1], f2.xs[-1])

    k3 = [ L ]
    for x in f1.xs:
        if (L < x) and (x < U):
            k3.append (x)
    for x in f2.xs:
        if (L < x) and (x < U):
            k3.append (x)
    k3.sort ()    
    k3.append (U)

    # uniqify
    x_last = None
    k3u = []
    for x in k3:
        if x_last != x:
            k3u.append (x)
            x_last = x
    return k3u

def multiply (f1, f2): return _multiply (f1, f2)
    
cdef Pwfun _multiply (Pwfun f1, Pwfun f2):
    
#    print "K1 ", f1.xs
#    print "K2 ", f2.xs
#    print "K3U ", k3u
    k3u = compute_combined_knots (f1, f2)
    fns = []
    derivs = []
    for i in xrange(len (k3u)-1):
        x = k3u[i]
        j1 = f1._find_knot (x)
        j2 = f2._find_knot (x)
#        print "KNOT %.4f %d %d" % (x, j1, j2)
        pfn = ProdFn (f1.fns[j1], f2.fns[j2], f1.derivs[j1], f2.derivs[j2])
        fns.append (pfn.f1)
        derivs.append (pfn.df1)
        
    return Pwfun (k3u, fns, derivs)
        

def add (f1, f2): return _add (f1, f2)

cdef Pwfun _add (Pwfun f1, Pwfun f2):

    k3u = compute_combined_knots (f1, f2)

#    print "K1 ", f1.xs
#    print "K2 ", f2.xs
#    print "K3U ", k3u

    fns = []
    derivs = []
    for i in xrange(len (k3u)-1):
        x = k3u[i]
        j1 = f1._find_knot (x)
        j2 = f2._find_knot (x)
#        print "KNOT %.4f %d %d" % (x, j1, j2)
        pfn = SumFn (f1.fns[j1], f2.fns[j2], f1.derivs[j1], f2.derivs[j2])
        fns.append (pfn.f1)
        derivs.append (pfn.df1)

    return Pwfun (k3u, fns, derivs)
        
# doesn't quite work
cdef Pwfun _foo_multiply (Pwfun f1, Pwfun f2):
    k1 = f1.xs[:]
    k2 = f2.xs[:]
    
    fns = []
    derivs = []
    k3 = []
    
    xmax = min(k1[-1], k2[-1])
    
    i = j = 0
    x1 = k1.pop(0)
    x2 = k2.pop(0)
    
    while k1 or k2:
        print "%d %d %.4f %.4f" % (i,j,x1,x2)
        this_fn = ProdFn (f1.fns[i], f2.fns[j], f1.derivs[i], f2.derivs[j])
        fns.append (this_fn.f1)
        derivs.append (this_fn.df1)
        if x1 < x2:
            print "==> inc 1"
            k3.append (x1)
            x1 = k1.pop (0)
            if i+1 < len(f1.fns): i = i + 1
        elif x1 == x2:
            print "==> inc 12"
            k3.append (x1)
            x1 = k1.pop(0)
            x2 = k2.pop(0)
            if i+1 < len(f1.fns): i = i + 1
            if j+1 < len(f2.fns): j = j + 1
        else:
            print "==> inc 2"
            k3.append (x2)
            x2 = k2.pop(0)
            if j+1 < len(f2.fns): j = j + 1
        print "==> %.4f : %.4f %.4f" % (k3[-1],x1,x2)
            
    k3.append (xmax)
    print (k3)
    
    return Pwfun (k3, fns, derivs)
    
# Univariate piecewise linear function.  Contains some code that makes 
#   it easy to represent a pointwise minimum of linear functions.
cdef class Pwlin:
    """Class for maintaining a piecewise linear upper bound to a log-concave function."""

    def choose_tangent_points (self, L, U, fprime):
        if L == -numpy.inf and U == numpy.inf:
            return -10., 10   # arbitrary.  must be better ways to choose
        elif L == -numpy.inf:
            x1 = U - 1
            while abs(fprime(x1)) > 1e50: x1 = x1-1        
            x0 = x1 - 1.
        elif U == numpy.inf:
            x0 = L + 1
            while abs(fprime(x0)) > 1e50: x0 = x0+1        
            x1 = x0 + 1
        else:
            # whew. base case
            x0 = (2*L + U)/3
            x1 = (L + 2*U)/3
        return x0,x1
        
    def __init__ (self, L, U, f, fprime):
        
        # don't use L, U for tangent points in case f is undefined at boundaries
        
        x0,x1 = self.choose_tangent_points (L,U,fprime)
        f_L = f(x0)
        f_U = f(x1)
        df_L = fprime(x0)
        df_U = fprime(x1)
        z = intersection (x0, x1, f_L, f_U, df_L, df_U)
        
        self.f = f
        self.fprime = fprime
        self.x_tan = [ x0, x1 ]
        self.heights = [ f_L, f_U ]
        self.derivs = [ df_L, df_U ]
        self.xs = [ L, z, U ]
        
        if any(numpy.isnan(self.xs)) or any(numpy.isnan(self.derivs)) or any(numpy.isnan(self.heights)):
            raise Exception ("NaN in pwfun: %s" % self.dump())

    def __call__ (self, x):
        if x > self.U():
            raise Exception ("Value %s greater than upper bound %s" % (x, self.U()))
        if x < self.L():
            raise Exception ("Value %s less than lower bound %s" % (x, self.L()))            
        i = find_gt (self.xs, x)-1
        if i >= len(self.heights): i = len(self.heights) - 1 # special handling for x == U
        return self.heights[i] + self.derivs[i]*(x - self.x_tan[i])

    def __repr__ (self):
        return "[PWLIN: %s]" % zip(self.x_tan, self.derivs, self.heights)

    cdef double value (self, double x):
        i = find_gt (self.xs, x)-1
        if i >= len(self.heights): i = len(self.heights) - 1 # special handling for x == U
        return self.heights[i] + self.derivs[i]*(x - self.x_tan[i])
    
    def knots (self): return self.xs 
    def U (self): return self.xs[-1]
    def L (self): return self.xs[0]
    
    def tangent_points (self): return self.x_tan
    def section_derivs (self): return self.derivs
    def section_heights (self): return self.heights
    def section_deriv (self, i): return self.derivs[i]
    
    def dump (self): return "<PWLIN> KNOTS: %s\n TAN POINTS: %s \nDERIVS: %s\nHEIGHTS: %s</PWLIN>" % (self.xs, self.x_tan, self.derivs, self.heights)
        
    def dump_table (self):
        ret = "Pwlin: L %.17g  U: %.17g\n" % (self.L(), self.U())
        for x in self.xs:
            ret += "%-20.17g%-20.17g\n" % (x, self(x))
        ret += "================\n"
        for tup in zip(self.x_tan, self.derivs, self.heights):
            ret += "%-20.17g%-20.17g%-20.17g\n" % tup
        ret += "================\n"
        return ret
        
    # override
    def integrate_one_section (self, i):
        x1 = self.xs[i]
        x2 = self.xs[i+1]
        xmid = self.x_tan[i]
        df = self.derivs[i]
        f1 = self.heights[i]
#        print "x: [ %.4f .. %.4f ] mid: %.4f df: %.4f  f1: %.4f" % (x1,x2,xmid,df,f1)
        return (f1 - df*xmid)*(x2 - x1) + df/2*(x2*x2 - x1*x1)

    def integrate_exp (self, i):
        """Returns integral of exp(f) from knots[i] to knots[i+1]."""
        cdef double x1, x2, df, f1, sign
        x1 = self.xs[i]
        x2 = self.xs[i+1]
        xmid = self.x_tan[i]
        df = self.derivs[i]
        f1 = self.heights[i] 
        
        if abs(df) < 1e-10:
            ret = exp(f1)*(x2-x1)
        else:
            ret = (exp(df*(x2 - xmid) + f1) - exp(df*(x1 - xmid) + f1)) / df
        if numpy.isnan(ret):
            raise Exception ("ERROR in integrate_exp(%i) for %s\n  f1: %s  df %s x1 %s mid %s x2 %s" % (i, self.dump(), f1, df, x1, xmid, x2) )
        return ret

    def integrate_exp_sections (self):
        """Returns list of integrals for all sectinons"""
        return map (self.integrate_exp, xrange(len(self.xs) - 1))
     
    def argmin (self, x0): 
        """Modifies this object so that it represents min(self, f),
            where f is the linear function f(x) = c + df*(x - x0).
           For efficiency, assumes that this Pwlin object is an argmin of tangents to a concave function."""
        self.c_argmin(x0)
    
    cdef c_argmin (self, double x0):
        cdef int i
        cdef double z1, z2, df, c
        
        df = self.fprime(x0)
        c = self.f (x0)
        
        i = find_ge (self.x_tan, x0)
        # Now i : index of smaller tangent point that is still bigger than x0
        
#        print "i = %d  x0: %.4f  x_tan: %s" % (i, x0, self.x_tan)

        # Modify tangent points structures
        self.x_tan.insert (i, x0)
        self.derivs.insert (i, df)
        self.heights.insert (i, c)
        # Now i --> index of new point in tangent points structure

        # Handle the knots.  This is complicated
        #  N.B. self.xs[i] is the knot between x_tan[i-1] and x_tan[i]
        #   self.xs[0] = L , self.xs[-1] = U
        
#        print "K1: %s (idx == %d)" % (self.xs, i)
        
        # delete redundant knot (if needed)
        # the knot between old x_tan[i-1] and x_tan[i] is obsolete, because x0 now interposes.  Delete.
        if 0 < i and i < len(self.x_tan)-1:
            del self.xs[i]

#        print "K2: ", self.xs
            
        # handle knot xs[i]
        if i > 0:
            z1 = intersection (self.x_tan[i-1], x0, self.heights[i-1], c, self.derivs[i-1], df)
            self.xs.insert(i, z1)

#        print "K3: ", self.xs
                
        # handle knot xs[i+1] (between x_tan[i] and x_tan[i+1])
        if i+1 < len(self.x_tan):
            z2 = intersection (x0, self.x_tan[i+1], c, self.heights[i+1], df, self.derivs[i+1])
#            print "z2 %.4f x0 %.4f x_next %.4f c %.4f h[i+1] %.4f df %.4f df[i+1] %.4f" % (z2, x0, self.x_tan[i+1], c, self.heights[i+1], df, self.derivs[i+1])
            self.xs.insert(i+1, z2)

#        print "K4: ", self.xs

    def sample_from_section_exp (self, i):
        cdef double unif, k0, v1, df, x0
        if i >= len(self.x_tan) or len(self.x_tan)+1 != len(self.xs):
            print self.dump()
            raise Exception ("Error sampling from section %d" % i)            
            
        unif = rand()
        k0 = self.xs[i]
        k1 = self.xs[i+1]
        df = self.derivs[i]
        
        sgn = sign(df)
        
        if k1-k0 < 1e-10: return k0
        
        if df > 1e-10:
            x = k1 + log (unif + (1-unif) * exp(df*(k0-k1))) / df
        elif df < -1e-10:
            x = k0 + log (unif*exp(df*(k1-k0)) + (1-unif)) / df
        else:            
            x = unif*k0 + (1-unif)*k1

#        print "i %d x %.4f  unif=%.4f; df=%.4f; k0=%.4f; k1=%.4f; x0=%.4f;" % (i, unif, x, df, k0, k1, x0)
#        print self.dump_table()
        
        assert not numpy.isnan(x), "ERROR: Integrated value is NaN\n  i %d   unif=%.4f; df=%s; k0=%s; k1=%s;" % (i, unif, df, k0, k1)
        if k0 > x:
             raise Exception ("ERROR: Sampled value %s (should be in %s ... %s)\n  i %d x %.4f  unif=%.4f; df=%.4f; k0=%.4f; k1=%.4f; x0=%.4f;" % (x, k0, k1, i, unif, x, df, k0, k1, x0))
        if x > k1:
             raise Exception ("ERROR: Sampled value %s (should be in %s ... %s)\n  i %d x %.4f  unif=%.4f; df=%.4f; k0=%.4f; k1=%.4f; x0=%.4f;" % (x, k0, k1, i, unif, x, df, k0, k1, x0))

        return x
        

        
    def integrate_sections (self):
        """Returns list of integrals for all sectinons"""
        return map (self.integrate_one_section, xrange(len(self.xs)-1))

    def integrate (self):
        """Returns total integral"""
        return sum(self.integrate_sections())

    
        
    
cdef double intersection (double x_i, double x_i2, double f_i, double f_i2, double df_i, double df_i2):
    cdef double ret
    if abs(df_i - df_i2) < 1e-10:
        ret = 0.5*x_i + 0.5*x_i2 # arbitrary
    else:
        ret = (f_i2 - f_i + x_i * df_i - x_i2 * df_i2) / (df_i - df_i2)
    return ret      
    
cdef int find_ge (xs, x):
    """Returns the index of the first double in xs that is greather than or equal to x. 
        If X greater than all XS, returns len(xs)"""
    for i from 0 <= i < len(xs):
        if x <= xs[i]:
            return i            
    return len(xs)
    
    
cdef int find_gt (xs, x):
    """Returns the index of the first double in xs that is greather than  x. 
        If X greater than all XS, returns len(xs)"""
    for i from 0 <= i < len(xs):
        if x < xs[i]:
            return i            
    return len(xs)

# not in python 2.4
def any (l):
  for x in l:
    if x: 
      return True
  return False

# perhaps move into utils file
cdef double logminusexp (double x1, double x2):
    if x1 >= x2:
        return x1 + log (1 - exp(x2-x1))
    else:
        raise ValueError

cdef double logsumexp (double x1, double x2):
    if x1 < -50 or x2 < -50:
        return max (x1, x2)
    if x1 < x2:
        return x2 + log (1.0 + exp(x1-x2))
    else:
        return x1 + log (1.0 + exp(x2-x1))

def logsumexp0 (lst):
    m = max(lst)
    if m == -numpy.inf: return -numpy.inf
    return m + numpy.log (numpy.sum(numpy.exp (lst - m)))

cdef double sign (double x):
    if x > 0: return 1
    else: return -1
