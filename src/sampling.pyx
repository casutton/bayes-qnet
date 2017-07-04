from numpy import random, isinf
import numpy
import traceback
import sys
from pwfun import Pwlin

cimport pwfun
from pwfun cimport Pwfun, Pwlin


cdef extern from "math.h":
    double sqrt(double x)
    double exp(double x)
    double abs(double x)
    enum: INFINITY

DEBUG = 1
USE_ARMS = 0

cdef extern from "randomkit.h":
    ctypedef struct rk_state:
        unsigned long key[624]
        int pos
        int has_gauss
        double gauss

    void rk_seed(unsigned long seed, rk_state *state)
    double rk_double(rk_state *state)
    unsigned long rk_interval(unsigned long max, rk_state *state)

cdef extern from "cdist.h":
    double rk_uniform(rk_state *state, double loc, double scale)
    double rk_exponential(rk_state *state, double scale)

# global statistics

cdef int total_ars_values, total_ars_iter, num_ars_rejects, num_mh_rejects
cdef int slice_calls, slice_evals

def reset_statistics ():
    global total_ars_values, total_ars_iter, num_ars_rejects, num_mh_rejects
    total_ars_values = 0
    total_ars_iter = 0
    num_mh_rejects = 0
    num_ars_rejects = 0
    slice_evals = 0
    slice_calls = 0
    
def statistics ():
    return { 'values_generated': total_ars_values, 'ars_iter': total_ars_iter, 
              'ars_rejects': num_ars_rejects, 'mh_rejects': num_mh_rejects,
              'slice_calls': slice_calls, 'slice_evals': slice_evals }
    
reset_statistics()


# external interface

def arms (fn, fprime, L, U, x_initial=None, n=1, num_points=0):
    cdef int i
    cdef double x_cur
    
    l = []

    if not x_initial: x_initial = choose_initial (L, U)
    
    fu = Pwlin (L, U, fn, fprime)
    for i in xrange(num_points):
        knots = fu.knots()
        i = random.randint (len(knots) - 2)
        fu.argmin (knots[i+1])
        
    x_cur = x_initial
    
    for i from 0 <= i < n:
        x_cur = _arms (fn, fprime, L, U, x_cur, fu)
        l.append (x_cur)
        
    print "ARS done.  Function %s "  % fu.dump()
    
    return l

def arms_pwfun (Pwfun fn, x_initial=None):
    if not x_initial: 
        x_initial = choose_initial (fn.L(), fn.U())

    return _arms_pwfun (fn, x_initial)
        

cdef double _arms_pwfun (Pwfun fn, double x_initial) except *:
    """Uses ARS to sample from a Pwfun."""
    if fn.U() - fn.L() < 1e-10:
        return fn.L()
    
    areas = fn.integrate_exp_sections ()
    s_idx = _roll_die_unnorm (areas)

    knots = fn.knots()
    if s_idx+1 >= len(knots):
        print "ERROR: Sampled value off the end"
        print knots
        print areas
        print s_idx
        print fn.dump_table()
        print "F: %s ... %s" % (fn(knots[0]), fn(knots[-1] - 1e-10))

    L = knots[s_idx]
    U = knots[s_idx+1]
 
    sub_f = fn.fns[s_idx]
    sub_fprime = fn.derivs[s_idx]

    return _arms (fn, fn.fprime, L, U, x_initial, Pwlin(L, U, sub_f, sub_fprime))
    
    
cdef double _arms (fn, fprime, double L, double U, double x_old, Pwlin fu) except *:
    
    global total_ars_values, total_ars_iter, num_ars_rejects, num_mh_rejects
    
    cdef double x, C, unif, df_i, f_i, k_i, u_x, fn_x, ratio
    cdef int i
    cdef int iter, nrej
    iter = 0
    nrej = 0

    # propose x
    i = select_section (fu)
    
    # ARS reject loop
    for iter from 1 <= iter < 100:
            
        x = fu.sample_from_section_exp (i)

        # accept/reject
        unif = random.rand()
        u_x = fu.value(x)        
        fn_x = fn(x)
        ratio = exp(fn_x - u_x)
        
#        if (ratio > 1.0 + 1e-3) and abs((u_x - fn_x) / fn_x) > 1e-14:
#            print "iter %d :: i: %d U: %.25g x: %.25g f(x): %.32g u(x): %.32g diff: %s ratio: %s" % (iter, i, unif, x, fn_x, u_x, (fn_x - u_x), ratio)
#            print fu.dump_table()
#            raise Exception ("ARS: Function not log-convex")
            
        if DEBUG:
            print "iter %d :: i: %d U: %.4f x: %.4f f(x): %.4f u(x): %.4f ratio: %.4f" % (iter, i, unif, x, fn(x), u_x, ratio)
        
        if unif < ratio:
            # accept
            break
        else:
            # reject.  add a knot
            nrej += 1
            num_ars_rejects += 1
            
            if ratio < 0.9:
                fu.argmin (x)
                
            i = select_section (fu)
            cs = fu.integrate_exp_sections()
            i = _roll_die_unnorm (cs)

    if iter == 100:
        print "WARNING: ARMS max iterations reached, x is %.4f\n  L %.17f  U %.17f" % (x, L, U)
        print fu
        
    # Metropolis rejec
    cdef double mh_ratio, fn_old, fu_old
    
    if USE_ARMS:
        # This is broken in the case that we're using arms_pwfun, and
        #   x_old comes from a different region than we're sampling from.
        unif = random.rand()
        fn_old = fn(x_old)
        fu_old = fu.value(x_old)
    
        mh_ratio = ((fn_x + min(fn_old, fu_old)) - (fn_old + min(fn_x, u_x)))
    
        if DEBUG: 
            print "x_old: %.17g  x: %.17g fn_x: %.17g  ratio: %.17g" % (x_old, x, fn_x, mh_ratio)
            print fu.dump()
            print fu.dump_table()
            print fn.dump_table()

        if mh_ratio < 0: #  if not, accep
            if unif > exp(mh_ratio):
                num_mh_rejects += 1
                x = x_old   # MH reject
            else:
                pass        # MH accept
                
    total_ars_iter += iter
    total_ars_values += 1
    
    return x
    

cdef int select_section (fu) except -1:
    cs = fu.integrate_exp_sections()
    
    i = _roll_die_unnorm (cs)

    if len(cs) <= i:
        print "ERROR: Cannot sample function"
        print fu.dump()
        print fu.dump_table()
        raise Exception ("Drew value off end. value %d  probs %s" % (i, cs))
        
    return i
    
    
def roll_die_unnorm (p):
    return _roll_die_unnorm (p)
    
cdef int _roll_die_unnorm (p):
    cdef double Z, u, S
    cdef int i
    
    Z = sum(p)
    if Z < 1e-50:
        return random.randint (0, len(p))
        
    u = Z*random.rand()
    S = 0
    
    for i from 0 <= i < len(p):
        S = S + p[i]
        if u < S: break
        
    return i
    
def random_elt (lst):
    return lst[random.randint (0, len(lst))]
    
def choose_initial (L, U):
    if isinf(U) and isinf(L):
        x_initial = 0
    elif isinf(U):
        x_initial = L + 1
    elif isinf(L):
        x_initial = U - 1
    else:
        x_initial = (U + L) / 2
    return x_initial

    


### Slice sampling


# Transcribed from Radford Neal's R code


def slice (g, x0, N=1, thin=2, w=1, m=numpy.inf, lower=-numpy.inf, upper=+numpy.inf):
    s = numpy.zeros(N)
    cdef double x_cur = float(x0)
    for i from 0 <= i < N:
        for j from 0 <= j < thin:
            x_cur = _slice (x_cur, g, lower, upper)
        s[i] = x_cur
    return s
    
cdef double w0 = 1
cdef unsigned long m = 1024

cdef double _slice (double x0, g, double lower, double upper) except *:

  global slice_calls, slice_evals

  cdef unsigned long J, K
  cdef double logy, L, R, gx0, gx1
  cdef double w, u
  cdef int num_steps = 0
  
  w = w0
  
  # Keep track of the number of calls made to this function.
  slice_calls += 1

  # Find the log density at the initial point
  slice_evals += 1
  gx0 = g(x0)

  # Determine the slice level, in log terms.
  logy = gx0 - rk_exponential (the_state, 1.0)

  # robustness wrt point masses
  if logy == -INFINITY:
      print "WARNING: Point mass in slice(): x0 %.10f  g(x0) %.10f logy %.10f" % (x0, g(x0), logy)
      return x0
  
  # Find the initial interval to sample from.

  if lower > -INFINITY and upper < INFINITY:
      # user gave us initial interval
      L = lower
      R = upper
  else:
      # sample an interval
      u = rk_uniform (the_state, 0, w)
      L = x0 - u
      R = x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff

      # Expand the interval until its ends are outside the slice, or until
      # the limit on steps is reached.

      J = rk_interval (m, the_state)
      K = (m - 1) - J

      while J > 0:
          if L<=lower: break
          slice_evals += 1
          num_steps += 1
          if g(L) <= logy: break
          #        print "J: %d L: %.4f g(L) %.4f logy %.4f" % (J, L, g(L), logy)
          L = L - w
          J = J - 1
          w *= 2
          
      while K > 0:
          if R >= upper: break
          slice_evals += 1
          num_steps += 1
          if g(R) <= logy: break
#        print "K: %d R: %.4f g(R) %.4f logy %.4f" % (K, R, g(R), logy)
          R = R + w
          K = K - 1
          w *= 2

      if num_steps > 25:
          print "Warning: _slice: excessive number of steps: %d\n  L0 %f  R0 %f  L %f  R %f  w %f  m %d" % (num_steps, x0 - u, x0 + (w-u), L, R, w, m)
      
      
      # Shrink interval to lower and upper bounds.

      if L < lower: L = lower
      if R > upper: R = upper

  # Sample from the interval, shrinking it on each rejection.

  cdef int i = 0
  
  while 1:
    x1 = rk_uniform (the_state, L, R-L)
#    x1 = random.uniform(L,R)

    slice_evals += 1
    gx1 = g(x1)

    if DEBUG: print "SLICE %d %d %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f " % (slice_calls, i, x0, x1, gx1, logy, L, R, lower, upper)

    if gx1 >= logy: break

    if x1 > x0: 
        R = x1
    else:
        L = x1

    if R-L < 1e-10:
        print "slice: quitting: interval too small"
        return L
    
    i += 1
    
  return x1

cdef double _slice_pwfun (Pwfun fn, double x_initial) except *:
      """Uses slice sampling to sample from a Pwfun."""
      if fn.U() - fn.L() < 1e-10:
          return fn.L()

      areas = fn.integrate_exp_sections ()
      s_idx = _roll_die_unnorm (areas)

      knots = fn.knots()
      if s_idx+1 >= len(knots):
          print "ERROR: Sampled value off the end"
          print knots
          print areas
          print s_idx
          print fn.dump_table()
          print "F: %s ... %s" % (fn(knots[0]), fn(knots[-1] - 1e-10))
    
      cdef double L, U
      
      L = knots[s_idx]
      U = knots[s_idx+1]

      sub_f = fn.fns[s_idx]
      sub_fprime = fn.derivs[s_idx]

      # suppose i pick a different partition than the one that x was in
      if x_initial < L or U < x_initial:
          x_initial = choose_initial (L, U)
          
      # use thinning 2
      cdef double x
      x = _slice (x_initial, sub_f, L, U)
      x = _slice (x, sub_f, L, U)
      
      return x


# simple rejection sampler

def rejection (Pwfun fn):
    x0 = fn.L()
    x1 = fn.U()
    if numpy.isinf (x0) or numpy.isinf (x1):
        raise Exception ("Can't handle infinite range.")

    # stupid assumption: no mode (works for exponentials)
    M = max(exp(fn(x0)), exp(fn(x1)))
    while True:
        xm = uniform (x0, x1)
        u = rand()
        A = exp(fn(xm)) / M
        if u < A:
            break

    return xm

    
    

# randomkit interface

cdef rk_state internal_state 
cdef rk_state *the_state = &internal_state

cdef rk_state *_state(): return &internal_state

def set_seed (seed):
    rk_seed (seed, &internal_state)
    numpy.random.seed (seed)

set_seed (13)

cdef double rand ():
   return rk_double (the_state)

cdef double uniform (double L, double U):
   return rk_uniform (the_state, L, U-L)

cdef double exponential (double scale):
   return rk_exponential (the_state, scale)
