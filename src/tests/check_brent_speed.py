#!/usr/bin/python

from scipy.optimize import *
from scipy.integrate import *
import distributions
import mytime

ln = distributions.LogNormal (0, 1)

N = 100000
tmr = mytime.timeit()
for i in xrange(N):
    root = brenth (ln.dx_lpdf, 0.1, 25)
    
print "root was %.10f" % root

elapsed = tmr.total("Time for %d root-findings" % N)
print "Roots per second = %.4f" % (N / elapsed)

N = 10000
tmr = mytime.timeit()
for i in xrange(N):
    the_max = bisect (ln.dx_lpdf, 0.1,25)
    
print "max was %.10f" % the_max

elapsed = tmr.total("Time for %d maxes" % N)
print "Maxes per second = %.4f" % (N / elapsed)


N = 10000
tmr = mytime.timeit()
for i in xrange(N):
    integral = quad (ln.lpdf, 0.1, 3)
    
print "integral was %.10f" % integral[0]

elapsed = tmr.total("Time for %d integrations" % N)
print "Integrations per second = %.4f" % (N / elapsed)

