#!/usr/bin/python

from scipy.optimize import *
from scipy.integrate import *
import distributions
import misc
import mytime
import numpy

ln = distributions.LogNormal (0, 1)

N = 100
tmr = mytime.timeit()
allint = []

for i in xrange(N):
    allint.append (misc.mcint (ln.lpdf, 0, 3))
    
elapsed = tmr.total("Time for %d integrations" % N)
print "Integrations per second = %.4f" % (N / elapsed)

print "Integrals: mean %.4f  sd %.4f" % (numpy.mean(allint), numpy.std(allint))
print "From quad:", quad(ln.lpdf, 0, 3)


ln2 = distributions.LogNormal (-3.75, 1)

N = 10000
tmr = mytime.timeit()
for i in xrange(N):
    integral = quad (ln2.lpdf, 0, 3)
    
print "integral was %.10f" % integral[0]

elapsed = tmr.total("Time for %d integrations" % N)
print "Integrations per second = %.4f" % (N / elapsed)
