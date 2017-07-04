#!/usr/bin/python

import distributions
import misc
import mytime

ln = distributions.LogNormal (0, 1)

N = 10000000
tmr = mytime.timeit()
allint = []

sample = ln.sample(N)
for i in xrange(N):
    ln.lpdf (sample[i])
    
elapsed = tmr.total("Time for %d evaluations" % N)
print "Evaluations per second = %.4f" % (N / elapsed)
