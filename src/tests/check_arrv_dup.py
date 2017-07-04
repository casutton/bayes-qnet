#!/usr/bin/python

import qnet
import qnetu
import sys
import mytime

f = open (sys.argv[1])
net = qnetu.qnet_from_text(f)

arrv = net.sample(1000)

N = 250
tmr = mytime.timeit()
for i in xrange(N):
    dup = arrv.duplicate()
    
elapsed = tmr.total("Time for %d duplications" % N)
print "Dups per second = %.4f" % (N / elapsed)
print "Event dups per second = %.4f" % (N*arrv.num_events() / elapsed)

