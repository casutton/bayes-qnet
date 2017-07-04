import unittest
import qnet
import qnetu
import numpy
import mytime
import netutils
import estimation
import yaml

import pwfun
import sampling
import arrivals
import qstats
import queues
import test_qnet

from scipy import integrate
from numpy import random

import gc
import arrivals

import sys
from math import sqrt

class TestDelayStation (unittest.TestCase):  

    def test_delay_read (self):
        self.assertEquals (2, self.del1q.num_queues())

        q1 = self.del1q.queue_by_name ("I0")
        q2 = self.del1q.queue_by_name ("Q")
        self.assertTrue (q1 != None)
        self.assertTrue (q2 != None)

        self.assertTrue (isinstance (q1, queues.QueueGGk))
        self.assertTrue (isinstance (q2, queues.DelayStation))

    def test_sample (self):
        N = 100
        arrv = self.del1q.sample (N)

        mu_task = qstats.mean_task_time (arrv)
        mus = qstats.mean_service (arrv)
        mur = qstats.mean_response_time (arrv)

        self.assertAlmostEquals (mus[1], mur[1], 10)
        self.assertTrue (abs(mur[1] - 2.0) * sqrt(N) * 0.5 < 1.96, "MUR looks wrrong, expected 2.0 was %.5f" % mur[1])

    def test_arrival_deriv (self):
        N = 10
        arrv = self.del1q.sample(N)
        task5 = arrv.events_of_task (5)
        for e in task5: e.obs_a = e.obs_d = 0
        for e in task5: print e

        e = task5[1]
        fn = e.queue().pyArrivalLik (arrv, e)

        L,U = fn.range()   # N.B. not infinite, b/c arrival
        xs = [ L + i*(U-L)/10 for i in range(10) ]

        check_derivative (self, fn, fn.fprime, xs)


    def test_departure_deriv (self):
        N = 10
        arrv = self.del1q.sample(N)
        task5 = arrv.events_of_task (5)
        for e in task5: e.obs_a = e.obs_d = 0

        e = task5[1]
        fn = e.queue().pyDepartureLik (arrv, e)

        L,U = fn.range()
        U = L + 5.0
        xs = [ L + i*(U-L)/10 for i in range(10) ]

        check_derivative (self, fn, fn.fprime, xs)


    def test_delay_sem (self):
        nreps = 5
        N = 500
        pct = 0.5
        net = self.del1q

        theta_tot = numpy.zeros(len(net.parameters))

        for rep in range(nreps):
            p_orig = net.parameters[:]
            arrv = net.sample (N)
            obs = arrv.subset_by_task (pct)
            a0 = net.gibbs_initialize (obs)
            estimation.sem (net, a0, 0, 100)
            print "PARAMS_FINAL ", net.parameters
            theta_tot += net.parameters

        theta_avg = theta_tot / nreps
        print "THETA_AVG (GOLD, EST):"
        for p_gold, p_hat in zip(net.parameters, list(theta_avg)):
            print p_gold, p_hat
        for p_gold, p_hat in zip(net.parameters, list(theta_avg)):
            self.assertTrue (abs(p_hat - p_gold) < 0.25)

    def setUp (self):
        self.del1q = qnetu.qnet_from_text (del1q_text)
    
EPS = 1e-10
def check_derivative (testcase, fx, dfx, xs):
    for x in xs:
        df_analytic = dfx(x)
        df_numerical = (fx(x+EPS) - fx(x)) / EPS
        testcase.assertTrue (abs (df_analytic - df_numerical) < 1e-5, "Derivative mismatch: X %.5f  DFX: %.5f  f(X): %.15f  f(X+eps) %.15f  df_numerical %.15f" % (x, df_analytic, fx(x), fx(x+EPS), df_numerical))
    
del1q_text = """
states:
 - name: I0
   queues: [I0]
   successors: [S]
 - name: S
   queues: [Q]
queues:
  - { name: I0, service: [M, 1.0] }
  - { name: Q, service: [M, 2.0], type: DELAY }
"""


def main():
    if len(sys.argv) > 1:
        for test_name in sys.argv[1:]:
            suite = unittest.TestLoader().loadTestsFromName("test_delay_station.TestDelayStation.%s" % (test_name,))
            unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()

