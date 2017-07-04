import unittest
import qnet
import qnetu
import queues
import numpy
import numpy.random
import mytime
import netutils
import estimation
import yaml
import arrivals
import qstats
import math

import test_qnet

import sys
import sampling
from math import sqrt

# Test class for stochastic EM on singe-processor queues
#  Main tests compare the M-H sampler to a full slice sampler

class TestSem1 (unittest.TestCase):  

    def test_mm1_sem (self):
        net = self.mm1
        nt = 1000
        pct = 0.05
        sampling.set_seed (324121)
        mu_orig = net.parameters[:]

        queues.set_proposal (0)
        estimation.set_sampler ("SLICE")
        
        arrv = net.sample (nt)
        subset = arrv.subset_by_task (pct)
        initial = net.gibbs_initialize (subset)

        print "GOLD SERVICE (ALL) ", qstats.mean_service (arrv)
        print "GOLD SERVICE (OBS) ", qstats.mean_obs_service (arrv)

        mus, arrvl = estimation.sem (net, initial, 0, 1000, report_fn=estimation.mean_service_fn)
        mu_avg = map(numpy.mean, zip(*mus))

        for mu_star, mu_hat in zip(mu_orig, mu_avg):
            print "%.5f  %.5f" % (mu_star, mu_hat)

        print "SERVICE_OBS", qstats.mean_obs_service (subset)

        for mu_star, mu_hat in zip(mu_orig, mu_avg):
            self.assertTrue (abs (mu_star - mu_hat) < 0.1)
        

    # Tests the slice-only sampler in the same circumstance
    def test_mm1_slice (self):
        net = self.mm1
        nt = 1000
        pct = 0.02
        sampling.set_seed (324121)
        mu_orig = net.parameters[:]

        queues.set_proposal (0)
        
        arrv = net.sample (nt)
        subset = arrv.subset_by_task (pct)
        initial = net.gibbs_initialize (subset)

        resampled = estimation.slice_resample (net, initial, 100, report_fn=estimation.mean_service_fn)
        mu_avg = qstats.mean_service(resampled[-1])

        for mu_star, mu_hat in zip(mu_orig, mu_avg):
            print "%.5f  %.5f" % (mu_star, mu_hat)

        print "SERVICE_OBS", qstats.mean_obs_service (subset)

        for mu_star, mu_hat in zip(mu_orig, mu_avg):
            self.assertTrue (abs (mu_star - mu_hat) < 0.1)

    def test_mm1_slice_stationary (self):
        net = self.mm1
        nt = 100
        nrep = 100
        niter = 1
        sampling.set_seed (324121)
        self.do_test_stationary (net, nt, nrep, niter)

    def test_gamma1_slice_stationary (self):
        net = self.gamma1
        nt = 100
        nrep = 100
        niter = 1
        sampling.set_seed (324121)
        self.do_test_stationary (net, nt, nrep, niter)

    def test_ln1_slice_stationary (self):
        net = self.ln1
        nt = 100
        nrep = 100
        niter = 1
        sampling.set_seed (324121)
        self.do_test_stationary (net, nt, nrep, niter)

    def do_test_stationary (self, net, nt, nrep, niter):
    
        mu_tot = numpy.zeros(len (net.parameters))
        mu_orig = net.parameters[:]

        for ri in range(nrep):
            arrv = net.sample (nt)
            subset = arrv.subset_by_task (0.0, adapt_fn=test_qnet.copy_evt)

            resampled = net.slice_resample (subset, niter)
            net.estimate (resampled)

#            resampled = estimation.slice_resample (net, subset, niter, report_fn=estimation.gen_graphing_report_fn(1.0, "stat%d_" % ri))

            mu_avg = net.parameters[:]
            print "SUBSET\n", subset
            print "RESAMPLED\n", resampled[-1]
            print "THETA ", mu_avg
            mu_tot += numpy.array(mu_avg)

            net.parameters = mu_orig[:]

        mu_tot = mu_tot / nrep

        for mu_star, mu_hat in zip(net.parameters, mu_tot):
            print "%.5f  %.5f" % (mu_star, mu_hat)

        for mu_star, mu_hat in zip(net.parameters, mu_tot):
            self.assertTrue (abs (mu_star - mu_hat) < 0.1)
        

    def setUp (self):
        self.mm1 = qnetu.qnet_from_text (mm1_text)
        self.ln1 = qnetu.qnet_from_text (ln1_text)
        self.gamma1 = qnetu.qnet_from_text (gamma1_text)
        self.oneq = qnetu.qnet_from_text (oneq_text)

mm1_text = """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: TIER1
  queues: [ TIER1_0 , TIER1_1 ]
  successors: [TIER2]
- name: TIER2
  queues: [ TIER2_0 , TIER2_1 ]
  successors: [TIER3]
- name: TIER3
  queues: [ TIER3_0 , TIER3_1 ]
queues:
- { name: INITIAL, service: [M, 1.0 ] }
- { name: TIER1_0, service: [M, 1.3333 ] }
- { name: TIER1_1, service: [M, 1.3333 ] }
- { name: TIER2_0, service: [M, 1.3333 ] }
- { name: TIER2_1, service: [M, 1.3333 ] }
- { name: TIER3_0, service: [M, 1.3333 ] }
- { name: TIER3_1, service: [M, 1.3333 ] }
"""

oneq_text = """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: TIER1
  queues: [ TIER1_0 ]
queues:
- { name: INITIAL, service: [M, 1.5 ] }
- { name: TIER1_0, service: [M, 1.0 ] }
"""

gamma1_text = """
states:
- name: INITIAL
  queues: [INITIAL]
  initial: TRUE
queues:
- { name: INITIAL, service: [G, 2, 1.5 ] }
"""

ln1_text = """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: TIER1
  queues: [ TIER1_0 , TIER1_1 ]
  successors: [TIER2]
- name: TIER2
  queues: [ TIER2_0 , TIER2_1 ]
  successors: [TIER3]
- name: TIER3
  queues: [ TIER3_0 , TIER3_1 ]
queues:
- { name: INITIAL, service: [G, 1.33333, 1.5] }
- { name: TIER1_0, service: [LN, 0.9, 0.5 ] }
- { name: TIER1_1, service: [LN, 0.9, 0.5 ] }
- { name: TIER2_0, service: [LN, 0.9, 0.5 ] }
- { name: TIER2_1, service: [LN, 0.9, 0.5 ] }
- { name: TIER3_0, service: [LN, 0.9, 0.5 ] }
- { name: TIER3_1, service: [LN, 0.9, 0.5 ] }
"""    

def main():
    if len(sys.argv) > 1:
        test_name = sys.argv[1]
        suite = unittest.TestLoader().loadTestsFromName("test_sem1.TestSem1.%s" % (test_name,))
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()

