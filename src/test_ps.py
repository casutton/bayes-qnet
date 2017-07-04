import unittest
import qnet
import qnetu
import numpy
import mytime
import netutils
import estimation
import yaml
import StringIO

import pwfun
import distributions
import sampling
import arrivals
import qstats
import queues
import test_qnet

from scipy import integrate
from numpy import random

import arrivals

import sys
from math import sqrt

class TestPS (unittest.TestCase):  

   def test_sample_small (self):
      sampling.set_seed (2341243)
      net = self.twoq
      arrv = net.sample (5)
      print arrv
      arrv.validate()

   def test_sample_validate (self):
      sampling.set_seed (2341243)
      net = self.twoq
      nr = 10
      nt = 100
      for ri in range(nr):
         arrv = net.sample (nt)
         arrv.validate()
         mu = qstats.mean_service (arrv)
         expected = [ 1.0, 0.5, 0.5 ]
         for mu0, exp0 in zip(mu,expected):
            sd = 1 / (exp0 * sqrt(nt))
            self.assertTrue (abs (mu0 - exp0) < 3*sd, "Mismatch (SD: %.5f)\nTRU %s\nEXP %s" % (sd,mu,expected))

   def test_read_multif (self):
      sampling.set_seed (2341243)
      net = self.twoq
      nr = 10
      nt = 100
      for ri in range(nr):
         arrv = net.sample (nt)
         print "ORIG", arrv
         arrv.validate()
         qnetu.write_multif_to_prefix ("ps_test_sample_validate", arrv)
         arrv2 = qnetu.read_multif_of_prefix ("ps_test_sample_validate", net)
#         print "COPY", arrv2
         arrv2.validate()

   def test_read_multif2 (self):
      sampling.set_seed (2341243)
      net = self.twoq
      nr = 10
      nt = 100
      for ri in range(nr):
         arrv = net.sample (nt)
         arrv.validate()
         obs = arrv.subset_by_task (0.5)
#         print "ORIG", arrv
         qnetu.write_multif_to_prefix ("ps_test_sample_validate2", obs)
         arrv2 = qnetu.read_multif_of_prefix ("ps_test_sample_validate2", net)
#         print "COPY", arrv2
         arrv2.validate()

   def test_initialize (self):
      sampling.set_seed (2341243)
      net = self.twoq
      nr = 10
      nt = 100
      for ri in range(nr):
         arrv = net.sample (nt)
         obs = arrv.subset_by_task (0.5)
         ini = net.gibbs_initialize (obs)
         print "TRUE", arrv
         print "INI", ini
         ini.validate()

   def test_initialize_for_ps (self):
      sampling.set_seed (2341243)
      net = self.twoq
      nr = 10
      nt = 100
      for ri in range(nr):
         arrv = net.sample (nt)
         obs = arrv.subset_by_task (0.5)
         ini = qnet.gibbs_initialize_for_ps (net, obs)
         ini.validate()
         
   def test_sem (self):
      sampling.set_seed (67826)
#      net = self.twoq
      net = self.oneq
      nr = 1
      nt = 50
      theta0 = net.parameters[:]
      for ri in range(nr):
         arrv = net.sample (nt)
         obs = arrv.subset_by_task (0.25)
         ini = net.gibbs_initialize (obs)
         estimation.sem (net, ini, 0, 100)
         print "MU  ", net.parameters
         net.parameters = theta0[:]
         print "TRU ", theta0
         
   def test_bayes (self):
      sampling.set_seed (67826)
#      net = self.twoq
      net = self.oneq
      nr = 1
      nt = 100
      theta0 = net.parameters[:]

      def reporter (net, arrv, iter, lp):
         lp_scratch = net.log_prob(arrv)
         assert abs(lp - lp_scratch) < 1e-10, \
            "Mismatch LP. Running total %.10f  from scratch %.10f" % (lp, lp_scratch)
         if 0 == (iter % 10):
            f = open ("ps_test_bayes_%d.txt" % iter, "w")
            arrv.write_csv(f)
            f.close()

      for ri in range(nr):
         arrv = net.sample (nt)
         obs = arrv.subset_by_task (0.25)
         ini = net.gibbs_initialize (obs)
         estimation.bayes (net, ini, 100, report_fn=reporter)
         print "MU  ", net.parameters
         net.parameters = theta0[:]
      print "TRUE ", theta0

   def test_ps_stationary (self):
      nr = 50
      nt = 50
      net = self.twoq
      allmu = numpy.zeros (len(net.parameters))
      allmax = numpy.zeros (len(net.parameters))
      for i in range (nr):
         arrv = net.sample (nt)
         obs = arrv.subset_by_task (0.0)
         net.slice_resample (arrv, 10)
         mu = numpy.array (qstats.mean_service (arrv))
         this_max = numpy.array (qstats.max_service (arrv))
         print "MU", mu
         allmu += mu
         allmax += this_max
      avg = allmu / nr
      print "AVG", avg
      print "TRU", net.parameters
      print "MAX", allmax / nr

   def test_ps_likelihood (self):
      sampling.set_seed (23134)
      net = self.oneq
      net.parameters = [ 10.0, 0.1 ]
      nt = 10
      arrv = net.sample(nt)
      tid = 3
      evts = arrv.events_of_task (tid)
      e1 = evts[1]
      e1.obs_d = 0
      lp0 = net.log_prob (arrv)
      print arrv
      print "LP0", lp0

      dexp = distributions.Exponential (net.parameters[1])
      gibbs = qnet.GGkGibbs (net, arrv, e1, lp0)
      l = e1.a
      u = e1.a + 3.0
      diff = gibbs.inner_dfun(l) - dexp.lpdf(0)
      for i in range(10):
         x = l + 0.1*i*(u-l)
         gval = gibbs.inner_dfun(x)
         print "%.10f  %.10f  %.10f %.10f" % (x, gval, gval - diff, dexp.lpdf(x-l))

   def test_zero_s (self):
      sampling.set_seed (23134)
      net = self.oneq

      arrv = net.sample (1)
      e1 = arrv.event (1)
      q1 = e1.queue()

      print net.parameters
      print arrv
      
      e1.d = e1.a
      e1.s = 0
      lp0 = net.log_prob (arrv)
      self.assertAlmostEquals (-1.47313356106, lp0, 5)

      e1.d = e1.a + 1.
      e1.s = 1.
      dl = q1.pyDiffListForDeparture (e1, e1.a)
      lp1 = net.log_prob (arrv)
      dlik = q1.likelihoodDelta(arrv, dl)
      print arrv
      print lp1, dlik
      self.assertAlmostEquals (-1.47313356106, lp1 + dlik, 5)
      
      
   def test_likelihood_delta (self):
      sampling.set_seed (23134)
      net = self.oneq
      net.parameters = [ 10.0, 5.0 ]
      nt = 10
      tid = 3
      
      arrv = net.sample(nt)
      evts = arrv.events_of_task (tid)
      e1 = evts[1]
      e1.obs_d = 0
      q1 = e1.queue()
      
      lp0 = net.log_prob (arrv)
      print arrv

      d0 = e1.d
      deltas = [-0.005, 0.0, 0.1, 0.5, 1.0]
      for delta in deltas:
         d_new = d0 + delta
         dl = q1.pyDiffListForDeparture (e1, d_new)
         dlik = q1.likelihoodDelta(arrv, dl)
         lik_a = lp0 + dlik

         a2 = arrv.duplicate()
         dl = q1.pyDiffListForDeparture (e1, d_new)
         a2.applyDiffList (dl)
         lik_b = net.log_prob (a2)

         print a2
         print lik_a, lik_b
         
         self.assertTrue (abs(lik_b - lik_a) < 1e-5)
         

   # mainly for profiling atm
   def test_ps_em (self):
      nt = 600
      niter = 2
      net = self.twoq

      arrv = net.sample (nt)
      obs = arrv.subset_by_task (0.5)
      estimation.sem (net, obs, 0, niter)

      print net.parameters

   def setUp (self):
       self.oneq = qnetu.qnet_from_text (oneq_text)
       self.twoq = qnetu.qnet_from_text (twoq_text)

oneq_text = """
states:
 - name: I0
   queues: [I0]
   successors: [S1]
 - name: S1
   queues: [Q0]
queues:
  - { name: I0, service: [M, 1.0] }
  - { name: Q0, service: [M, 0.5], type: PS }
"""

twoq_text = """
states:
 - name: I0
   queues: [I0]
   successors: [S1]
 - name: S1
   queues: [Q1]
   successors: [S2]
 - name: S2
   queues: [Q2]
queues:
  - { name: I0, service: [M, 1.0] }
  - { name: Q1, service: [M, 0.5], type: PS }
  - { name: Q2, service: [M, 0.5], type: PS }
"""

def main():
    if len(sys.argv) > 1:
        for test_name in sys.argv[1:]:
            suite = unittest.TestLoader().loadTestsFromName("test_ps.TestPS.%s" % (test_name,))
            unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()

