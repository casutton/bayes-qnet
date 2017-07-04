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

class TestGG1R (unittest.TestCase):  

   def test_validate (self):
      nt = 100
      N = 10
      for i in range(N):
         arrv = self.twoq.sample (nt)
         arrv.validate()

   def test_order_modify (self):
      sampling.set_seed (23432)
      arrv =self.twoq.sample(2)
      print arrv

      i0 = self.twoq.queue_by_name ("I0")
      q1 = self.twoq.queue_by_name ("Q1")
      q2 = self.twoq.queue_by_name ("Q2")

      evt  = arrivals.Event (90, 99, 0, i0, 0.0, 3.11, 3.11, 0)
      arrv.append(evt)
      evt.queue().recompute_service (evt)

      e_star = arrivals.Event (91, 99, 1, q1, 3.11, 7.0, 10.11, 1)
      arrv.append (e_star)
      e_star.queue().recompute_service (e_star)

      evt =  arrivals.Event (92, 99, 2, q2, 10.11, 0.0, 17.0, 2)
      arrv.append (evt)
      evt.queue().recompute_service (evt)

      evt  = arrivals.Event (93, 100, 0, i0, 0.0, 3.12, 3.12, 0)
      arrv.append(evt)
      evt.queue().recompute_service (evt)

      e2 = arrivals.Event (94, 100, 1, q1, 3.12, 0.0, 10.12, 1)
      arrv.append (e2)
      e2.queue().recompute_service (e2)

      evt =  arrivals.Event (95, 100, 2, q2, 10.12, 0.0, 19.0, 2)
      arrv.append (evt)
      evt.queue().recompute_service (evt)

      for evt in arrv:
         if evt.tid == 99:
            evt.obs_a = evt.obs_d = 0
         else:
            evt.obs_a = evt.obs_d = 1

      arrv.validate()
      print arrv
      self.assertEquals(e2, arrv.final_arrival_by_queue (q1))

      arrvl = self.twoq.slice_resample (arrv, 1)
      print arrvl[-1]

      arrvl[-1].validate()
      self.assertEquals (arrvl[-1].final_arrival_by_queue (q1), e_star)

   def test_departure_random (self):
       N = 100
       arrv = self.twoq.sample (N)
       evts = arrv.events_of_qid (1)[:]
       arrivals.sort_by_arrival(evts)  # now they're sorted by arrival time

       bigger = [0, 0]
       for e0,e1 in zip(evts, evts[1:]):
           if e0.d < e1.d:
               bigger[0] += 1
           else:
               bigger[1] += 1

       arrv.validate()

       # EXPECTED: bigger[0] dist Binom (0.5, N)
       self.assertEquals (bigger[0] + bigger[1], N-1) # sanity
       self.assertTrue (abs(0.5*(N-1) - bigger[1]) < 2*0.5*0.5*(N-1), "Weird # of inversions, expected %d was %d" % (0.5*(N-1), bigger[1])) # 2 SD

   def test_gg1r_initialize (self):
       N = 6
       pct = 0.5
       arrv = self.twoq.sample (N)
       sub = arrv.subset_by_task (pct)
       print "ARRV ", arrv
       print "SUB ", sub
       initial = self.twoq.gibbs_initialize (sub)
       initial.validate()

   def test_multif (self):
       arrv = self.do_read_multif ()
       print arrv
       arrv.validate()
       self.assertEquals (3, arrv.num_tasks())
       self.assertEquals (9, arrv.num_events())

   def test_duplicate (self):
       arrv = self.do_read_multif ()
       arrv.validate()
       dup = arrv.duplicate()
       print arrv
       print dup
       dup.validate()
       self.assertEquals (9, arrv.num_events())
       self.assertEquals (9, dup.num_events())

       e0 = dup.event(0)
       e3 = dup.event(3)
       self.assertTrue (e0 != None)
       self.assertTrue (e3 != None)
       self.assertEquals (e3, e0.next_by_queue())
       self.assertEquals (e0, e3.previous_by_queue())

   def do_read_multif (self):
        statef = StringIO.StringIO (multif_state)
        af = StringIO.StringIO (multif_a)
        df = StringIO.StringIO (multif_d)
        return qnetu.read_multifile_arrv (self.twoq, statef, af, df)

   def test_gg1r_sem (self):
        sampling.set_seed (3121323)

        nreps = 1
        N = 100
        pct = 0.5
        net = self.rho85  #  self.twoq

        theta_tot = numpy.zeros(len(net.parameters))

        for rep in range(nreps):
            p_orig = net.parameters[:]
            arrv = net.sample (N)
            obs = arrv.subset_by_task (pct)

            a0 = net.gibbs_initialize (obs)
            a0.validate()
            print "INITIAL", a0
#            self.assertEquals (3*N, a0.num_events())

            mus,lst = estimation.sem (net, a0, 0, 50)
            print "PARAMS_FINAL ", net.parameters
            theta_tot += net.parameters

            f = open("gg1r_sem_arrv.txt", "w")
            lst[-1].write_csv (f)
            f.close()

            print "RESPONSE_GOLD     ", qstats.mean_response_time(arrv)
            print "RESPONSE_RESAMPLED", qstats.mean_response_time(a0)
            net.parameters = p_orig[:]

            lst[-1].validate()

        theta_avg = theta_tot / nreps
        print "THETA_AVG (GOLD, EST):"
        for p_gold, p_hat in zip(net.parameters, list(theta_avg)):
            print p_gold, p_hat
        for p_gold, p_hat in zip(net.parameters, list(theta_avg)):
            self.assertTrue (abs(p_hat - p_gold) < 0.5)
    
   def test_gibbs_kernel (self):
      sampling.set_seed (23412832)

      net = self.ini_only
      nt = 4
      N = 5000

      arrv_star = net.sample (nt)
      for evt in arrv_star: 
         if evt.next_by_queue():
            evt.obs_a = evt.obs_d = 0
      print arrv_star

      # OK, so only the last task is observed.  True p(v) is gamma shape=2

      allLp = []
      arrv0 = arrv_star.duplicate()
      for i in range(N):
         arrv1 = net.slice_resample (arrv0.duplicate(), 1)[-1]
         print "======"
         print arrv0
         print arrv1 #GGG
         jp = net.log_prob (arrv1)
         kp = net.gibbs_kernel (arrv0, arrv1)
         print "JP %.5f  KP %.5f" % (jp, kp)
         allLp.append (jp - kp)
         net.slice_resample (arrv0, 1)

      print allLp
      gamma = distributions.Gamma (nt, net.parameters[0])
      ievts = arrv_star.initial_events()
      expected = gamma.lpdf (ievts[-1].d)

#      actual = pwfun.logsumexp0 (numpy.array (allLp)) - numpy.log(N)
      actual = numpy.log (numpy.mean ([ numpy.exp(lp) for lp in allLp ]))

      self.assertTrue (abs(expected - actual) < 0.01, "Mismatch: log prob %.10f  R-B %.10f" % (expected, actual))

   def test_ninq (self):
      sampling.set_seed (123421)
      net = self.rho85
      arrv = net.sample (3) #(100)
      q = net.queue_by_name("Q1")
      qid = 1
#      print arrv

      self.do_check_ninq (net, arrv, q, qid)
      for evt in arrv: 
         evt.obs_d = 0
         if evt.qid != 0: evt.obs_a = 0
      arrvl = net.slice_resample (arrv, 5)
#      print arrvl[-1] 
      self.do_check_ninq (net, arrvl[-1], q, qid)

   def do_check_ninq (self, net, arrv, q, qid):
      ninq = arrv.ninq_of_queue (q)
      evts = arrv.events_of_queue (q)
      alla = [ (e.a, 1) for e in evts ]
      alla.extend  ([ (e.d, -1) for e in evts ])
      alla.sort()

      N = 0
      eps = 1e-20
      for t,delta in alla:
         N += delta
         self.assertEquals (N, ninq.N(t), "Time %.5f  Brute force %d  Cached %d" % (t, N, ninq.N(t)))
         self.assertEquals (N, ninq.N(t+eps))

   def test_dfn_initialize(self):
#      qnet.set_initialization_type("DFN1")
#      qnet.set_initialization_type(qnet.gibbs_initialize_via_order)
      qnet.set_initialization_type(qnet.gibbs_initialize_for_gg1r)
      net = self.m2one
      arrv = net.sample(500)
      obs = arrv.subset_by_task(0.5)
      ini = net.gibbs_initialize(obs)
      ini.validate()

      mu1 = qstats.mean_service(arrv)
      mu_ini = qstats.mean_service(ini)
      muo_ini = qstats.mean_obs_service(ini)
      for m1, m2, m3 in zip(mu1, mu_ini, muo_ini):
         print "%.5f %.5f %.5f" % (m1,m2,m3)
      print ini
#      for evt in ini:
#         if not evt.qid == 0:
#            print "QID", evt.qid, "S", evt.s, "WAIT", evt.wait(), evt


   def reporter (self, net, arrv, iter, lp):
      lp_true = net.log_prob (arrv)
      print "SLICE ", iter, qstats.mean_service(arrv)
      self.assertTrue (abs (lp_true - lp) < 0.01, "No match: LP %.5f  cached %.5f" % (lp_true, lp))

   def test_cached_lp (self):
      net = self.m2one
      arrv = net.sample(250)
      obs = arrv.subset_by_task(0.1)
      ini = net.gibbs_initialize(obs)
      ini.validate()
      net.slice_resample (ini, 25, report_fn=self.reporter)

   def test_gamma_bayes (self):
      net = self.gm2one
      arrv = net.sample(100)
      obs = arrv.subset_by_task(0.25)
      ini = net.gibbs_initialize(obs)
      ini.validate()
#      estimation.bayes (net, ini, 100)
      estimation.sem(net, ini, 0, 100)

   def test_invalid1 (self):
      statef = StringIO.StringIO (test_validate_state)
      af = StringIO.StringIO (test_validate_a)
      df = StringIO.StringIO (test_validate_d)
      arrv = qnetu.read_multifile_arrv (self.oneq, statef, af, df)
      print arrv
      try:
         arrv.validate()
         passed = False
      except AssertionError, e:
         # assertion failure expected
         passed = True
      self.assertTrue (passed, "arrv should not have validated:\n%s" % arrv)

   def test_invalid2 (self):
      arrv = qnetu.read_table_from_string (self.oneq, test_validate2_arrv)
      print arrv
      try:
         arrv.validate()
         passed = False
      except AssertionError, e:
         # assertion failure expected
         passed = True
      self.assertTrue (passed, "arrv should not have validated:\n%s" % arrv)

   def test_slice_stays_valid (self):
      sampling.set_seed (234123)
      nt = 100
      niter = 10
      pct = 0.5
      net = self.twoq
      nreps = 10

      for ri in range(nreps):
         gold = net.sample(nt)      
         arrv0 = gold.subset_by_task (pct)
         arrv = net.gibbs_initialize (arrv0)
         print arrv
         print "INITIAL LIK= ", net.log_prob (arrv)

         try:
            arrv.validate()
         except AssertionError:
            print "Initialization didn't validate on rep ", ri
            continue
         
         for ii in range(niter):
            print "Iter ", ii
            arrv.validate()
            mul, arrvl = estimation.sample_departures (net, arrv, 1)
            arrv = arrvl[-1]

   # a simple example that tests the arrival diff list
   def test_arrival_diff_list (self):
        net = self.dfla
        
        statef = StringIO.StringIO (test_dlfa_1_state)
        af = StringIO.StringIO (test_dlfa_1_a)
        df = StringIO.StringIO (test_dlfa_1_d)
        arrv1 = qnetu.read_multifile_arrv (net, statef, af, df)

        statef = StringIO.StringIO (test_dlfa_2_state)
        af = StringIO.StringIO (test_dlfa_2_a)
        df = StringIO.StringIO (test_dlfa_2_d)
        arrv2 = qnetu.read_multifile_arrv (net, statef, af, df)

        print net.queues[0]
        print arrv1
        print arrv2

        lp1 = net.log_prob (arrv1)
        lp2 = net.log_prob (arrv2)        
        print "A1 LP ", lp1
        print "A2 LP ", lp2

        evt = arrv1.event (4)
        lp1b = qnet.apply_departure (arrv1, evt, 1.5, lp1)
        print arrv1
        print "A1 LP b ", lp1b

        self.assertTrue ((lp1b - lp2) < 1e-10)
        
   def setUp (self):
       self.oneq = qnetu.qnet_from_text (oneq_text)
       self.twoq = qnetu.qnet_from_text (twoq_text)
       self.ini_only = qnetu.qnet_from_text (ini_only_text)
       self.rho85 = qnetu.qnet_from_text (rho85_text)
       self.m2one = qnetu.qnet_from_text (m2one_text)
       self.gm2one = qnetu.qnet_from_text (gm2one_text)
       self.dfla = qnetu.qnet_from_text (dlfaq_text)
       
oneq_text = """
states:
 - name: I0
   queues: [I0]
   successors: [S1]
 - name: S1
   queues: [Q0]
queues:
  - { name: I0, service: [M, 1.0], type: GGK, processors: 1 }
  - { name: Q0, service: [M, 2.0], type: GG1R }
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
  - { name: Q1, service: [M, 2.0], type: GG1R }
  - { name: Q2, service: [M, 2.0], type: GG1R }
"""

ini_only_text = """
states:
 - name: I0
   queues: [I0]
queues:
  - { name: I0, service: [M, 1.0] }
"""

rho85_text = """
states:
 - name: I0
   queues: [I0]
   successors: [S1]
 - name: S1
   queues: [Q1]
queues:
  - { name: I0, service: [M, 1.0] }
  - { name: Q1, service: [M, 0.85], type: GG1R }
"""

m2one_text = """
states:
 - name: I0
   queues: [I0]
   successors: [S1]
 - name: S1
   queues: [Q1,Q2,Q3]
   successors: [S2]
 - name: S2
   queues: [Q11]
queues:
  - { name: I0, service: [M, 1.0] }
  - { name: Q1, service: [M, 2.0], type: GG1R }
  - { name: Q2, service: [M, 2.0], type: GG1R }
  - { name: Q3, service: [M, 2.0], type: GG1R }
  - { name: Q11, service: [M, 0.75], type: GG1R }
"""

gm2one_text = """
states:
 - name: I0
   queues: [I0]
   successors: [S1]
 - name: S1
   queues: [Q1,Q2,Q3]
   successors: [S2]
 - name: S2
   queues: [Q11]
queues:
  - { name: I0, service: [M, 1.0] }
  - { name: Q1, service: [G, 4.0, 0.5], type: GG1R }
  - { name: Q2, service: [G, 4.0, 0.5], type: GG1R }
  - { name: Q3, service: [G, 4.0, 0.5], type: GG1R }
  - { name: Q11, service: [M, 0.75], type: GG1R }
"""

multif_state = """
0 0 I0 0 1
0 1 Q1 1 2
0 2 Q2 2 -1
1 3 I0 0 4
1 4 Q1 1 5
1 5 Q2 2 -1
2 6 I0 0 7
2 7 Q1 1 8
2 8 Q2 2 -1
"""

multif_a = """
0 0
1 10
2 13
3 0
4 15
5 19
6 0
7 20
8 22
"""

multif_d = """
0 10
1 13
2 17
3 15
4 19
5 27
6 20
7 22
8 24
"""

# this guy is invalid
test_validate_state = """
0 0 I0 0 1
0 1 Q0 1 -1
1 2 I0 0 3
1 3 Q0 1 -1
2 4 I0 0 5
2 5 Q0 1 -1
"""
test_validate_a = "0 0\n1 1\n2 0\n3 3\n4 0\n5 3.25\n"
test_validate_d = "0 1\n1 2\n2 3\n3 4\n4 3.25\n5 3.5\n"

# this guy is also invalid
test_validate2_arrv = """EID TID QID A    D    S    W   STATE PROC OBS_A OBS_D CMP
0   0   0   0.00 1.00 1.00 0.0 0 0 1 1 NA
1   0   1   1.00 1.50 0.50 0.0 1 0 1 1 NA
2   1   0   0.00 2.50 1.50 1.0 0 0 1 1 NA
3   1   1   2.50 3.25 0.25 0.0 1 0 1 1 NA
4   2   0   0.00 2.75 0.25 2.5 0 0 1 1 NA
5   2   1   2.75 3.00 0.25 0.0 1 0 1 1 NA
"""

# this guy is also invalid
test_validate2_arrv = """EID TID QID A    D    S    W   STATE PROC OBS_A OBS_D CMP
0   0   0   0.00 1.00 1.00 0.0 0 0 1 1 NA
1   0   1   1.00 1.50 0.50 0.0 1 0 1 1 NA
2   1   0   0.00 2.50 1.50 1.0 0 0 1 1 NA
3   1   1   2.50 3.25 0.25 0.0 1 0 1 1 NA
4   2   0   0.00 2.75 0.25 2.5 0 0 1 1 NA
5   2   1   2.75 3.00 0.25 0.0 1 0 1 1 NA
"""


dlfaq_text = """
states:
 - name: I0
   queues: [I0]
   successors: [S1]
 - name: S1
   queues: [Q0]
queues:
  - { name: I0, service: [M, 1.0], type: DELAY }
  - { name: Q0, service: [M, 1.0], type: GG1R }
"""
test_dlfa_1_state = """
0 0 I0 0 1
0 1 Q0 1 -1
1 2 I0 0 3
1 3 Q0 1 -1
2 4 I0 0 5
2 5 Q0 1 -1
"""
test_dlfa_1_a = "0 0.0\n1 1.0\n2 0.0\n3 2.0\n4 0.0\n5 3.0"
test_dlfa_1_d = "0 1.0\n1 2.5\n2 2.0\n3 3.5\n4 3.0\n5 4.5"
test_dlfa_2_state = """
0 0 I0 0 1
0 1 Q0 1 -1
1 2 I0 0 3
1 3 Q0 1 -1
2 4 I0 0 5
2 5 Q0 1 -1
"""
test_dlfa_2_a = "0 0.0\n1 1.0\n2 0.0\n3 2.0\n4 0.0\n5 1.5"
test_dlfa_2_d = "0 1.0\n1 2.5\n2 2.0\n3 3.5\n4 1.5\n5 4.5"

def main():
    if len(sys.argv) > 1:
        for test_name in sys.argv[1:]:
            suite = unittest.TestLoader().loadTestsFromName("test_gg1r.TestGG1R.%s" % (test_name,))
            unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()

