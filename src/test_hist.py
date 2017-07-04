#!/usr/bin/python

import sys
import unittest
import qnetu
import sampling
import numpy
import netutils

class TestHistogram (unittest.TestCase):  
    
    def test_1proc_conditional (self):
        sampling.set_seed(2334)
        N = 5000
        gi = 2500
        eps = 0.01
        nevts = 2
        self.do_test_conditional (self.net_small, N, gi, eps, nevts)
    
    def test_2proc_conditional (self):
        sampling.set_seed(2334)
        N = 5000
        gi = 5000
        eps = 0.01
        nevts = 4
        self.do_test_conditional (self.net2, N, gi, eps, nevts, True)

    def do_test_conditional (self, net, N, gi, eps, nevts, print_all=False):
        q = net.queue_by_name ("Q")
        pick_evt = lambda arrv,i: arrv.events_of_queue(q)[i]

        a2 = []
        all_arrv = []
        for i in range(N):
            arrv = net.sample (nevts)
            evt = pick_evt (arrv,1)
            all_arrv.append (arrv)
            a2.append (evt.a)
        a_mean = numpy.mean (a2)

        eid_to_fix = nevts / 2

        # collect arrv close to mean

        arrv_eps = [ arrv for arrv in all_arrv if abs(pick_evt(arrv, eid_to_fix).a - a_mean) < eps ]
        print "Mean of queue event a_%d  = %.4f" % (eid_to_fix, a_mean)
        print "Number of selected histories = ", len(arrv_eps)
        a0_sampled = [ pick_evt (arrv,0).a for arrv in arrv_eps ]
        netutils.print_hist (a0_sampled)

        if print_all:
            i = 0
            for arrv in arrv_eps:
                f = open ("th_arrv_sampled_%d" % i, "w")
                f.write (str (arrv))
                f.close()
                i += 1

        # now run the gibbs sampling
        a0_gibbs = []
        for arrv in arrv_eps:
            evts = arrv.events_of_queue(q)
            for i in range(nevts/2):
                evt0 = evts[i]
                evt0.obs_a = 0
                evt0.previous_by_task().obs_d = 0
            arrv_out = net.gibbs_resample (arrv, 0, gi, return_arrv = False)
            a0_gibbs.append (pick_evt(arrv_out[-1],0).a)
        print "====="
        netutils.print_hist (a0_gibbs)

        if print_all:
            i = 0
            for arrv in arrv_eps:
                f = open ("th_arrv_gibbs_%d" % i, "w")
                f.write (str (arrv))
                f.close()
                i += 1
            
        f = open("test_hist_sampled.txt", "w")
        f.write ("\n".join(map (str, a0_sampled)))
        f.close()
        f = open("test_hist_gibbs.txt", "w")
        f.write ("\n".join(map (str, a0_gibbs)))
        f.close()

    net_small_text = """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: STATE
  queues: [ Q ]
queues:
- { name: INITIAL, service: [M, 0.5] }
- { name: Q, service: [M, 0.5] }
"""
    net2_text = """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: STATE
  queues: [ Q ]
queues:
- { name: INITIAL, service: [M, 0.5] }
- { name: Q, service: [M, 1.0], processors: 2 }
"""

    net3_text = """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: STATE
  queues: [ Q ]
queues:
- { name: INITIAL, service: [M, 0.5] }
- { name: Q, service: [M, 1.5], processors: 3 }
"""

    def setUp (self):
        self.net_small = qnetu.qnet_from_text (TestHistogram.net_small_text)
        self.net2 = qnetu.qnet_from_text (TestHistogram.net2_text)
        self.net3 = qnetu.qnet_from_text (TestHistogram.net3_text)

    
def main():
    if len(sys.argv) > 1:
        for test_name in sys.argv[1:]:
            suite = unittest.TestLoader().loadTestsFromName("test_hist.TestHistogram.%s" % (test_name,))
            unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()
