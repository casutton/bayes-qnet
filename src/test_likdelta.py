import unittest
import qnet
import qnetu
import numpy
import numpy.random
import mytime
import netutils
import estimation
import yaml
import arrivals
import sampling

import qstats

import sys

class TestLikDelta (unittest.TestCase):  

    def test_mm1_delta (self):
        sampling.set_seed(68310)    
        nreps = 10
        ntasks = 10
        pct = 0.5
        net = self.mm1
        self.do_test_delta_internal (net, nreps, ntasks, pct)

    def test_mmk_delta (self):
        sampling.set_seed(68310)    
        nreps = 10
        ntasks = 10
        pct = 0.5
        net = self.mmk
        self.do_test_delta_internal (net, nreps, ntasks, pct)

    def test_mmrss_delta (self):
        sampling.set_seed(68310)    
        nreps = 10
        ntasks = 100
        pct = 0.25
        net = self.mmrss
        self.do_test_delta_internal (net, nreps, ntasks, pct)

    def do_test_delta_internal (self, net, nreps, ntasks, pct):
        for ri in range(nreps):
            arrv = net.sample (ntasks)
            obs = arrv.subset_by_task (pct)
            samples = net.slice_resample (obs, 0, 5)
            arrv_from = samples[len(samples)-1]
            print "Computing LIK0"
            lik0 = net.log_prob (arrv_from)
            for e in arrv_from:
                if not e.obs_d:
#                    print "Testing evt ", e
                    dfn = qnet.GGkGibbs(net, arrv_from, e, lik0).dfn()
                    d0 = e.d
                    d_test = [ d0+delta for delta in [ -0.5, -0.1, 0.1, 0.5, 1.0, 1.5, 3.0 ] ]
                    for d1 in d_test:
#                        print "Testing departure ", d1
                        lik_incremental = dfn(d1)
                        if numpy.isinf (lik_incremental): continue # probably right
                        lik_true = self.compute_full_lik (net, arrv_from, e, d1)
                        print "%d %.4f %.4f %.4f %.4f" % (e.eid, d0, d1, lik_incremental, lik_true)
                        if numpy.isinf(lik_true):
                            self.assertTrue (numpy.isinf(lik_incremental))
                        else:
                            self.assertAlmostEquals (lik_true, lik_incremental, 5)

    def compute_full_lik (self, net, arrv0, evt, d):
        arrv = arrv0.duplicate()
        dl0 = evt.queue().pyDiffListForDeparture (evt, d)
        evt_next = evt.next_by_task()
        if evt_next:
            dl0.extend (evt_next.queue().pyDiffListForArrival(evt_next, d))
        arrv.applyDiffList (dl0, 0)
        return net.log_prob (arrv)

    
    mm1_text = """
states:
  - name: INITIAL
    queues: [ INITIAL ]
    successors: [ TIER1 ]
    initial: TRUE 
  - name: TIER1
    queues: [ WEB1, WEB2 ]
    successors: [ TIER2 ]
  - name: TIER2
    queues: [ APP1 ]
queues:
  - { name: INITIAL, service: [M, 10.0]  }
  - { name: WEB1, service: [M, 3.0] }
  - { name: WEB2, service: [M, 3.0] }
  - { name: APP1, service: [M, 8.0] }
"""

        
    mmk_text = """
states:
  - name: INITIAL
    queues: [ INITIAL ]
    successors: [ TIER1 ]
    initial: TRUE 
  - name: TIER1
    queues: [ WEB1, WEB2 ]
    successors: [ TIER2 ]
  - name: TIER2
    queues: [ APP1 ]
queues:
  - { name: INITIAL, service: [M, 5.0]  }
  - { name: WEB1, service: [M, 3.0], processors: 3 }
  - { name: WEB2, service: [M, 3.0], processors: 4 }
  - { name: APP1, service: [M, 8.0], processors: 2 }
"""
        
    mmrss_text = """
states:
  - name: INITIAL
    queues: [ INITIAL ]
    successors: [ TIER1 ]
    initial: TRUE 
  - name: TIER1
    queues: [ WEB1, WEB2 ]
    successors: [ TIER2 ]
  - name: TIER2
    queues: [ APP1, APP2 ]
    successors: [ TIER3 ]
  - name: TIER3
    queues: [ DB1, DB2 ]    
queues:
  - { name: INITIAL, service: [M, 10.0]  }
  - { name: WEB1, service: [M, 15.0], type: GG1R }
  - { name: WEB2, service: [M, 17.0], type: GG1R }
  - { name: APP1, service: [M, 10.0], type: GG1R }
  - { name: APP2, service: [M, 10.0], type: GG1R }
  - { name: DB1, service: [M, 7.0], type: GG1R }
  - { name: DB2, service: [M, 7.0], type: GG1R }
"""

    def setUp (self):
        self.mm1 = qnetu.qnet_from_text (TestLikDelta.mm1_text)
        self.mmk = qnetu.qnet_from_text (TestLikDelta.mmk_text)
        self.mmrss = qnetu.qnet_from_text (TestLikDelta.mmrss_text)

def main():
    if len(sys.argv) > 1:
        for test_name in sys.argv[1:]:
            suite = unittest.TestLoader().loadTestsFromName("test_likdelta.TestLikDelta.%s" % (test_name,))
            unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()

