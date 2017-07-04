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

class TestArrivals (unittest.TestCase):  

    def test_subset_by_time (self):
        net = self.mm1
        sampling.set_seed (2324)
        arrv = net.sample (10)
        maxTime = 50

        def task_in_time (arrv, tid):
            evtl = arrv.events_of_task (tid)
            return min(e.d for e in evtl) < maxTime

        subset = arrv.subset_by_task_fn (task_in_time)

        for e in subset.initial_events():
            self.assertTrue (e.d < maxTime)
            old_e = arrv.event (e.eid)
            self.assertTrue (old_e is not None)
            self.assertTrue (e.a == old_e.a)
            self.assertTrue (e.d == old_e.d)
        for e in arrv.initial_events():
            if e.d < maxTime:
                self.assertTrue (subset.event(e.eid) is not None)
        

    mm1_text = """
states:
  - name: INITIAL
    queues: [ INITIAL ]
    successors: [ TIER1 ]
    initial: TRUE 
  - name: TIER1
    queues: [ WEB1, WEB2 ]
queues:
  - { name: INITIAL, service: [M, 10.0]  }
  - { name: WEB1, service: [M, 3.0] }
  - { name: WEB2, service: [M, 3.0] }
"""


    def setUp (self):
        self.mm1 = qnetu.qnet_from_text (TestArrivals.mm1_text)

def main():
    if len(sys.argv) > 1:
        for test_name in sys.argv[1:]:
            suite = unittest.TestLoader().loadTestsFromName("test_arrivals.TestArrivals.%s" % (test_name,))
            unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()

