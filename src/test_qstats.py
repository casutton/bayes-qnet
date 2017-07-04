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

import qstats

import sys

class TestQstats (unittest.TestCase):  

    def test_utilization1 (self):
        u = qstats.utilization (self.arrv_mm1)
        expected  = [0.61937, 0.10858, 0.32649, 0.72082]
        self.check_lists (u, expected)

    def test_mean_wait (self):
        u = qstats.mean_wait (self.arrv_mm1)
        expected = [8.020012, 0.0, 2.16130, 2.58173]
        self.check_lists (u, expected)

    def test_mean_response (self):
        u = qstats.mean_response_time (self.arrv_mm1)
        expected = [12.185512, 3.65106, 4.90602, 7.4294440]
        self.check_lists (u, expected)

    def test_mean_wait_pct (self):
        u = qstats.mean_wait_percentage (self.arrv_mm1)
        expected = [ 0.57401, 0.0, 0.352484, 0.3654092 ]
        self.check_lists (u, expected)

    def test_mean_qstats_arriving (self):
        u = qstats.mean_qsize_arriving (self.arrv_mm1)
        expected = [ 2.0, 0, 0.75, 0.8 ]
        self.check_lists (u, expected)

    def test_mean_task_time (self):
        u = qstats.mean_task_time (self.arrv_mm1)
        self.assertAlmostEquals (12.084472, u, 5)

    def test_mean_task_time_in_range (self):
        u = qstats.mean_task_time_in_range (self.arrv_mm1, 5, 10.5)
        self.assertAlmostEquals (11.4137, u, 3)
        u = qstats.mean_task_time_in_range (self.arrv_mm1, 11, 25)
        self.assertAlmostEquals (13.09059, u, 3)

    def test_num_tasks_in_range (self):
        print self.arrv_mm1
        self.assertEquals (3, qstats.num_tasks_in_arrival_range (self.arrv_mm1, 5, 10.5))
        self.assertEquals (2, qstats.num_tasks_in_arrival_range (self.arrv_mm1, 11, 25))

    def test_mean_obs_response (self):
        u = qstats.mean_obs_response (self.arrv_mm1)
        print u
        print self.arrv_mm1
        expected = [12.185512, 3.65106, 5.198517, 8.43051]
        self.check_lists (u, expected)

    def test_mean_obs_service (self):
        u = qstats.mean_obs_service (self.arrv_mm1)
        expected = [4.165502, 3.65106, 3.172583, 5.775833]
        self.check_lists (u, expected)

    def test_write_arrv (self):
        qstats.write_arrv (sys.stdout, self.arrv_mm1)


    
    def check_lists (self, u, expected):
       self.assertEquals (len(u), len(expected))
       for a,b in zip(u,expected):
           self.assertAlmostEquals (a, b, 4)
 
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

    mm3_text = """
    states:
      - name: INITIAL
        queues: [ INITIAL ]
        successors: [ TIER1 ]
        initial: TRUE 
      - name: TIER1
        queues: [ WEB1 ]
        successors: [ TIER2 ]
      - name: TIER2
        queues: [ APP1 ]
    queues:
      - { name: INITIAL, service: [M, 10.0]  }
      - { name: WEB1, processors: 3, service: [M, 30.0] }
      - { name: APP1, processors: 4, service: [M, 30.0] }
    """

    mm1_arrvls = """
!Arrivals
events:
- !Event {arrival: 0.0, departure: 5.0573782176061455, obs: 1, queue: INITIAL, state: 0,
  task: 0}
- !Event {arrival: 0.0, departure: 9.5433894717630476, obs: 1, queue: INITIAL, state: 0,
  task: 1}
- !Event {arrival: 0.0, departure: 10.100174658936853, obs: 1, queue: INITIAL, state: 0,
  task: 2}
- !Event {arrival: 0.0, departure: 15.399120173836607, obs: 1, queue: INITIAL, state: 0,
  task: 3}
- !Event {arrival: 0.0, departure: 20.827496579245874, obs: 1, queue: INITIAL, state: 0,
  task: 4}
- !Event {arrival: 5.0573782176061455, departure: 6.1519936034830538, obs: 1, queue: WEB2,
  state: 1, task: 0}
- !Event {arrival: 6.1519936034830538, departure: 12.941934313997622, obs: 1, queue: APP1,
  state: 2, task: 0}
- !Event {arrival: 9.5433894717630476, departure: 16.177979244871636, obs: 1, queue: WEB2,
  state: 1, task: 1}
- !Event {arrival: 10.100174658936853, departure: 17.966515520657275, obs: 1, queue: WEB2,
  state: 1, task: 2}
- !Event {arrival: 15.399120173836607, departure: 19.427648828212721, obs: 0, queue: WEB2,
  state: 1, task: 3}
- !Event {arrival: 16.177979244871636, departure: 22.911136513937258, obs: 0, queue: APP1,
  state: 2, task: 1}
- !Event {arrival: 17.966515520657275, departure: 23.089052981871603, obs: 0, queue: APP1,
  state: 2, task: 2}
- !Event {arrival: 19.427648828212721, departure: 28.781177557554468, obs: 1, queue: APP1,
  state: 2, task: 3}
- !Event {arrival: 20.827496579245874, departure: 24.478559020698714, obs: 1, queue: WEB1,
  state: 1, task: 4}
- !Event {arrival: 24.478559020698714, departure: 33.626615314865091, obs: 1, queue: APP1,
  state: 2, task: 4}

"""
    def setUp (self):
        self.mm1 = qnetu.qnet_from_text (TestQstats.mm1_text)
        self.mm3 = qnetu.qnet_from_text (TestQstats.mm3_text)
        self.arrv_mm1 = qnetu.load_arrivals (self.mm1, TestQstats.mm1_arrvls)

def main():
    if len(sys.argv) > 1:
        for test_name in sys.argv[1:]:
            suite = unittest.TestLoader().loadTestsFromName("test_qstats.TestQstats.%s" % (test_name,))
            unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()

