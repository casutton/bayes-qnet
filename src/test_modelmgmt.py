import unittest
import qnetu
import modelmgmt
import sys

import numpy

class TestModelMgmt (unittest.TestCase):  

   def test_split_state (self):
      net1c, converter = modelmgmt.split_state (self.net1, "S1", "COPY0", "COPY1", [ "Q1" ], ["Q2", "Q3"])
      print net1c
      self.assertEquals (5, net1c.num_states())
      self.assertEquals (6, net1c.num_queues())

      print net1c
      print net1c.fsm.a
      print net1c.fsm.o

      nt = 1200
      arrv = net1c.sample (nt)
      self.assertEquals (arrv.num_events(), nt*3)

      nbyq,nbys = self.queue_state_counts (net1c, arrv)
      s_expected = [1200, 400, 600, 600, 800]
      q_expected = [1200, 400, 400, 400, 600, 600]
      self.assertArraysAlmostEquals (nbys, s_expected)
      self.assertArraysAlmostEquals (nbyq, q_expected)

      arrv_old = self.net1.sample(1200)
      nbyq,nbys = self.queue_state_counts (self.net1, arrv_old)
      s_expected = [1200, 1200, 600, 600]
      q_expected = [1200, 400, 400, 400, 600, 600]
      self.assertArraysAlmostEquals (nbys, s_expected)
      self.assertArraysAlmostEquals (nbyq, q_expected)

      arrv_new = converter(arrv_old)
      nbyq,nbys = self.queue_state_counts (net1c, arrv_new)
      s_expected = [1200, 1200, 600, 600]
      q_expected = [1200, 400, 400, 400, 600, 600]
      self.assertArraysAlmostEquals (nbys, s_expected)
      self.assertArraysAlmostEquals (nbyq, q_expected)

   def test_add_hidden_queue (self):
      nt = 100
      qtext = "{ name: HIDDEN, service: [M, 2.0], type: GG1R }"
      arrv = self.net1.sample (nt)
      net2,arrv2 = modelmgmt.add_hidden_queue (self.net1, arrv, "S1", "S1_HIDDEN", qtext)

      self.assertTrue (5, net2.num_states())
      self.assertTrue (7, net2.num_queues())
      arrv3 = net2.sample (1000)
      nbyq,nbys = self.queue_state_counts (net2, arrv3)
      s_expected = [1000, 1000, 500, 500, 1000]
      q_expected = [1000,  333, 333, 333, 500, 500, 1000]
      self.assertArraysAlmostEquals (nbys, s_expected)
      self.assertArraysAlmostEquals (nbyq, q_expected)

      # check arrivals object
      print arrv2
      self.assertEquals (4*nt, arrv2.num_events())
      arrv2.validate()

   def queue_state_counts (self, net1c, arrv):
      nbyq = [0] * net1c.num_queues()
      nbys = [0] * net1c.num_states()
      for evt in arrv:
         nbyq[evt.qid] += 1
         nbys[evt.state] += 1
      return nbyq,nbys

   def assertArraysAlmostEquals (self, nbys, s_expected):
      for act,exp in zip(nbys, s_expected):
         self.assertTrue (abs(exp-act) <= 0.1*exp, "Missed! exp %d  act %d" % (exp,act))

   def test_for_model_search (self):
      net0 = self.net1
      arrv0 = net0.sample (15)
      
      sname = "S1"
      qs = net0.queues_of_state (1)
      q0 = qs[0] 
      q1 = qs[1]

      qlist = map(lambda q: q.name, qs)
      qlist.remove(q0.name)
      qlist.remove(q1.name)

      sbad = "_%s_ABNORMAL" % sname
      sgood = "_%s_OK" % sname
      qtext = "{ name: HIDDEN, service: [M, 2.0], type: GG1R }"

      net1, converter = modelmgmt.split_state (net0, sname, sbad, sgood, [q0.name, q1.name], qlist)
      arrv1 = converter(arrv0)
      arrv1.validate()

      net2, arrv2 = modelmgmt.add_hidden_queue (net1, arrv1, sbad, "HIDDEN", qtext)

      print net0.as_yaml()
      print net1.as_yaml()
      print net2.as_yaml()
      print arrv1

      self.assertEquals (6, net2.num_states())

      # check that states have been changed in ARRV1
      q3 = net1.queue_by_name ("Q3")
      si = net1.sid_by_name ("_S1_OK")
      evts3 = arrv1.events_of_queue (q3)
      for e in evts3:
         self.assertEquals (si, e.state)

      # the way I've set things up, if it goes through Q3,
      #  it shouldn't go through hidden
      q3 = net2.queue_by_name ("Q3")
      qhid = net2.queue_by_name ("HIDDEN")
      evts3 = arrv2.events_of_queue (q3)
      for e in evts3:
         next = e.next_by_task()
         self.assertTrue (next.queue() is not qhid)

   def setUp (self):
      self.net1 = qnetu.qnet_from_text (self.net1_text)

   net1_text = """
states:
 - name: I0
   queues: [I0]
   successors: [S1]
 - name: S1
   queues: [Q1, Q2, Q3]
   successors: [S2,S3]
 - name: S2
   queues: [Q20]
 - name: S3
   queues: [Q30]
queues:
  - { name: I0, service: [M, 1.0] }
  - { name: Q1, service: [M, 2.0], type: GG1R }
  - { name: Q2, service: [M, 2.0], type: GG1R }
  - { name: Q3, service: [M, 2.0], type: GG1R }
  - { name: Q20, service: [M, 2.0], type: GG1R }
  - { name: Q30, service: [M, 2.0], type: GG1R }
"""

def main():
    if len(sys.argv) > 1:
        for test_name in sys.argv[1:]:
            suite = unittest.TestLoader().loadTestsFromName("test_modelmgmt.TestModelMgmt.%s" % (test_name,))
            unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()

