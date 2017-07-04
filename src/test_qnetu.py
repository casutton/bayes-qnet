import unittest, sys
import qnet
import qnetu
import yaml
import sampling 
    
oneq_text = """
states:
  - name: INITIAL
    queues: [ INITIAL ]
    successors: [ STATE1 ]
    initial: TRUE 
  - name: STATE1
    queues: [ THE_Q ]
queues:
  - { name: INITIAL, service: [M, 0.1 ] }
  - { name: THE_Q, service: M }
"""

non_fifo_arrv_text = """
!Arrivals
events:
- !Event { arrival: 0, departure: 0.75, queue: INITIAL, state: 0, task: 123 }
- !Event { arrival: 0, departure: 1.25, queue: INITIAL, state: 0, task: 555 }
- !Event { arrival: 0, departure: 1.4, queue: INITIAL, state: 0, task: 7 }
- !Event { arrival: 0, departure: 2.0, queue: INITIAL, state: 0, task: 14 }
- !Event { arrival: 0.75, departure: 1.05, queue: THE_Q, state: 1, task: 123 }
- !Event { arrival: 1.25, departure: 1.60, queue: THE_Q, state: 1, task: 555 }
- !Event { arrival: 1.4, departure: 1.55, queue: THE_Q, state: 1, task: 7 }
- !Event { arrival: 2.0, departure: 2.3, queue: THE_Q, state: 1, task: 14 }
"""

multif_state = """
0 0 INITIAL 0 1
0 1 THE_Q 1 -1
1 2 INITIAL 0 3
1 3 THE_Q 1 -1
2 4 INITIAL 0 5
2 5 THE_Q 1 -1
"""

multif_a = """
0 0
1 10
2 0
3 20
4 0
5 15
"""

multif_d = """
0 10
1 17
2 20
3 22
4 15
5 19
"""

    
twoq_text = """
states:
  - name: INITIAL
    queues: [ INITIAL ]
    successors: [ S1 ]
    initial: TRUE 
  - name: S1
    queues: [ Q1 ]
    successors: [ S2 ]
  - name: S2
    queues: [ Q2 ]
queues:
  - { name: INITIAL, service: [M, 0.1 ] }
  - { name: Q1, service: M }
  - { name: Q2, service: M }
"""

hidden_multif_state = """
0 0 INITIAL 0 1
0 1 Q1 1 2
0 2 Q2 2 -1
1 3 INITIAL 0 4
1 4 Q1 1  5
1 5 Q2 2 -1
2 6 INITIAL 0 7
2 7 Q1 1 8
2 8 Q2 2 -1
"""

# EVT 4 --> 5 hidden
hidden_multif_a = """
0 0
1 10
2 15
3 0
4 15
6 0
7 30
8 32 
"""

hidden_multif_d = """
0 10
1 15
2 17
3 15
5 23
6 30
7 32
8 35
"""

import StringIO

class TestQnetu (unittest.TestCase):  

    def test_dump_arrv (self):
        sampling.set_seed (343)
        net = self.oneq
        arrv = net.sample (5)
        arrv.validate()
        
        yml = yaml.dump (arrv)
        print yml
        
        arrv2 = qnetu.load_arrivals (net, yml)
        arrv2.validate()
        
        print arrv
        print arrv2
        
        for evt in arrv2:
            evt_old = arrv.event (evt.eid)
            self.assertEquals (evt.a, evt_old.a)
            self.assertEquals (evt.d, evt_old.d)
            self.assertAlmostEqual (evt.s, evt_old.s)
    
#    def test_non_fifo (self):
#        print self.non_fifo_arrv
#        self.assertEquals (6, self.non_fifo_arrv.num_events ())
#        self.non_fifo_arrv.validate()

    def test_multif (self):
        statef = StringIO.StringIO (multif_state)
        af = StringIO.StringIO (multif_a)
        df = StringIO.StringIO (multif_d)

        arrv = qnetu.read_multifile_arrv (self.oneq, statef, af, df)
        print arrv

        self.assertEquals (6, arrv.num_events())
        self.assertEquals (22.0, arrv.event(3).d)
        self.assertEquals (2.0, arrv.event(3).s)
        arrv.validate()

    # Did this ever work?
    def fixme_test_hidden_multif (self):
        statef = StringIO.StringIO (hidden_multif_state)
        af = StringIO.StringIO (hidden_multif_a)
        df = StringIO.StringIO (hidden_multif_d)

        arrv = qnetu.read_multifile_arrv (self.twoq, statef, af, df)
        print arrv
        self.assertEquals (3, arrv.num_tasks())
        self.assertEquals (9, arrv.num_events())

        e4 = arrv.event (4)
        e5 = arrv.event (5)
        self.assertTrue (e4.next_by_task() == e5)
        self.assertTrue (e5.previous_by_task() == e4)

        for evt in arrv:
            prev_byt = evt.previous_by_task()
            print "::::::::::\n  %s\n  %s" % (prev_byt, evt)
            if prev_byt:
                assert prev_byt.next_by_task() == evt

        qnet.gibbs_initialize_via_minilp (self.twoq, arrv)
        arrv.validate()

    # shouldn't technically go here
    def test_as_yaml1 (self):
        text = self.twoq.as_yaml()
        expected = """states:
  - name: INITIAL
    queues: [ INITIAL ]
    successors: [ S1 ]
  - name: S1
    queues: [ Q1 ]
    successors: [ S2 ]
  - name: S2
    queues: [ Q2 ]
queues:
  - { name: INITIAL, type: GGk, processors: 1, service: [ M, 0.1 ] }
  - { name: Q1, type: GGk, processors: 1, service: [ M, 1.0 ] }
  - { name: Q2, type: GGk, processors: 1, service: [ M, 1.0 ] }
"""
        print text
        self.assertEquals (expected, text)


    def setUp (self):
        self.oneq = qnetu.qnet_from_text (oneq_text)
        self.twoq = qnetu.qnet_from_text (twoq_text)
#        self.non_fifo_arrv = qnetu.load_arrivals (self.oneq, non_fifo_arrv_text)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        test_name = sys.argv[1]
        suite = unittest.TestLoader().loadTestsFromName("test_qnetu.TestQnetu.%s" % (test_name,))
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()
