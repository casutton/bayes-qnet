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

class TestMixtureTemplate (unittest.TestCase):  

    def test_net_output (self):
        print self.net1.as_yaml()

    def test_mean (self):
        q1 = self.net1.queue_by_name ("Q")        
        self.assertAlmostEqual(2.662180, q1.service.mean(), 5)

    def test_std (self):
        q1 = self.net1.queue_by_name ("Q")        
        self.assertAlmostEqual(2.813840, q1.service.std(), 5)
        
    def setUp (self):
        self.net1 = qnetu.qnet_from_text (TestMixtureTemplate.net1_text)

    net1_text = """
states:
  - name: INITIAL
    queues: [ INITIAL ]
    successors: [ TIER1 ]
  - name: TIER1
    queues: [ Q ]
queues:
  - { name: INITIAL, service: [M, 10.0]  }
  - { name: Q, service: [ MIX, 0.75, [M, 3.0], 0.25, [LN, 0.0, 1.0] ] }
"""


def main():
    if len(sys.argv) > 1:
        for test_name in sys.argv[1:]:
            suite = unittest.TestLoader().loadTestsFromName("test_mixture_template.TestMixtureTemplate.%s" % (test_name,))
            unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()

