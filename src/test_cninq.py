import unittest
import cninq
import sys
import numpy

class TestCninq (unittest.TestCase):

   def test_sorted_double (self):
      lst = cninq.SortedDouble (100)
      lst.add_time (30)
      lst.add_time (10)
      lst.add_time (12.5)
      lst.add_time (32.3)
      lst.add_time (3.0)
      lst.add_time (12.5)
      self.assertEquals (3.0, lst.item (0))
      self.assertEquals (32.3, lst.item (5))

   def test_bisect1 (self):
      lst = cninq.SortedDouble (100)
      lst.add_time (1)
      lst.add_time (3)
      lst.add_time (5)
      self.assertEquals (1.0, lst.item (0))
      self.assertEquals (3.0, lst.item (1))
      self.assertEquals (5.0, lst.item (2))
      self.assertEquals (0, lst.num_lte (0.5))
      self.assertEquals (1, lst.num_lte (1.0))
      self.assertEquals (1, lst.num_lte (1.5))
      self.assertEquals (3, lst.num_lte (7.0))

   def test_move_time (self):
      lst = cninq.SortedDouble (100)
      lst.add_time (1)
      lst.add_time (27)
      lst.add_time (3)
      lst.add_time (5)
      self.check_lst (lst, [1,3,5,27])
      lst.move_time (5,2)
      self.check_lst (lst, [1,2,3,27])
      lst.move_time (27, 0.5)
      self.check_lst (lst, [0.5,1,2,3])
      lst.move_time (0.5,7.2)
      self.check_lst (lst, [1,2,3,7.2])
      lst.move_time (2,4)
      self.check_lst (lst, [1,3,4,7.2])

   def test_add (self):
      ninq = cninq.Ninq(100)
      ninq.add_birth_death (1.0, 2.0)
      ninq.add_birth_death (1.75, 2.25)
      ninq.add_birth_death (1.8, 1.81)
      self.assertEquals (ninq.knots(), 
                         zip ([1.0, 1.75, 1.8, 1.81, 2.0, 2.25],
                              [1, 2, 3, 2, 1, 0]))
      ninq.add_birth_death (0.5, 2.5)
      self.assertEquals (ninq.knots(), 
                         zip ([0.5, 1.0, 1.75, 1.8, 1.81, 2.0, 2.25, 2.5],
                              [1, 2, 3, 4, 3, 2, 1, 0]))
      ninq.add_birth_death (0.25, 2.1)
      self.assertEquals (ninq.knots(),
                         zip ([0.25, 0.5, 1.0, 1.75, 1.8, 1.81, 2.0, 2.1, 2.25, 2.5],
                              [1, 2, 3, 4, 5, 4, 3, 2, 1, 0]))

   def test_add2 (self):
      ninq = cninq.Ninq(100)
      ninq.add_birth_death (1.0, 2.0)
      ninq.add_birth_death (1.75, 2.25)
      ninq.add_birth_death (2.1, 2.15)
      self.assertEquals (ninq.knots(),
                         zip([1.0, 1.75, 2.0, 2.1, 2.15, 2.25],
                             [  1,    2,   1,   2,    1,    0]))

   def test_N (self):
      self.assertEquals (0, self.ninq.N(0.1))
      self.assertEquals (2, self.ninq.N(0.75))
      self.assertEquals (2, self.ninq.N(2.1))
      self.assertEquals (0, self.ninq.N(3))



   def test_overlay1 (self):
      overlay = cninq.Overlay(self.ninq) 
      overlay.move_arrival (0.25, 0.1)      
      overlay.move_arrival (0.5, 0.7)      
      overlay.move_departure (1.81, 3.1)
      self.assertEquals (0, self.ninq.N(0.1))
      self.assertEquals (2, self.ninq.N(0.75))
      self.assertEquals (2, self.ninq.N(2.1))
      self.assertEquals (0, self.ninq.N(3))
      self.assertEquals (0, overlay.N(0.01))
      self.assertEquals (1, overlay.N(0.1))
      self.assertEquals (1, overlay.N(0.25))
      self.assertEquals (1, overlay.N(0.5))
      self.assertEquals (2, overlay.N(0.7))
      self.assertEquals (2, overlay.N(0.75))
      self.assertEquals (5, overlay.N(1.81))
      self.assertEquals (1, overlay.N(3.0))

   def test_iterator (self):
      iterator = self.ninq.interval_iterator(0,3)
      expected = [ (0,0.25,0), (0.25, 0.5, 1), (0.5, 1., 2), (1., 1.75, 3), (1.75, 1.8, 4), (1.80, 1.81, 5), (1.81, 2, 4), (2., 2.1, 3), (2.1, 2.25, 2), (2.25, 2.5, 1), (2.5, 3., 0)  ] 
      kts = self.iterator2knots (iterator)
      for t1,t2 in zip (expected, kts):
         self.assertEquals (t1,t2)

   def test_overlay_iterator (self):
      overlay = cninq.Overlay (self.ninq)
      overlay.move_arrival (0.25, 0.33)
      overlay.move_departure (2.5, 2.75)
      iterator = overlay.interval_iterator(0,3)
      
      expected = [ (0,0.33,0), (0.33, 0.5, 1), (0.5, 1., 2), (1., 1.75, 3), (1.75, 1.8, 4), (1.80, 1.81, 5), (1.81, 2, 4), (2., 2.1, 3), (2.1, 2.25, 2), (2.25, 2.75, 1), (2.75, 3., 0)  ]
      
      kts = self.iterator2knots (iterator)
      self.assertEquals (0, kts[0][0])
      self.assertEquals (3, kts[-1][1])

      N_expected = numpy.sum ( (b-a)*N for a,b,N in expected )
      N_actual = numpy.sum ( (b-a)*N for a,b,N in kts )
      self.assertEquals (N_expected, N_actual)

   def test_overlay_iterator2 (self):
      ninq = cninq.Ninq(100)
      ninq.add_birth_death (50.0, 55.0)
      ninq.add_birth_death (52.0, 53.0)
      overlay = cninq.Overlay (ninq)
      overlay.move_departure(55.0, 50.0)
      iterator = overlay.interval_iterator (50.0, 55.0)
      expected = [(50.0, 52.0, 0), (52.0, 53.0, 1), (53.0, 55.0, 0)]
      self.assertEquals (expected, self.iterator2knots (iterator))
                         
   def iterator2knots (self, i):
      ret = []
      while i.has_next():
         ret.append ( (i.T0(), i.T1(), i.N()) )
         i.advance()
      return ret
   
   def test_move_arrival2 (self):
      ninq = cninq.Ninq(10)
      ninq.add_birth_death (1.0, 2.0)
      ninq.add_birth_death (1.5, 3.5)
      print ninq
      ninq.move_arrival (1.5, 2.5)
      print ninq
      self.assertEquals (1, ninq.N(1.0))
      self.assertEquals (1, ninq.N(1.5))
      self.assertEquals (0, ninq.N(2.25))
      self.assertEquals (1, ninq.N(3.0))

   def setUp (self):
      self.ninq = cninq.Ninq(100)
      self.ninq.add_birth_death (1.0, 2.0)
      self.ninq.add_birth_death (1.75, 2.25)
      self.ninq.add_birth_death (1.8, 1.81)
      self.ninq.add_birth_death (0.5, 2.5)
      self.ninq.add_birth_death (0.25, 2.1)

   def check_lst (self, lst, expected):
      self.assertEquals (len(expected), len(lst))
      for i in range(len(expected)):
         self.assertEquals (expected[i], lst.item(i))
   
def main():
    if len(sys.argv) > 1:
        for test_name in sys.argv[1:]:
            suite = unittest.TestLoader().loadTestsFromName("test_cninq.TestCninq.%s" % (test_name,))
            unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()

