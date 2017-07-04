import unittest
import ninqueue
import sys

class TestNinqueue (unittest.TestCase):  

   def test_add (self):
      ninq = ninqueue.NinQueue()
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
      print ninq

   def test_add2 (self):
      ninq = ninqueue.NinQueue()
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

   def test_move_arrival1 (self):
      ninq = ninqueue.NinQueue()
      ninq.add_birth_death (1.0, 2.0)
      ninq.add_birth_death (1.75, 2.25)
      print ninq
      ninq.move_arrival (1.0, 0.5)
      self.assertEquals (ninq.knots(), 
                         zip ([0.5, 1.75, 2.0, 2.25 ],
                              [1, 2, 1, 0]))

   def test_move_arrival2 (self):
      ninq = ninqueue.NinQueue()
      ninq.add_birth_death (1.0, 2.0)
      ninq.add_birth_death (1.75, 2.25)
      ninq.add_birth_death (2.1, 2.15)
      print ninq
      ninq.move_arrival (1.0, 0.5)
      self.assertEquals (ninq.knots(), 
                         zip([0.5, 1.75, 2.0, 2.1, 2.15, 2.25 ],
                             [1, 2, 1, 2, 1, 0]))


   def test_overlay1 (self):
      overlay = ninqueue.Overlay(self.ninq) 
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

   def test_move_arrival2 (self):
      ninq = ninqueue.NinQueue()
      ninq.add_birth_death (1.0, 2.0)
      ninq.add_birth_death (1.5, 3.5)
      print ninq
      ninq.move_arrival (1.5, 2.5)
      print ninq
      self.assertEquals (1, ninq.N(1.0))
      self.assertEquals (1, ninq.N(1.5))
      self.assertEquals (0, ninq.N(2.25))
      self.assertEquals (1, ninq.N(3.0))
      
   def test_zero (self):
      ninq = ninqueue.NinQueue()
      ninq.add_birth_death (1.0, 1.0)
      eps = 1e-10
#      self.assertEquals (1, ninq.N(1.0), ninq)
      self.assertEquals (0, ninq.N(1.0 - eps))
      self.assertEquals (0, ninq.N(1.0 + eps))

   # tests fix to a bug where knots_in_range couldn't handle simultaneous arrivals
   def test_knots_in_range (self):
      ninq0 = ninqueue.NinQueue()
      ninq0.add_birth_death (1.0, 2.0)
      ninq0.add_birth_death (1.25, 1.75)
      ninq0.add_birth_death (1.25, 1.3)

      knots_expected = [(1.0, 1), (1.25, 2), (1.25, 3), (1.3, 2), (1.75, 1), (2.0, 0)]
      self.assertEquals (ninq0.knots(), knots_expected)
      
      knots_nonredundant = [(0,0), (1.0, 1), (1.25, 3), (1.3, 2), (1.75, 1), (2.0, 0), (5.0, 0)]
      self.assertEquals (ninq0.knots_in_range(0, 5.0), knots_nonredundant)

      knots_small = [(1.1,1), (1.25, 3), (1.3, 2), (1.5, 2)]
      self.assertEquals (ninq0.knots_in_range (1.1, 1.5), knots_small)

   def setUp (self):
      self.ninq = ninqueue.NinQueue()
      self.ninq.add_birth_death (1.0, 2.0)
      self.ninq.add_birth_death (1.75, 2.25)
      self.ninq.add_birth_death (1.8, 1.81)
      self.ninq.add_birth_death (0.5, 2.5)
      self.ninq.add_birth_death (0.25, 2.1)

def main():
    if len(sys.argv) > 1:
        for test_name in sys.argv[1:]:
            suite = unittest.TestLoader().loadTestsFromName("test_ninqueue.TestNinqueue.%s" % (test_name,))
            unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()

