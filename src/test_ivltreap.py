import unittest
import ivltreap
import sys
import random

class TestIvltreap (unittest.TestCase):  

   def test_intersect (self):
      keys = dict()
      def adder(n): keys[n.other] = 1
      self.treap1.intersect (4,4, adder)
      self.assertEquals (2, len(keys), "Found: %s" % keys)
      for k in ["A", "E"]:
         self.assertTrue (k in keys)

   def test_random_intersect (self):
      nrep = 100
      ninterval = 100
      for ri in range(nrep):
         treap = ivltreap.IntervalTreap()
         brutef= []
         for ii in range(ninterval):
            s = random.uniform (0, 10*ninterval)
            e = s + random.uniform (0, 10)
            treap.insert (s,e,"")
            brutef.append ((s,e))
         self.assertEquals (ninterval, len(treap))
         treap.validate()
         s0 = random.uniform (0, 10*ninterval)
         e0 = random.uniform (0, 10)
         nfound = [0]
         def inner (node): nfound[0] += 1
         treap.intersect (s0, e0, inner)
         n1 = 0
         for s1,e1 in brutef:
            if s1 < e0 and s0 < e1: n1 += 1
         self.assertEquals (n1, nfound[0])

   def test_random_find (self):
      nrep = 100
      ninterval = 100
      for ri in range(nrep):
         treap = ivltreap.IntervalTreap()
         brutef= []
         for ii in range(ninterval):
            s = random.uniform (0, 10*ninterval)
            e = s + random.uniform (0, 10)
            treap.insert (s,e,"")
            brutef.append ((s,e))
         s0,e0 = random.choice (brutef)
         node = treap.find (s0, e0, "")
         self.assertTrue (node != None)
         self.assertTrue (node.start == s0)
         self.assertTrue (node.end == e0)

   def test_remove (self):
      def adder(n): keys[n.other] = 1

      # first delete irrelevant interval
      keys = dict()
      self.treap1.remove (1,2,"B")
      self.treap1.intersect (4,4, adder)
      self.assertEquals (2, len(keys), "Found: %s" % keys)
      for k in ["A", "E"]:
         self.assertTrue (k in keys)

      # first delete relevant interval
      keys = dict()
      self.treap1.remove (3,5,"A")
      self.treap1.intersect (4,4, adder)
      self.assertEquals (1, len(keys), "Found: %s" % keys)
      self.assertTrue ("E" in keys)


   def test_random_remove (self):
      nrep = 100
      ninterval = 10
      for ri in range(nrep):
         # add random intervals
         treap = ivltreap.IntervalTreap()
         for ii in range(ninterval):
            s = random.uniform (0, 5*ninterval)
            e = s + random.uniform (0, 10)
            treap.insert (s,e,"")
         self.assertEquals (len(treap), ninterval)

         # remove random intervals
         killme = []
         def do_select (foo): 
            r = random.uniform (0.0, 1.0)
            if r < 0.25:
               killme.append(foo)
         treap.traverse (do_select)
         print "#KILL ", len(killme)
         print treap
         Nk = 0
         for node in killme:
            print "KILL [%s,%s)" % (node.start, node.end)
            treap.remove (node.start, node.end, node.other)
            print treap
            Nk += 1
            self.assertEquals (len(treap), ninterval - Nk)

         # compute bf list
         brutef = treap.all_intersect(0, 5*ninterval)
         print "#BRUTEF", len(brutef)
         print brutef[0]
         self.assertEquals (len(brutef) + len(killme), ninterval, 
                            "Huh? L(KILL) %d L(BRUTEF) %d expected sum %d" % (len(killme), len(brutef), ninterval))

         # check that a random intersection agrees
         s0 = random.uniform (0, 10*ninterval)
         e0 = random.uniform (0, 10)
         nfound = [0]
         def inner (node): nfound[0] += 1
         treap.intersect (s0, e0, inner)
         n1 = 0
         for other in brutef:
            if other.start < e0 and s0 < other.end: n1 += 1
         self.assertEquals (n1, nfound[0])

   def test_remove2 (self):
      for i in range(100):
         treap1 = ivltreap.IntervalTreap()
         treap1.insert (1,5,"A")
         treap1.insert (1,2,"B")
         treap1.insert (1,8,"C")
         treap1.remove (1,2,"B")
         self.assertEquals (None, treap1.find (1,2,"B"))
         treap1.validate()

   def setUp (self):
      self.treap1 = ivltreap.IntervalTreap()
      self.treap1.insert (3,5,"A")
      self.treap1.insert (1,2,"B")
      self.treap1.insert (6,8,"C")
      self.treap1.insert (4.5,6.5,"D")
      self.treap1.insert (2,7,"E")

def main():
    if len(sys.argv) > 1:
        for test_name in sys.argv[1:]:
            suite = unittest.TestLoader().loadTestsFromName("test_ivltreap.TestIvltreap.%s" % (test_name,))
            unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()

