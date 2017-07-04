## A treap for overlapping-interval search
## from bx-python  http://bx-python.trac.bx.psu.edu/browser/trunk/lib/bx/intervals/operations/quicksect.py

import numpy
import math
import time
import sys
import random

class IntervalTreap( object ):
    def __init__ ( self ):
        self.root = None
    def insert ( self, start, end, payload ):
        if self.root:
            self.root = self.root.insert ( start, end, other=payload )
        else:
            self.root = IntervalNode ( start, end, other=payload )
#        self.report() #GGG
    def intersect( self, start, end, report_func ):
        if self.root:
            self.root.intersect (start, end, report_func)
    def all_intersect (self, start, end):
        result = []
        if self.root:
            self.root.intersect (start, end, lambda node: result.append (node))
        return result

    def traverse (self, fn):
        if self.root:
            self.root.traverse(fn)

    def __len__ (self):
        N = [0]
        def inner (node): N[0] += 1
        if self.root: self.root.traverse (inner)
        return N[0]
    def find (self, start, end, payload, trace=False):
        """Finds interval containing object PAYLOAD in tree.
           (START, END) is any interval that overlaps the one with payload.
           (Including this allows us to avoid a brute-force search."""
#        print "FIND", start, end, payload 
        result = list(); result.append (None)
        def inner (node):
            if trace: print "FIND_INNER ", node, "MATCHES: START", (node.start == start), "END", (node.end == end), "PAYLOAD", (node.other is payload)
            if (node.other is payload) and (node.start == start) and (node.end == end):
                assert result[0] is None, \
                   "Double appearance in treap\n  %s\n  %s\n  %s" % (result[0], node, self)
                result[0] = node       
        if self.root:
#            self.root.traverse (inner) #SLOW
            self.root.intersect (start - 1e-10, end + 1e-10, inner, trace=trace)
        return result[0] 
    def remove (self, start, end, payload):
        node = self.find (start, end, payload)
        if node is None: 
            print "ERROR IN", payload
            self.find (start, end, payload, trace=True)
            print self
            raise KeyError ("Can't remove non-existent interval [%s,%s) %s" % (start, end, payload ) )
        self.root = self.root.remove (node)

        # GGG
        node2 = self.find(start,end,payload)
        assert node2 is None, \
           "Not gone\n%s\n%s\n ==found %s\nTREAP\n%s" % (node, node2, node is node2, self)

    def validate (self):
        def inner (node):
            if node.left:
                assert node.left.priority <= node.priority, \
                   "Left priority mismatch\n  %s\n  %s" % (node, node.left)
                assert node.left.start <= node.start, \
                   "Left key mismatch\n  %s\n  %s" % (node, node.left)
            if node.right:
                assert node.right.priority <= node.priority, \
                   "Right priority mismatch\n  %s\n  %s" % (node, node.right)
                assert node.right.start >= node.start, \
                   "Right key mismatch\n  %s\n  %s" % (node, node.right)
            node2 = self.find (node.start, node.end, node.other)
            assert node2 is node, \
                "Find mismatch\n  %s\n  %s" % (node, node2)
            
        self.traverse (inner)

    def report (self):
        N = list([0])
        def inner1(node): N[0] += 1
        if self.root: self.root.traverse (inner1)
        depths = eval_depths (self.root) if self.root else []
        print "TREAP #Nodes %d  #Leaves %d  Avg depth %.5f  Max depth %d" % (N[0], len(depths), numpy.mean(depths), max(depths))

    def __repr__ (self):
        if self.root: return self.root.display()
        else: return "[IntervalTreap: empty]"

def eval_depths (node, depth=0):
    if node.left is None and node.right is None:
        return list([depth])
    ell = list()
    if node.left: ell.extend (eval_depths (node.left, depth+1))
    if node.right: ell.extend (eval_depths (node.right, depth+1))
    return ell

class IntervalNode( object ):
    def __init__( self, start, end, linenum=0, other=None ):
        # Python lacks the binomial distribution, so we convert a
        # uniform into a binomial because it naturally scales with
        # tree size.  Also, python's uniform is perfect since the
        # upper limit is not inclusive, which gives us undefined here.
#cas        self.priority = math.ceil( (-1.0 / math.log(.5)) * math.log( -1.0 / (random.uniform(0,1) - 1)))
        self.priority = random.uniform(0,1)
        self.start = start
        self.end = end
        self.maxend = self.end
        self.minend = self.end
        self.left = None
        self.right = None
        self.linenum = linenum
        self.other = other
    def insert( self, start, end, linenum=0, other=None ):
        root = self
        if start > self.start:
            # insert to right tree
            if self.right:
                self.right = self.right.insert( start, end, linenum, other )
            else:
                self.right = IntervalNode(start, end, linenum, other )
            # rebalance tree
            if self.priority < self.right.priority:
                root = self.rotateleft()
        else:
            # insert to left tree
            if self.left:
                self.left = self.left.insert( start, end, linenum, other )
            else:
                self.left = IntervalNode(start, end, linenum, other )
            # rebalance tree
            if self.priority < self.left.priority:
                root = self.rotateright()
        root.recompute_minmaxend()
        return root

    def remove (self, node):        
        if self is node:
            return self.removeSelf()
        elif self.start < node.start:
            # go right, young man
            if self.right: self.right = self.right.remove (node)
            self.recompute_minmaxend()
            return self
        elif self.start > node.start:
            if self.left: self.left = self.left.remove (node)
            self.recompute_minmaxend()
            return self
        else: # hack for ==
            if self.right: self.right = self.right.remove (node)
            if self.left: self.left = self.left.remove (node)
            self.recompute_minmaxend()
            return self
            

    def removeSelf (self):
        if not self.left: return self.right
        if not self.right: return self.left
        if self.left.priority < self.right.priority:
            root = self.rotateleft()
            root.left = self.removeSelf()
            return root
        else:
            root = self.rotateright()
            root.right = self.removeSelf()
            return root
        root.recompute_minmaxend()
        return root

    def recompute_minmaxend (self):
        if self.right and self.left: 
            self.maxend = max( self.end, self.right.maxend, self.left.maxend )
            self.minend = min( self.end, self.right.minend, self.left.minend )
        elif self.right: 
            self.maxend = max( self.end, self.right.maxend )
            self.minend = min( self.end, self.right.minend )
        elif self.left:
            self.maxend = max( self.end, self.left.maxend )
            self.minend = min( self.end, self.left.minend )

    def rotateright( self ):
        root = self.left
        self.left = self.left.right
        root.right = self
        self.recompute_minmaxend()
        root.recompute_minmaxend()
        return root
        
    def rotateleft( self ):
        root = self.right
        self.right = self.right.left
        root.left = self
        self.recompute_minmaxend()
        root.recompute_minmaxend()
        return root

    def intersect( self, start, end, report_func, trace=False ):
        # use nonstrict inequalities.  they can only matter if one
        #  of the intervals has length 0.  Then we'll call it an intersect,
        if trace: print "INTERSECT checking %.15f %.15f\n%s\n%s" % (start, end, self, self.other.dump())
        if start <= self.end and end >= self.start: report_func( self )
        if self.left and start <= self.left.maxend:
            self.left.intersect( start, end, report_func )
        if self.right and end >= self.start:
            self.right.intersect( start, end, report_func )

    def traverse( self, func ):
        if self.left: self.left.traverse( func )
        func( self )
        if self.right: self.right.traverse( func )

    def __repr__ (self):
        return "[%s %s] %s (PRI:%s) (CACHE:%.5f,%.5f) \n" % (self.start, self.end, self.other, self.priority, self.minend, self.maxend)

    def display (self):
        result = str(self)
        if self.left:
            result += "  LEFT: [%s %s] %s\n" % (self.left.start, self.left.end, self.left.other)
        if self.right:
            result += "  RIGHT: [%s %s] %s\n" % (self.right.start, self.right.end, self.right.other)
        if self.left: result += self.left.display()
        if self.right: result += self.right.display()
        return result

def main():
    test = None
    intlist = []
    for x in range(20000):
        start = random.randint(0,1000000)
        end = start + random.randint(1, 1000)
        if test: test = test.insert( start, end )
        else: test = IntervalNode( start, end )
        intlist.append( (start, end) )
    starttime = time.clock()
    for x in range(5000):
        start = random.randint(0, 10000000)
        end = start + random.randint(1, 1000)
        result = []
        test.intersect( start, end, lambda x: result.append(x.linenum) )
    print "%f for tree method" % (time.clock() - starttime)
    starttime = time.clock()
    for x in range(5000):
        start = random.randint(0, 10000000)
        end = start + random.randint(1, 1000)
        bad_sect( intlist, start, end)
    print "%f for linear (bad) method" % (time.clock() - starttime)

def test_func( node ):
    print "[%d, %d), %d" % (node.start, node.end, node.maxend)

def bad_sect( lst, int_start, int_end ):
    intersection = []
    for start, end in lst:
        if int_start < end and int_end > start:
            intersection.append( (start, end) )
    return intersection

if __name__ == "__main__":
    main()


