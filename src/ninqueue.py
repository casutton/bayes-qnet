import bisect

class NinQueue:
    def __init__ (self):
        self.a = ArrivalProcess()
        self.d = ArrivalProcess()

    def __repr__ (self):        
        ret = ""
        for t,N in self.knots():
            ret += "%.15f %d %d\n" % (t, N, self.N(t))
        return ret

    def knots (self):
        N = 0
        lst = []
        ret = []
        for t in self.a.times: lst.append ((t,1))
        for t in self.d.times: lst.append ((t,-1))
        lst.sort(cmp_knots)
        for t,delta in lst:
            N += delta
            ret.append ((t, N))
        return ret

    def add_birth_death (self, t_a, t_d):
        self.a.add_time (t_a)
        self.d.add_time (t_d)
        
    def move_arrival (self, a0, a1):
        self.a.move_time (a0, a1)
        
    def move_departure (self, a0, a1):
        self.d.move_time (a0, a1)            

    def N (self, t):
        na = self.a.nprev (t)
        nd = self.d.nprev (t)
        return na - nd

    def contains_arrival (self, t):
        return self.a.contains (t)

    def contains_departure (self, t):
        return self.d.contains (t)

    def knots_in_range (self, l, r):
        knots = dict()
        for a in self.a.in_range (l,r):
            knots[a] = 1
        for d in self.d.in_range (l,r):
            knots[d] = 1
        knots[l] = 1
        knots[r] = 1

        times = knots.keys()
        times.sort()        
        return [ (t, self.N(t)) for t in times ]

    def test_knots (self):
        lst = self.knots()
        for t,n in lst:
            assert n >= 0, ("NINQ less than 0\n%s" % lst)
            
# works just like NinQueue
class Overlay:
    def __init__(self, ninq):
        self.inner = ninq
        self.plus = ArrivalProcess()
        self.minus = ArrivalProcess()

    def __repr__(self):
        xs = [ tup[0] for tup in self.inner.knots() ]
        xs.extend (self.plus.times)
        xs.extend (self.minus.times)
        xs.sort()

        foo = ""
        for x in xs: foo += "%.5f %d\n" % (x, self.N(x))

        return "OVERLAY\n%s\n---\n%s\n---\n%s\n---\n%s/////////" % (foo, self.inner, self.plus, self.minus)

    def move_arrival (self, a0, a1):
        self.plus.add_time (a1)
        self.minus.add_time (a0)

    def move_departure (self, a0, a1):
        self.minus.add_time (a1)
        self.plus.add_time (a0)

    def N (self, t):
        delta = self.plus.nprev(t) - self.minus.nprev(t)
        return self.inner.N(t) + delta

    def knots_in_range (self, l, r):
        knots = dict()
        for a in self.inner.a.in_range (l,r):
            knots[a] = 1
        for d in self.inner.d.in_range (l,r):
            v = knots.get(d, 0)
            knots[d] = v - 1        
        for a in self.plus.in_range (l,r):
            v = knots.get(a, 0)
            knots[a] = v + 1        
        for d in self.minus.in_range (l,r):
            v = knots.get(d, 0)
            knots[d] = v - 1          
        if r not in knots: knots[r] = 0
        lst = knots.items()
        lst.sort()
        n0 = self.N(l) - knots.get(l,0)
        return deltas_to_cumsum (lst, l, n0)


class ArrivalProcess:
        
    def __init__ (self):
        self.times = []
        self.t2ix = None
        self.t2ix_stale = True

    def __repr__ (self): return str(self.times) + "\n" + str(self.t2ix)

    def add_time (self, t):
        bisect.insort (self.times, t)
    
    def move_time (self, t0, t1):
        # works no matter t0 < t1 or reverse: first delete, then insert
        # first delete
        ix = bisect.bisect (self.times, t0)
        assert abs(self.times[ix-1] - t0) < 1e-13, "Off moving %.10f --> %.10f\n%s" % (t0, t1, self)
        del self.times[ix-1]
        # then insert
        self.add_time (t1)
        
    def contains (self, t):
        ix = bisect.bisect (self.times, t)
        return (ix > 0) and (self.times[ix-1] == t)

    def nprev (self, t):
        return bisect.bisect_right (self.times, t)
        
    def in_range (self, l, r):
        ix0 = bisect.bisect_left (self.times, l)
        ix1 = bisect.bisect_right (self.times, r)
        return self.times[ix0:ix1]

class SortedList:

    def __init__ (self):
        self.k = []
        self.k2e = dict()

    def __repr__ (self):
        return "SORTED LIST\n%s\n%s\n//////////\n" % ("\n".join(map(str, self.k)), self.k2e)

    def add (self, k, evt):
        if k in self.k2e:
            self.k2e[k].append (evt)
        else:
            self.k2e[k] = [evt]
        bisect.insort (self.k, k)

    def in_range (self, l, r):
        i0 = bisect.bisect_left (self.k, l)
        i1 = bisect.bisect_right (self.k, r)
        lst = []
        for k0 in self.k[i0:i1]:
            lst.extend (self.k2e[k0])
        return lst

    def remove (self, k, eid):
        assert isinstance(eid, int), "Bad EID %s" % eid
        ix = bisect.bisect (self.k, k)
        assert self.k[ix-1] == k, "Tried to remove C %.15f from\n%s" % (k, self)
        eixl = [ e.eid for e in self.k2e[k] ]
        if eid in eixl:
            eix = eixl.index (eid)
            del self.k2e[k][eix]
            del self.k[ix-1]
            if len(self.k2e[k]) == 0: del self.k2e[k]
#         for k0 in self.k: 
#             assert k0 in self.k2e, "Mismatch deleting %s (ix was %d)\n %s\n %s" % (k, ix, self.k, self.k2e)

    def contains_key (self, k):
        return k in self.k2e

    def contains_kv (self, k, v):
        return (k in self.k2e) and (v in self.k2e[k])

    def __len__ (self): return len(self.k)

class EventSet:
    def __init__ (self):
        self.eids = dict()
        self.evts = []
    def add (self, evt):
        self.eids[evt.eid] = 1
        self.evts.append (evt)
    def contains(self, evt):
        return evt.eid in self.eids
    def items (self): return self.evts



# utilitites

def deltas_to_cumsum (lst, t0, N0):        
    Ncur = N0
    result = [ (t0, N0) ]
    for t,delta in lst:
        Ncur += delta
        result.append ( (t,Ncur) )
    return result

def cmp_knots (k0, k1):
    return cmp(k0[0], k1[0]) or - cmp(k0[1], k1[1])
