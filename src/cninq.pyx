# Faster version of ninqueue.py

cdef extern from "math.h":
    enum: INFINITY

cdef extern from "Python.h":
    void *PyMem_Malloc(unsigned int n)
    void PyMem_Free(void *p)

cdef extern from "string.h":
    void *memmove(void *s1, void *s2, size_t n)



cdef class SortedDouble
cdef class OverlayIterator
cdef class Overlay

cdef class Ninq:
 
    def __init__ (self, Nmax):
        self.a = SortedDouble (Nmax)
        self.d = SortedDouble (Nmax)
        self.Nmax = Nmax

    def __repr__ (self):
        ret = "NINQ:\n"
        kts = self.knots()
        for t,N in kts:
            ret += "%.5f %d\n" % (t,N)
        return ret
    
    def knots (self):
        N = 0
        lst = []
        ret = []
        for 0 <= i < len(self.a): lst.append (( self.a.item(i), 1 ))
        for 0 <= i < len(self.d): lst.append (( self.d.item(i), -1 ))
        lst.sort(cmp_knots)
        for t,delta in lst:
            N += delta
            ret.append ((t, N))
        return ret

    cpdef add_birth_death (self, double t_a, double t_d):
        self.a.add_time (t_a)
        self.d.add_time (t_d)
        
    cpdef move_arrival (self, double a0, double a1):
        self.a.move_time (a0, a1)
        
    cpdef move_departure (self, double a0, double a1):
        self.d.move_time (a0, a1)

    cpdef int N (self, double t):
        cdef int na = self.a.num_lte (t)
        cdef int nd = self.d.num_lte (t)
        return na - nd
    
    cpdef OverlayIterator interval_iterator (self, double l, double r):
        cdef Overlay ovl = Overlay(self)
        return ovl.interval_iterator(l, r)

    cpdef int contains_arrival (self, double v):
        cdef int idx = self.a.bisect (v)-1
        return self.a.item (idx) == v

    cpdef int contains_departure (self, double v):
        cdef int idx = self.d.bisect (v)-1
        return self.d.item (idx) == v

# works just like NinQueue
cdef class Overlay:

    def __init__(self, Ninq ninq):
        self.inner = ninq
        self.plus = SortedDouble(ninq.Nmax)
        self.minus = SortedDouble(ninq.Nmax)

    def knots (self):
        N = 0
        lst = []
        ret = []
        for 0 <= i < len(self.inner.a): lst.append (( self.inner.a.item(i), 1 ))
        for 0 <= i < len(self.inner.d): lst.append (( self.inner.d.item(i), -1 ))
        for 0 <= i < len(self.plus): lst.append (( self.plus.item(i), 1 ))
        for 0 <= i < len(self.minus): lst.append (( self.minus.item(i), -1 ))
        lst.sort(cmp_knots)
        for t,delta in lst:
            N += delta
            ret.append ((t, N))
        return ret

    def __repr__(self):
        return "OVERLAY\n"+ "\n".join(map(str, self.knots()))

    cpdef move_arrival (self, double a0, double a1):
        self.plus.add_time (a1)
        self.minus.add_time (a0)

    cpdef move_departure (self, double a0, double a1):
        self.minus.add_time (a1)
        self.plus.add_time (a0)

    cdef OverlayIterator _interval_iterator (self, double l, double r):
        return OverlayIterator (self, l, r)

    cpdef OverlayIterator interval_iterator (self, double l, double r):
        return OverlayIterator (self, l, r)

    cpdef int N (self, double t):
        cdef int n0 = self.inner.N (t)
        cdef int plus0 = self.plus.num_lte (t)
        cdef int minus0 = self.minus.num_lte (t)
        return n0 + plus0 - minus0

cdef class SortedDouble:

    def __init__ (self, cap):
        self.val = <double*> PyMem_Malloc (cap * sizeof(double))
        self.capacity = cap

    def __dealloc__ (self):
        if self.val != NULL:
            PyMem_Free (self.val)
            self.val = NULL

    def __repr__ (self):
        ret = "[SL: "
        for i in range(self.N):
            ret += str(self.val[i])
            ret += ", "
        ret += "]"
        return ret

    def __len__ (self): return self.N

    cpdef double item (self, int i):
        if i < self.N:
            return self.val[i]
        else:
            return INFINITY
    
    cpdef add_time (self, double x):
        cdef int i
        if self.N >= self.capacity:
            raise Exception ("Can't increment SortedDouble past capacity %d" % self.capacity)
        cdef int idx = self.bisect (x)
        for i in range(self.N, idx, -1):
            self.val[i] = self.val[i-1]
        self.N += 1
        self.val[idx] = x

    cpdef move_time (self, double x0, double x1):
        cdef int i
        cdef int idx0 = self.bisect (x0) - 1
        cdef int idx1 = self.bisect (x1)
        if self.val[idx0] != x0: raise Exception ("Can't find %.5f in list\n%s" % (x0, self))
        if idx0 < idx1:
            # need to move stuff backward
            for i in range(idx0, idx1):
                self.val[i] = self.val[i+1]
            self.val[idx1-1] = x1
        elif idx1 < idx0:
            # move other stuff forward; 
            for i from idx0 >= i > idx1:
                self.val[i] = self.val[i-1]
            self.val[idx1] = x1
        else:
            # else equal, don't do memmove
            self.val[idx1] = x1

    cpdef int num_lte (self, v):
        return self.bisect (v)

    # Returns index i so that self.val[i-1] <= v < self.val[i]
    cdef int bisect (self, double v):
        cdef int lo = 0
        cdef int hi = self.N
        cdef int mid
        cdef double midval
        while lo < hi:
            mid = (lo+hi)//2
            midval = self.val[mid]
            if midval < v:
                lo = mid+1
            elif midval > v: 
                hi = mid
            else:
                lo = mid
                break
        # invariant of above: lo <= correct_answer < hi
        # special case if values in self.val equal        
        while (lo < self.N) and (self.val[lo] <= v): lo += 1
        return lo
 

cdef class OverlayIterator:

    def __init__ (self, Overlay ovl, double l, double r):
        self.r = r
        self.t1 = l

        self.a = ovl.inner.a
        self.d = ovl.inner.d
        self.plus = ovl.plus
        self.minus = ovl.minus
        
        self.i_a = self.a.bisect (l)-1
        self.i_d = self.d.bisect (l)-1
        self.i_plus = self.plus.bisect (l)-1
        self.i_minus = self.minus.bisect (l)-1
        self.is_a = 0
        self.is_d = 0
        self.is_plus = 0
        self.is_minus = 0
        
        self.advance()
        
    cpdef int has_next (self):
        return (self.t0 < self.r)

    cpdef double T0 (self): return self.t0
    cpdef double T1 (self): return self.t1
    
    cpdef int N (self):
        return self.i_a - self.i_d + self.i_plus - self.i_minus

    cpdef advance (self):
        # new t1 == t0
        self.t0 = self.t1
        # advance i
        self.i_a += self.is_a
        self.i_d += self.is_d
        self.i_plus += self.is_plus
        self.i_minus += self.is_minus        
        self.is_plus = self.is_minus = self.is_a = self.is_d = 0
        # new t1 :: either a_next, d_next or r
        cdef double a_next = self.a.item(self.i_a + 1)
        cdef double d_next = self.d.item(self.i_d + 1)
        cdef double plus_next = self.plus.item(self.i_plus + 1)
        cdef double minus_next = self.minus.item(self.i_minus + 1)
        if (self.r < a_next) and (self.r < d_next) and (self.r < plus_next) and (self.r < minus_next):
            self.t1 = self.r
        elif (a_next < d_next) and (a_next < plus_next) and (a_next < minus_next):
            self.is_a = 1
            self.t1 = a_next
        elif (plus_next < minus_next) and (plus_next < d_next):
            self.is_plus = 1
            self.t1 = plus_next
        elif (minus_next < d_next):
            self.is_minus = 1
            self.t1 = minus_next
        else:
            self.is_d = 1
            self.t1 = d_next
        
# utilitites
def cmp_knots (k0, k1):
    return cmp(k0[0], k1[0]) or - cmp(k0[1], k1[1])
