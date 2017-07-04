##  Data structure for managing Arrivals.
##   This needs to be really fast

import warnings
import numpy
import sys
import traceback
import mytime
import pwfun
import hmm
import netutils

from queues cimport Queue
from qnet cimport Qnet

import queues
import cninq
import ninqueue
import ivltreap

cimport cython

# arrival orders 
ARRV_BY_A = 0
ARRV_BY_D = 1

cdef extern from "string.h":
     void *memcpy(void *s1, void *s2, long n)

cdef extern from "Python.h":
    void *PyMem_Malloc(unsigned int n)
    void PyMem_Free(void *p)

cdef class Vars:
    def __init__ (self, map):
        self.map = map.copy()
        
    def duplicate(self):
        return Vars (self.map)
        
    def get (self, var):
        return self.map[var]
        
    def set (self, var, val): 
        self.map[var] = val
        
    def num_vars (self): return len(self.map)

    def allvars (self):
        return self.map.keys()
        
# This is basically a typed tuple for speed
cdef class Event:
    
    def __init__ (self, eid, tid, qid, q, a, s, d, state, obs_a=0, obs_d=0, proc=0, aux=None):
        self.eid = eid
        self.tid = tid
        self.qid = qid
        self.q = q
        self.a = a
        self.s = s
        self.c = max(a,d-s)
        self.d = d
        self.obs_a = obs_a
        self.obs_d = obs_d
        self.state = state
        self.proc = proc if proc else 0
        self.arrv = None
        if aux:
            self.auxiliaries = aux
        else:
            self.auxiliaries = Vars({})
        
    def __repr__ (self):
        return "[EVT %5d: %5.5f .. %5.5f S: %5.5f TID: %5d Q: %s ST: %d PID: %d OBS: %d%d]" % (self.eid, self.a, self.d, self.s, self.tid, self.q.name, self.state, self.proc, self.obs_a, self.obs_d)

    def __dealloc__ (self):
        if self.d_proc_prev != NULL:
            PyMem_Free (self.d_proc_prev)

    def dump (self):
        foo = "[EVT %5d: %23.17g .. %23.17g S: %23.17g TID: %5d Q: %s ST: %d PID: %d  OBS: %d%d " % (self.eid, self.a, self.d, self.s, self.tid, self.q.name, self.state, self.proc, self.obs_a, self.obs_d)
        if self.auxiliaries:
            for k in self.auxiliaries.allvars():
                foo += "%s: %s " % (k, self.auxiliaries.get(k))
        cdef int i
        if self.d_proc_prev != NULL:
            foo += "[TIMES "
            for i from 0 <= i < self.num_proc:
                foo += "%.5f " % self.d_proc_prev[i]
            foo += "]"
        else: foo += " [NO  TIMES]"
        foo += "]"
        return foo
    
    def as_csv (self):
         return "%d %d %d  %.5f %.5f %.5f %.5f %d %d %d %d" % (self.eid, self.tid, self.qid, self.a, self.d, self.s, self.wait(), self.state, self.proc, self.obs_a, self.obs_d)

    def arrival (self): return self.a
    def service (self): return self.s
    def departure (self): return self.d
    def queue (self): return self.q

    # a faster alternative to evt.s = v
    cdef set_service (self, double v):
        self._s = v
        c_old = self.c
        self.reset_commencement()
        if self.arrv: 
            self.arrv.update_commencement_cache (self, c_old, self.c)


    # a faster alternative to evt.s = v
    cdef set_departure (self, double d):
        c_old = self.c
        d_old = self.d
        self.d = d
        self.reset_commencement()
        if self.arrv: 
            self.arrv.update_departure_caches (self, d_old, d)
            self.arrv.update_commencement_cache (self, c_old, self.c)

    def set_state (self, sid): self.state = sid

    cpdef Event duplicate (self):
        return self._dup_inner (1)

    cdef Event _dup_inner (self, int dup_aux):
        cdef Event evt

        aux = self.auxiliaries.duplicate()

        evt = Event(self.eid, self.tid, self.qid, self.q, self.a, self.s, self.d, self.state, obs_a=self.obs_a, obs_d=self.obs_d, proc=self.proc, aux=aux)

        cdef int i
        if self.d_proc_prev != NULL:
            evt.num_proc = self.num_proc
            evt.d_prev = self.d_prev
            evt.init_dtimes()
            memcpy (evt.d_proc_prev, self.d_proc_prev, self.num_proc * sizeof(double))

        return evt

    cdef Event dup_with_structure (self):
        cdef Event evt = self._dup_inner (0)
        evt.prev_byq = self.prev_byq
        evt.next_byq = self.next_byq
        evt.prev_byt = self.prev_byt
        evt.next_byt = self.next_byt
        return evt

    cdef init_dtimes (self):
        self.d_proc_prev = <double*> PyMem_Malloc (self.num_proc * sizeof(double))

    def previous_by_queue (self): return self.prev_byq
    def next_by_queue (self): return self.next_byq
    def previous_by_task (self): return self.prev_byt
    def next_by_task (self): return self.next_byt
    def set_next_by_task (self, Event e): 
        self.next_byt = e
        if e is not None: e.prev_byt = self
    def set_previous_by_task (self, Event e): 
        self.prev_byt = e
        if e is not None: e.next_byt = self
    def set_previous_by_queue (self, Event e): 
        self.prev_byq = e
        if e is not None: e.next_byq = self
    def set_next_by_queue (self, Event e): 
        self.next_byq = e
        if e is not None: e.prev_byq = self
    
    def variable (self, var): return self.auxiliaries.get (var)
    def set_variable (self, var, value): self.auxiliaries.set (var, value)
    def set_all_vars (self, Event e): self.auxiliaries = e.auxiliaries.duplicate()

    def copy_caches (self, Event other):
        self.copy_dtimes (other)

    cdef copy_dtimes (self, Event other):
       cdef int i
       if self.d_proc_prev == NULL: 
           self.num_proc = other.num_proc
           self.init_dtimes()
       for i from 0 <= i < other.num_proc:
           self.d_proc_prev[i] = other.d_proc_prev[i]

    def wait (self): return max(0, self.d - self.a - self.s)

    def commencement (self): return self.c
    def reset_commencement (self):
        # max b/c of numerical issues
        self.c = max(self.a, self.d - self.s)


    cdef int update_from_solution (self, lp_soln) except -1:
        d_new = round (lp_soln['d[%d]' % self.eid], 9)
        a_new = round (lp_soln['a[%d]' % self.eid], 9)
        
        if self.obs_d:
            assert abs(self.d - d_new) < 1e-8, "Error in initializing observed event %s...\n was [%s %s] got from LP [%s %s], diff %s" % (self, self.a, self.d, a_new, d_new, (d_new - self.d))

        if self.obs_a:
            assert abs(self.a - a_new) < 1e-8, "Error in initializing observed event %s...\n was [%s %s] got from LP [%s %s]" % (self, self.a, self.d, a_new, d_new)

        self.d = d_new
        self.a = a_new
        
        return 0
        
    def __cmp__ (self, other):
        c = cmp(self.a, other.a)
        if c == 0: c = cmp(self.d, other.d) 
        if c == 0: c = cmp(self.tid, other.tid)
        return c

    
def sort_by_arrival (evtl):
    evtl.sort(cmp=cmp_by_arrival)
    return evtl

def sort_by_departure (evtl):
    evtl.sort(cmp=cmp_by_departure)
    return evtl
        
def cmp_by_arrival (self, other):
    c = cmp(self.a, other.a)
    return cmp(self.d, other.d) if c == 0 else c
        
def cmp_by_departure (self, other):
    c = cmp(self.d, other.d)
    return cmp(self.a, other.a) if c == 0 else c
 

        
def _get_arrival (evt): return (evt.a, evt.d)

BIG_N = 10000

@cython.final
cdef class Arrivals:
    
    ## TODO: Are the python lists going to be fast enough?
    def __init__ (self, net, cap=BIG_N):
        """nq the number of queues in the net."""
        self.net = net
        nq = len (self.net.queues)
        self.byq = []
        self.final_byq = []
        for i in xrange(nq): 
            self.byq.append(None)
            self.final_byq.append(None)
        self.byt = dict()
        self.events = {}
        self._initial_events = None
        self.queue_orders = [ q.queue_order() for q in self.net.queues ]
        self.queue_orders[0] = ARRV_BY_D # special case for initial queue
        self.ninq = [ cninq.Ninq(cap) for q in range(nq) ]
        self.c2e_cache = [ ninqueue.SortedList() for q in range(nq) ]
        self.a2e_cache = [ ninqueue.SortedList() for q in range(nq) ]
        self.inq = self.setup_inq_cache()
        self.cap = cap

    def __iter__ (self): return self.events.values().__iter__()
    
    def __repr__ (self):
        ret = ""
        for qi in xrange(len(self.byq)):
            q = self.net.queues[qi]
            ret = ret + "QUEUE %s (%d)\n" % (q.name, qi)
            for evt in self.events_of_qid (qi):
                ret = ret + str(evt)
                ret = ret + "\n"
        return ret

    def display (self):
        for qi in xrange(len(self.byq)):
            q = self.net.queues[qi]
            print "QUEUE %s (%d)" % (q.name, qi)
            for evt in self.events_of_qid (qi):
                print evt

    def dump (self):
        ret = ""
        for qi in xrange(len(self.byq)):
            q = self.net.queues[qi]
            ret = ret + "QUEUE %s (%d)\n" % (q.name, qi)
            for evt in self.events_of_qid (qi):
                ret = ret + evt.dump()
                ret = ret + "\n"
        return ret
        
    def as_csv (self):
        ret = "EID TID QID A D S W STATE PROC OBS_A OBS_D CMP\n"
        cdef Event evt
        for evt in self:
            ret = ret + evt.as_csv()
            for k in evt.auxiliaries.allvars():
                ret += " %s" % (evt.auxiliaries.get(k))
            if evt.auxiliaries.num_vars() == 0:
                ret += " NA"
            ret += "\n"

        return ret

    def write_csv (self, f):
        f.write ("EID TID QID A D S W STATE PROC OBS_A OBS_D CMP\n")
        cdef Event evt
        for evt in self:
            f.write (evt.as_csv())
            for k in evt.auxiliaries.allvars():
                f.write (" %s" % (evt.auxiliaries.get(k)))
            if evt.auxiliaries.num_vars() == 0:
                f.write(" NA")
            f.write("\n")

    def all_events (self):
        return self.events.values()
    
    def events_of_qid (self, qid): 
        cdef Event e0 = self.byq[qid]
        result = []
        if not e0: return result

        visited = {}
        while e0.prev_byq: 
            if e0.eid in visited:
                print e0
                print visited.keys()
                raise Exception ("Cycle!!!")
            visited[e0.eid] = True
            e0 = e0.prev_byq
        self.byq[qid] = e0
        
        # add e0 and successors
        result = []
        cdef Event e = e0
        while e:
            result.append (e)
            e = e.next_byq
        
        return result

    def events_of_queue (self, queue):
        qid = self.net.qid_of_queue (queue)
        return self.events_of_qid (qid)
    
    def final_arrival_by_queue (self, queue):
        qid = self.net.qid_of_queue (queue)
        return self.final_arrival_by_qid (qid)

    def final_arrival_by_qid (self, qid):
        cdef Event final
        if self.byq[qid] is None:
            return None
        else:
            final = self.final_byq[qid]
            if final is None:
                final = self.byq[qid]
            while final.next_byq is not None:
                final = final.next_byq
            self.final_byq[qid] = final
            return final

    def first_arrival_by_queue (self, q):
        qid = self.net.qid_of_queue (q)
        return self.first_arrival_by_qid (qid)

    def first_arrival_by_qid (self, qid):
        cdef Event evt = self.byq[qid]
        if not evt: return None
        while evt.prev_byq:
            evt = evt.prev_byq
        return evt
        
    def events_of_task (self, tid):
        return self.byt[tid]
        
    def c2e_range (self, q, l, r):
        qid = self.net.qid_of_queue (q)
        return self.c2e_cache[qid].in_range(l,r)

    def arrivals_in_range (self, q, l, r):
        qid = self.net.qid_of_queue (q)
        return self.a2e_cache[qid].in_range(l,r)

    def evts_in_range (self, q, l, r):
        qid = self.net.qid_of_queue (q)
        nodes = self.inq[qid].all_intersect (l, r)
        return [ n.other for n in nodes ]

    cpdef inline int num_events (self): return len(self.events)
    
    def num_hidden_events (self):
        n = 0
        for e in self:
            if not e.obs_a or not e.obs_d:
                n += 1
        return n
        
    def num_tasks (self): return len(self.byt)
    
    def all_tasks (self): return self.byt.keys()

    def max_eidx (self):
        the_max = -1
        for evt in self: the_max = max (the_max, evt.eid)
        return the_max

    def max_tid (self):
        the_max = -1
        for evt in self: the_max = max (the_max, evt.tid)
        return the_max
        
    def initial_events (self):
        if not self._initial_events: 
            l = []
            for tid in self.byt:
                evts = self.byt[tid]
                if evts:
                    l.append (evts[0])
            self._initial_events = l
        return self._initial_events
        
        
    # Assumes event occurs after all events in q
    def append (self, Event evt):
        if evt.tid not in self.byt:
            self.byt[evt.tid] = []
            
        cdef Event prev_evt
        if self.byt[evt.tid]:
            prev_evt = self.byt[evt.tid][-1]
            prev_evt.next_byt = evt
            evt.prev_byt = prev_evt
            
        if self.byq[evt.qid]:
            prev_evt = self.final_arrival_by_queue (evt.q)
            prev_evt.next_byq = evt
            evt.prev_byq = prev_evt        

        self.byt[evt.tid].append (evt)
        self._update_final_byq (evt)
        self._update_ninq (None, evt)
        self.update_commencement_cache (evt, -1, evt.c)
        self.update_arrival_cache (evt, -1, evt.a)
        self.update_inq_cache (None, evt)

        self.events[evt.eid] = evt
        evt.arrv = self
        
    def insert (self, Event evt):
        if evt.tid not in self.byt:
            self.byt[evt.tid] = []
            
        cdef Event prev_evt
        if self.byt[evt.tid]:
            prev_evt = self.byt[evt.tid][-1]
            prev_evt.next_byt = evt
            evt.prev_byt = prev_evt
            
        cdef Event next_evt
        if self.byq[evt.qid]:
            next_evt = self.find_next_byq (evt.q, evt)
            prev_evt = next_evt.prev_byq if next_evt else self.final_arrival_by_queue(evt.q)

            # evt --> next
            evt.next_byq = next_evt
            if next_evt: next_evt.prev_byq = evt

            # prev --> evt
            evt.prev_byq = prev_evt        
            if prev_evt: prev_evt.next_byq = evt

        self.byt[evt.tid].append (evt)
        self._update_final_byq (evt)
        self._update_ninq (None, evt)
        self.update_commencement_cache (evt, -1, evt.c)
        self.update_arrival_cache (evt, -1, evt.a)
        self.update_inq_cache (None, evt)

        self.events[evt.eid] = evt
        evt.arrv = self

    def resort_byq (self, Event evt):
        cdef Event next_evt, prev_evt, old_prev, old_next

        # first remove evt
        old_prev = evt.prev_byq
        old_next = evt.next_byq
        if old_prev: old_prev.next_byq = old_next
        if old_next: old_next.prev_byq = old_prev
        if self.final_byq[evt.qid] is evt:
            self.final_byq[evt.qid] = old_prev
        if self.byq[evt.qid] is evt:
            if old_prev:
                self.byq[evt.qid] = old_prev
            else:
                self.byq[evt.qid] = old_next

        # can't use insert b/c it screws up the byt caches
        next_evt = self.find_next_byq (evt.q, evt)
        prev_evt = next_evt.prev_byq if next_evt else self.final_arrival_by_queue(evt.q)
        # evt --> next
        evt.next_byq = next_evt
        if next_evt: next_evt.prev_byq = evt
        # prev --> evt
        evt.prev_byq = prev_evt        
        if prev_evt: prev_evt.next_byq = evt
        # you'll have to regenerate teh caches after doing this


    def find_next_byq (self, q, Event evt):
#         cdef Event next_evt = self.byq[evt.qid]
#         order = self.queue_orders[evt.qid]
#         if order == ARRV_BY_D:
#             while next_evt and evt.d >= next_evt.d:
#                 next_evt = next_evt.next_byq
#         elif order == ARRV_BY_A:
# #            print "+++++++++"
#             while next_evt and evt.a >= next_evt.a:
#                 next_evt = next_evt.next_byq
        order = self.queue_orders[evt.qid]
        cdef Event evt0 = self.final_arrival_by_qid(evt.qid)
        cdef Event evt1 = None
        if order == ARRV_BY_D:
            while evt0 and evt0.d >= evt.d:
                evt1 = evt0
                evt0 = evt0.prev_byq
        elif order == ARRV_BY_A:
#            print "+++++++++"
            while evt0 and evt0.a >= evt.a:
                evt1 = evt0
                evt0 = evt0.prev_byq
        else:
            raise Exception ("Bad order %d for queue %s" % (order, q))

#        print evt
#        print next_evt

        return evt1

    def applyDiffList (self, diff_list, lik=None):
        cdef Event e_new, e_old
        cdef Event evt_prev, evt_next
        cdef Event e

        # must be done first
        # the dict() thing sucks
        if lik: 
            done = dict()
            for e in diff_list:
                if e.qid not in done:
                    delta = e.q.likelihoodDelta (self, diff_list)
                    done[e.qid] = 1
                    lik += delta

        for e_new in diff_list:
            e_old = self.event (e_new.eid)
            self._update_ninq (e_old, e_new)  # need to do this BEFORE mucking with e_old
            self.update_arrival_cache (e_old, e_old.a, e_new.a)
            self.update_inq_cache (e_old, e_new)

#            print "DL NEW ", e_new, e_new.c
#            print "DL OLD ", e_old, e_old.c

            e_old.a = e_new.a
            e_old.d = e_new.d
            e_old.s = e_new.s
            e_old.proc = e_new.proc
            e_old.copy_dtimes (e_new)

            # arrival order
            if e_new.prev_byq:
                # add evt_prev --> e_new
                evt_prev =  self.event (e_new.prev_byq.eid)
                e_old.prev_byq = evt_prev
                evt_prev.next_byq = e_old
            else:
                e_old.prev_byq = None

            if e_new.next_byq:
                evt_next = self.event (e_new.next_byq.eid)
                e_old.next_byq = evt_next
                evt_next.prev_byq = e_old
            else:
                e_old.next_byq = None
                
            e_old.q.update_caches_for (e_old)

            c_old = e_old.commencement()
            e_old.reset_commencement()
            self.update_commencement_cache (e_old, c_old, e_old.c)
            assert self.c2e_cache[e_old.qid].contains_kv (e_old.c, e_old)
            assert self.a2e_cache[e_old.qid].contains_kv (e_old.a, e_old)

            # TODO: auxiliaries
            
        # Can't do this sooner, b/c there may be temporary cycles in the byq 
        #   list while the above loop
        for e_new in diff_list:
            e_old = self.event (e_new.eid)
            self._update_final_byq(e_old)

        return lik

    def _update_final_byq (self, Event evt):
        cdef Event final
        if not self.byq[evt.qid]: 
            self.byq[evt.qid] = evt
            self.final_byq[evt.qid] = evt
        else:
            final = self.final_byq[evt.qid]
            while final.next_byq:
                final = final.next_byq
            self.final_byq[evt.qid] = final

    def _update_ninq (self, evt_old, evt):
        if evt_old:
            if evt_old.a != evt.a:
                self.ninq[evt.qid].move_arrival(evt_old.a, evt.a)
            if evt_old.d != evt.d:
                self.ninq[evt.qid].move_departure(evt_old.d, evt.d)
        else:
            self.ninq[evt.qid].add_birth_death(evt.a, evt.d)

    def update_commencement_cache (self, Event evt, double c0, double c1):
        cache = self.c2e_cache[evt.qid]
        if cache.contains_key(c0):
            cache.remove (c0, evt.eid)
        cache.add (c1, evt)

    def update_arrival_cache (self, Event evt, double a0, double a1):
        cache = self.a2e_cache[evt.qid]
        if cache.contains_key(a0):
            cache.remove (a0, evt.eid)
        cache.add (a1, evt)

    def setup_inq_cache (self):
        inq = list()
        for q in self.net.all_queues():
            if isinstance(q, queues.QueuePS):
                treap = ivltreap.IntervalTreap()
            else: treap = None
            inq.append (treap)
        return inq

    def update_inq_cache (self, Event evt_old, Event evt_new):        
        ivltr = self.inq[evt_new.qid]
        if ivltr is not None:  # if it's none, the inq cache is disabled
            if evt_old: 
#            print "REM %.5f %.5f %s" % (evt_old.a, evt_old.d, evt_old) #GG
                ivltr.remove (evt_old.a, evt_old.d, evt_old)
            # N.B. applyDiffList updates in place
            evt = evt_old if evt_old else evt_new
#        print "INQ %.5f %.5f %s" % (evt_new.a, evt_new.d, evt) #GG
            ivltr.insert (evt_new.a, evt_new.d, evt)

    def update_departure_caches (self, Event evt, double d0, double d1):
        self.ninq[evt.qid].move_departure (d0, d1)
        if self.inq[evt.qid] is not None:
            self.inq[evt.qid].remove (evt.a, d0, evt)
            self.inq[evt.qid].insert (evt.a, d1, evt)
#        self.inq[evt.qid].validate()

    def regenerate_all_caches (self):
        nq = self.num_queues()
        self.ninq = [ cninq.Ninq(self.cap) for q in range(nq) ]
        self.c2e_cache = [ ninqueue.SortedList() for q in range(nq) ]
        self.a2e_cache = [ ninqueue.SortedList() for q in range(nq) ]
        self.inq = self.setup_inq_cache()
        for evt in self:
            evt.reset_commencement()
            self._update_ninq (None, evt)
            self.update_commencement_cache (evt, -1, evt.commencement())
            self.update_arrival_cache (evt, -1, evt.a)
            self.update_inq_cache (None, evt)

    def ninq_of_queue (self, queue):
        qid = self.net.qid_of_queue (queue)
        return self.ninq[qid]

    # expensive
    def delete_task (self, tid):
        task_list = self.byt[tid][:]
        cdef Event evt, prev, next
        for evt in task_list:
            self.events.pop (evt.eid) # workaround cython bug??
            
            prev = evt.prev_byq
            next = evt.next_byq
            if prev: prev.next_byq = next
            if next: next.prev_byq = prev 

            exemplar = self.byq[evt.qid]
            if exemplar == evt:
                exemplar = evt.prev_byq
                if not exemplar: exemplar = evt.next_byq
                self.byq[evt.qid] = exemplar
            
        del self.byt[tid]
        
        self.recompute_all_service()
                
    def append_and_check (self, Event evt):
        """Adds the given event to this Arrivals object.
            Event must arrive after all other arrivals to this object.
            This method checks that adding evt would not break FIFO constraints.
            If it would, a warning is printed and the event is not added."""
        clear_time = evt.q.when_clear (self)
        if clear_time > evt.d:
            warnings.warn ("Arrivals: Disregarding event %s (breaks FIFO assumptions)" % evt)
        else:
            self.append (evt)

    def qnet (self): return self.net

    def num_queues (self): return len(self.byq)
    
    cpdef inline Event event (self, int eid): 
        return self.events[eid]

    def is_initial (self, evt):
        return self._is_initial (evt)
        
    cdef _is_initial (self, Event evt):
        return not evt.prev_byt
        
    def is_final (self, evt):
        return self._is_final (evt)
        
    cdef _is_final (self, Event evt):
        return not evt.next_byt
    
    def is_last_in_queue (self, evt):
        cdef Event e = evt
        return not e.next_byq
        
    def recompute_all_service (self):
        cdef Event evt
        for evt in self:
            evt.q.recompute_service (evt)
            evt.reset_commencement()
        # regenerate the commencement caches
        self.c2e_cache = [ ninqueue.SortedList() for q in range(self.num_queues()) ]  #
        for evt in self:
            self.update_commencement_cache (evt, -1, evt.commencement())

    def subset (self, fn, adapt_fn=None):
        cdef Arrivals a_new
        cdef Event evt, evt_new
        
        a_new = Arrivals (self.net, cap=self.cap)
        for tid in self.byt.keys():
            for evt in self.events_of_task (tid):
                evt_new = evt.duplicate()
                if not fn(self,evt):
                    evt_new.obs_a = evt_new.obs_d = 0
                a_new.insert (evt_new)

        return a_new
        
    def subset_by_task (self, pct, adapt_fn=None):
        """Return a new copy of this Arrivals object in which about PCT %% of the tasks are observed."""
        to_inc = dict()
        for tid in self.byt:
            if numpy.random.rand() > pct:
                to_inc[tid] = 1.0
        return self.subset (netutils.not_in_set2 (to_inc), adapt_fn=adapt_fn)

    # if function(arrv, tid) is true, put tid in return
    def subset_by_task_fn (self, fn, adapt_fn=None):
        cdef Arrivals a_new
        cdef Event evt, evt_new
        
        a_new = Arrivals (self.net, cap=self.cap)
        for tid in self.byt:
            if fn(self, tid):
                for e in self.events_of_task (tid):
                    a_new.insert (e.duplicate())

        return a_new
    
    # Updates all of the by_task pointer.  Use this if you've rearranged the task orders.
    def recompute_caches (self):
        cdef Event evt
        dct_new = dict()
        for k,v in self.byt.iteritems():
            evt = v[0]
            assert evt.tid == k
            while evt.prev_byt: evt = evt.prev_byt

            lst = [evt]
            while evt.next_byt:
                evt = evt.next_byt
                lst.append (evt)

            dct_new[k] = lst
        self.byt = dct_new

    def duplicate (self): return self.subset (netutils.all_true)

    def validate (self):
        cdef Event evt0, evt1
        try:
            evts_byq = dict()
            evts_byt = dict()
            for qi in range(len(self.byq)):
               qevts = self.events_of_qid (qi)
               for evt in qevts: evts_byq[evt.eid] = 1
               for evt0, evt1 in zip(qevts, qevts[1:]):
                    evts_byq[evt0.eid] = 1
                    evts_byq[evt1.eid] = 1
                    assert evt0.s >= 0.0, "Invalid (negative) service time %s for:\n  %s\n  %s" % (evt0.s, format_evt (evt0.prev_byq), format_evt(evt0))
                    assert evt1.s >= 0.0, "Invalid (negative) service time %s for:\n  %s\n  %s" % (evt1.s, evt1.prev_byq.dump(), evt1.dump())
                    assert evt0.next_byq == evt1, "Next_byq pointers don't match\n  %s\n ==> next_byq\n%s\nbut\n  %s\n  ==>prev_byq\n  %s" % (evt0, evt0.next_byq, evt1, evt1.prev_byq)
                    assert evt1.prev_byq == evt0, "Prev_byq pointers don't match:\n %s\n ==> prev_byq\n%s\nbut\n  %s\n  ==>next_byq\n  %s" % (evt1, evt1.prev_byq, evt0, evt0.next_byq)
                    assert evt0.eid in self.events, "Event not in event list:\n  %s" % evt0
                    assert self.events[evt0.eid] == evt0, "Event list mismatch\n  Expected: %s\n  Was: %s" % (evt0, self.events[evt0.eid])
            for tid, task_evts in self.byt.iteritems():
                for evt in task_evts: evts_byt[evt.eid] = 1
                for evt0, evt1 in zip(task_evts, task_evts[1:]):
                    assert evt0.next_byt == evt1, "Next_byt pointers don't match\nFROM   %s\nPOINTS  %s\nARRV   %s" % (evt0, evt0.next_byt, evt1)
                    assert evt1.prev_byt == evt0, "Prev_byt pointers don't match\n  %s\n  %s" % (evt0, evt1)
                    assert tid == evt0.tid, "TID pointers don't match (qnet has %s)\n  %s" % (tid, evt0)
                    assert tid == evt1.tid, "TID pointers don't match (qnet has %s)\n  %s" % (tid, evt1)
                    assert abs(evt0.d - evt1.a) < 1e-15 , "Mismatched arrivals (D_prev: %s  A: %s)\n  %s\n  %s" % (evt0.d, evt1.a, evt0, evt1)
            for evt0 in self:
                evt0.q.validate_service (evt0)
                if evt0.next_byt:
                    assert evt0.next_byt.tid == evt0.tid, "TID pointers don't match, but are linked\n  %s\n  %s" % (evt0, evt0.next_byt)
            for q in self.net.queues:
                q.validate_caches (self)
            for eid in self.events.iterkeys():
                assert eid in evts_byq, "Event in event list but not in by_q cache:\n  %s" % self.events[eid]
                assert eid in evts_byt, "Event in event list but not in by_t cache:\n  %s" % self.events[eid]
            for eid in evts_byq.iterkeys():
                assert eid in self.events, "Event in by_q cache but not events %s" % eid
            for eid in evts_byt.iterkeys():
                assert eid in self.events, "Event in by_t cache but not events %s" % eid
            for evt in self:
                ninq = self.ninq [evt.qid]
                assert ninq.contains_arrival (evt.a), \
                     "NINQ wrong, does not contain arrival\n%s\n%.15f\n%s" % (evt, evt.a, ninq)
                assert ninq.contains_departure (evt.d), \
                     "NINQ wrong, does not contain departure\n%s\n%.15f\n%s" % (evt, evt.d, ninq)
                if evt.s != 0:
                    assert ninq.N (evt.commencement()) > 0, "Can't find event in ninq!\nEVT %s\nE.C %.15f\nNINQ %s" % (evt, evt.commencement(), ninq)            
                if evt.state == self.net.fsm.initial_state():
                    assert evt.a == 0, "Initial event with bad arrival\n%s" % evt

            for inq in self.inq:
                if inq: inq.validate()

            check_commencement_cache (self)
            
        except AssertionError, e:
            self.write_csv (sys.stdout)
            print e
            traceback.print_exc()
            raise e

    def report (self):
        ne = self.num_events()
        nt = self.num_tasks()
        print "ARRIVALS NE %d  NT %d" % (ne, nt)
        for treap in self.inq:
            if treap: treap.report()

def format_evt (evt):
    if evt is None:
        return "None"
    else:
        return evt.dump()

EPS = 1e-5

def check_commencement_cache (Arrivals arrv):
    for qi in range(len(arrv.byq)):
        c2e = arrv.c2e_cache[qi]
        evtl = arrv.events_of_qid (qi)
        assert len(c2e) == len(evtl), "Huh? EVTL %s  C2E %s" % ("\n".join(map(str,evtl)), c2e)
        for e in evtl:
            c = e.commencement()
            assert c2e.contains_key (c)
            assert c2e.in_range (c-EPS, c+EPS).index(e) >= 0
