DEBUG = 0


cdef extern from "Python.h":
    void *PyMem_Malloc(unsigned int n)
    void PyMem_Free(void *p)

cdef extern from "math.h":
    double exp(double x)
    double log(double x)
    enum: INFINITY

EPS = 1e-14

import numpy
from numpy import random

import qnet # for has_tau
import arrivals
import sampling
from sampling cimport _roll_die_unnorm
from arrivals cimport Event
import distributions # for stub
import pwfun
import netutils
import ninqueue
from cninq cimport Ninq, Overlay, OverlayIterator

import traceback
import math

# control over proposal functions

PRP_NAMES = [ "twos", "uniform", "bar", "zoid" ]
TWOS_PROPOSAL = 0
UNIFORM_PROPOSAL = 1
BAR_PROPOSAL = 2
ZOID_PROPOSAL = 3
proposal = TWOS_PROPOSAL

def set_proposal(prp):
    global proposal
    proposal = prp
    print "MH proposal = ", PRP_NAMES[proposal]

def proposal_names(): return PRP_NAMES[:]

# end parameters

cdef class FactorTemplate
cdef class DistributionTemplate

cdef class Queue: 

    def __init__ (self, name, service):
        self.name = name
        self.service = service

    def __repr__ (self):
        return "[QUEUE: %s S: %s ]" % (self.name, self.service)

    def sample_service (self): return self.service.sample_service ()
    
    cdef Pwfun arrivalLik (self, Arrivals arrv, Event evt): raise NotImplementedError()    
    
    cdef Pwfun departureLik (self, Arrivals arrv, Event evt): raise NotImplementedError()

    def pyArrivalLik (self, arrv, evt):
        return self.arrivalLik (arrv, evt)

    def pyDepartureLik (self, arrv, evt):
        return self.departureLik (arrv, evt)

    
    cdef diffListForArrival (self, Event evt, double d): raise NotImplementedError()
    cdef diffListForDeparture (self, Event evt, double d): raise NotImplementedError()

    def pyDiffListForArrival (self, evt, d):
        return self.diffListForArrival (evt, d)
    def pyDiffListForDeparture (self, evt, d):
        return self.diffListForDeparture (evt, d)
    
    cpdef double likelihoodDelta (self, Arrivals arrv, diffList) except *:
        cdef Event e_new, e_old
        cdef double lik = 0, v_old, v_new
        for e_new in diffList:
            if e_new.q == self:
                if e_new.s < 0: 
                    return -numpy.inf
                e_old = arrv.event (e_new.eid)
                v_old = self.service_lpdf (e_old)
                v_new = self.service_lpdf (e_new)
                lik = lik + v_new - v_old
#                print "%d %.7f %.7f %.7f" % (e_new.eid, v_new, v_old, v_new - v_old) 
        return lik

    cdef double allEventLik (self, Arrivals arrv) except *: 
        cdef double result = 0, current
        cdef Event evt
#        print "////// aEL" 
        for evt in arrv.events_of_queue (self):
            if evt.s < 0:
                print "Arrv invalid: ", evt
                return -INFINITY
            current = self.service_lpdf(evt)  
#            print "%d %.7f %.7f" % (evt.eid, evt.s, current) 
            result += current
#        print "\\\\\\\\\\\\ aEL" 
        return result


    cpdef recompute_service (self, Event e): raise NotImplementedError()
    cpdef recompute_departure (self, Arrivals arrv, Event e): raise NotImplementedError()
    cdef double departure_of_service  (self, Event e, double s): raise NotImplementedError()
    
    cpdef double service_lpdf (self, Event e):
        return self.service.serviceLikelihood(e)(e.s)
        
    def validate_service (self, e): raise NotImplementedError()

    def previous_departing (self, e): raise NotImplementedError()
    def previous_arriving (self, Event e): return e.prev_byq

    def initialize_auxiliaries (self, e):
        self.service.initialize_auxiliaries (e)
        
    def resample_auxiliaries (self, Event e):
        self.service.resample_auxiliaries (e)

    def update_caches_for (self, Event e):
        pass

    def validate_caches (self, arrv): pass

    def queue_order (self):
        return arrivals.ARRV_BY_A

    # Assumes EVT_LIST is a list of Events that are sorted by arrival time for this queue,
    #  but have not been served yet.
    def select_job_for_service (self, evt_list): raise NotImplementedError

    def is_markov (self):
        f = (<DistributionTemplate>self.service).f
        return isinstance(f, distributions.Exponential)

# forward refs for queueGG1
cdef class ArrivalFn
cdef class DepartureFn
cdef class FinalDepartureFn



cdef class QueueGG1 (Queue):

    def __repr__ (self):
        return "<QUEUE TYPE='GG1' NAME='%s'>\n<SERVICE>%s</SERVICE>\n</QUEUE>" % (self.name, self.service)

    def as_yaml (self):
        return "  - { name: %s, type: GG1, service: %s }" % (self.name, self.service.as_yaml())

    def when_clear (self, Arrivals arrv):
        cdef Event last_evt = arrv.final_arrival_by_queue (self)
        return last_evt.d if last_evt else 0.0
        
    def previous_departing (self, Event e):
        return e.prev_byq # previous in arrival order
        
    cpdef recompute_service (self, Event e1):
        cdef Event e0, e2
        e0 = e1.prev_byq
        e2 = e1.next_byq
        
        if e0:
            e1.s = e1.d - max(e1.a, e0.d)
        else:
            e1.s = e1.d - e1.a
            
        if e2:
            e2.s = e2.d - max(e2.a, e1.d)
    
    cdef double departure_of_service  (self, Event e, double s):
        cdef Event byq = e.prev_byq
        cdef double d_prev = byq.d if byq else 0.0
        return s + max (e.a, d_prev)
    
    cpdef recompute_departure (self, Arrivals arrv, Event e1):
        cdef Event e_prev = e1.prev_byq
        cdef double d_prev = e_prev.d if e_prev else 0
        e1.set_departure (e1.s + max(e1.a, d_prev))
    
    def validate_service (self, evt):
        cdef Event e = evt
        d_prev = (e.prev_byq.d if e.prev_byq else 0)
        expected = e.d - max(e.a, d_prev)
        assert abs(e.s - expected) < 1e-10, "Invalid service time (expected %s was %s)\n  %s\n  %s" % (expected, e.s, e.prev_byq, e)
        if e.prev_byq:
            assert e.a >= e.prev_byq.a, "Invalid arrival times for events PREV: %s CURR: %s\n  %s\n  %s" % (e.prev_byq.a, e.a, e.prev_byq, e)
            
    cdef Pwfun arrivalLik (self, Arrivals arrv, Event evt): 
        cdef Event e_prev, e_next
        e_prev = evt.prev_byq
        e_next = evt.next_byq
        
        if e_prev:
            d_prev = e_prev.d
            L = e_prev.a
        else:
            L = 0
            d_prev = 0
            
        if e_next:
            U = min(e_next.a, evt.d)
        else:
            U = evt.d
            
        cdef ArrivalFn afn
        afn = ArrivalFn (self.service, evt, evt.d, d_prev) 
        if DEBUG: 
            print afn
            print "@ L(=%.17g)  (diff=%.17g) is %s or %s" % (L, L-afn.L(), afn.f1(L), afn.f2(L))
        
        if d_prev < L:
            return Pwfun ([L, U], [ afn.f2 ], [ afn.df2 ])
        elif U < d_prev:
            return Pwfun ([L, U], [ afn.f1 ], [ afn.df1 ])
        else:
            return Pwfun ([L, d_prev, U], [ afn.f1, afn.f2 ], [ afn.df1, afn.df2 ])
        

    cdef Pwfun departureLik (self, Arrivals arrv, Event evt):
        cdef Event e_prev, e_next
        cdef double L, A, U
        
        e_prev = evt.prev_byq
        e_next = evt.next_byq
        
        if e_prev: 
            L = max(evt.a, e_prev.d) 
        else: 
            L = evt.a            

        cdef DepartureFn dfun
        cdef FinalDepartureFn fdfun
        
        if e_next:             
            A = e_next.a
            U = e_next.d 
            dfun = DepartureFn (self.service, evt, e_next, L, A, U)
            if DEBUG: 
                print dfun
                print "@ L(=%.17g)  (diff=%.17g) is %s or %s" % (L, L-dfun.L(), dfun.f1(L), dfun.f2(L))
            
#            print "Creating departure fn:\n  %s\nL: %8.4f A: %8.4f U: %8.4f\ndfun: %s\n" % (format_six(evt), L, A, U, dfun)
            if A < L:
                return Pwfun([L, U], [dfun.f2], [dfun.df2])
            elif U < A:                
                return Pwfun([L, U], [dfun.f1], [dfun.df1])
            else:
                return Pwfun([L, A, U], [ dfun.f1, dfun.f2 ], [ dfun.df1, dfun.df2 ])
                
        else:
            fdfun = FinalDepartureFn (self.service, evt, L)
            max_x = numpy.inf
            return Pwfun ([L, max_x], [ fdfun.f ], [ fdfun.df ])
           
    cdef diffListForArrival (self, Event evt, double a):
        cdef Event e_new = evt.dup_with_structure()
        e_new.a = a
        self.recompute_service (e_new)
        return [e_new]
        
    cdef diffListForDeparture (self, Event evt, double d):
        cdef double d_prev = evt.prev_byq.d if evt.prev_byq else 0.0
        cdef Event e_new = evt.dup_with_structure()
        e_new.d = d
        e_new.s = d - max(e_new.a, d_prev)
        lst = [e_new]
        
        cdef Event byq = e_new.next_byq
        if byq:
            byq = byq.dup_with_structure()
            byq.prev_byq = e_new
            e_new.next_byq = byq
            byq.s = byq.d - max(byq.a, d)
            lst.append (byq)
            
        return lst      

    def select_job_for_service (self, evt_list):
        return evt_list[0]
          

# used by GG1
cdef class DepartureFn:
    cdef object sfn1
    cdef object sfn2
    cdef object dsfn1_dx
    cdef object dsfn2_dx
    cdef double b_e
    cdef double a_next
    cdef double d_next

    def __init__ (self, FactorTemplate tmpl, Event evt1, Event evt2, b_e, a_next, d_next):
        self.sfn1 = tmpl.serviceLikelihood (evt1)
        self.sfn2 = tmpl.serviceLikelihood (evt2)
        self.dsfn1_dx = tmpl.serviceDerivative (evt1)
        self.dsfn2_dx = tmpl.serviceDerivative (evt2)
        self.b_e = b_e
        self.d_next = d_next
        self.a_next = a_next

    # debugging
    cdef double L(self): return self.b_e

    def f1 (self, d):
        ret = self.sfn1 (d - self.b_e) + self.sfn2 (self.d_next - self.a_next)
#        print "depart_f1 %20.17f %20.17f  (%20.17g %20.17g %20.17g)" % (d, ret, self.b_e, (self.d_next - self.a_next), (d - self.b_e))      
        return ret

    def f2 (self, d):
        ret = self.sfn1 (d - self.b_e) + self.sfn2 (self.d_next - d)
#        print "depart_f2 %8.4f %8.4f (%8.4f %8.4f)" % (d, ret, self.b_e, self.d_next)      
        return ret

    def df1 (self, d):        
        return self.dsfn1_dx (d - self.b_e)

    def df2 (self, d):
        return self.dsfn1_dx (d - self.b_e) - self.dsfn2_dx (self.d_next - d)

    def __repr__ (self):
        return "<D_FN(x): pdf(x - %.17g) + pdf(%.17g - max(x, %.17g))>" % (self.b_e, self.d_next, self.a_next)



cdef class FinalDepartureFn:
    cdef object sfn
    cdef object dsfn_dx
    cdef double b_e

    def __init__ (self, FactorTemplate tmpl, Event evt, b_e):
        self.sfn = tmpl.serviceLikelihood (evt)
        self.dsfn_dx = tmpl.serviceDerivative (evt)
        self.b_e = b_e

    def f (self, d):        
        return self.sfn (d - self.b_e)

    def df (self, d):       
        return self.dsfn_dx (d - self.b_e)


cdef class ArrivalFn:
    cdef object sfn
    cdef object dsfn_dx
    cdef double d
    cdef double d_prev

    def __init__ (self, FactorTemplate tmpl, Event evt, d, d_prev):
        self.sfn = tmpl.serviceLikelihood (evt)
        self.dsfn_dx = tmpl.serviceDerivative (evt)
        self.d = d
        self.d_prev = d_prev

    cdef double L(self): return self.d

    # full function: afn(a) = f_s (d_e - max (a, d_prev))

    # case 1: a < d_prev 
    def f1 (self, a):
        ret = self.sfn (self.d - self.d_prev)
#        print "arrival_f1 %8.4f %s %8.4f" % (a, (self.d-self.d_prev), ret)      
        return ret

    def df1 (self, a): return 0.0

    # case 2: d_prev <= a
    def f2 (self, a):
        ret = self.sfn (self.d - a)
#        print "arrvl_f2 %8.4f %s %8.4f" % (a, self.d - a, ret)
#        print self.s
        return ret

    def df2 (self, a):
        return 0 # TODO
#        return -self.sfn.dx_lpdf (self.d - a)

    def __repr__ (self):
        return "<ARRIVAL_FN(a): pdf(%s - max(%s, a))  />" % (self.d, self.d_prev)


## Used by GGk.  Should probably backport these to GG1

cdef class ArrivalProposal:
    ## PWFUN for  g(a) = f_service (d_e - max(a, d_prev))
    ## Currently the discontinuity is not explicitly represented.
    
    cdef object s_pdf
    cdef object dpdf_dx
    cdef double d
    cdef double d_prev
    
    def __init__ (self, DistributionTemplate tmpl, Event evt, d_prev):
        self.d = evt.d
        self.d_prev = d_prev
        self.s_pdf = tmpl.serviceLikelihood (evt)    
        self.dpdf_dx = tmpl.serviceDerivative (evt)
        
    def f (self, a): 
#        print "pdf(%f - max(x:%f, %f)) = %f" % (self.d, a, self.d_prev, self.s_pdf (self.d - max(a, self.d_prev)))
        return self.s_pdf (self.d - max(a, self.d_prev))
        
    def df (self, a):
        if a < self.d_prev:
            return 0.0
        else:
            return -self.dpdf_dx(self.d - a)
    
    def proposal (self):
        return Pwfun ([0, self.d_prev, self.d], [self.f] * 2, [self.df] * 2)

cdef class DepartureProposal:
    ## PWFUN for g(d) = f_service (d - max(a_e, d_prev))
    
    cdef object s_pdf
    cdef object dpdf_dx
    cdef double b_e

    def __init__ (self, DistributionTemplate tmpl, Event evt, double b_e):
        self.b_e = b_e
        self.s_pdf = tmpl.serviceLikelihood (evt)
        self.dpdf_dx = tmpl.serviceDerivative (evt)
        
    def f(self, d):
#        print "pdf(x:%f - %f) = %f" % (d, self.b_e, self.s_pdf (d - self.b_e))
        return self.s_pdf (d - self.b_e)
        
    def df(self, d):
        return self.dpdf_dx (d - self.b_e)
    
    def proposal(self):
        return Pwfun ([self.b_e, numpy.inf], [self.f], [self.df])

cdef class QueueGGk (Queue):

    cdef readonly int k

    def __init__ (self, name, service, k):
        Queue.__init__ (self, name, service)
        self.k = k
        
    def __repr__ (self):
        return "<QUEUE TYPE='GGk' NAME='%s' K='%d'><SERVICE>%s</SERVICE>\n</QUEUE>" % (self.name, self.k, self.service)

    def as_yaml (self):
        return "  - { name: %s, type: GGk, processors: %d, service: %s }" % (self.name, self.k, self.service.as_yaml())

    cdef Pwfun arrivalLik (self, Arrivals arrv, Event evt):
        # IDEA: Return a Pwfun f_{service} (d_e - max{a_e, c_e})
        cdef Event pred = self.prev_by_proc (evt)
        cdef double d_prev = pred.d if pred else 0.0
        return ArrivalProposal (self.service, evt, d_prev).proposal()

    cdef Pwfun departureLik (self, Arrivals arrv, Event evt):
        # Return a Pwfun f_{service} (d_e - max{a_e, c_e})
        cdef Event pred = self.prev_by_proc (evt)
        cdef double d_prev = pred.d if pred else 0.0
        return DepartureProposal (self.service, evt, max(evt.a, d_prev)).proposal()

    cdef diffListForDeparture (self, Event evt, double d):
        lst = []

        cdef Event evt_new = evt.dup_with_structure()
        evt_new.d = d
        self.recompute_service (evt_new)
        lst.append (evt_new)

        cdef int has_changed = 1
        cdef int proc
        cdef int proc_min
        cdef double entry_time, d_prev

        cdef Event evt_prev, evt_cur = evt.next_byq

        while evt_cur != None and has_changed:
            evt_new = evt_cur.dup_with_structure()
            evt_prev = <Event>lst[-1]
            evt_new.prev_byq = evt_prev
            evt_prev.next_byq = evt_new

            # super-fast (hopefully) version of recompute service
            has_changed = 0
            entry_time = INFINITY
            for proc in range(self.k):
                d_prev = evt_prev.d_proc_prev[proc]
                if proc == evt_prev.proc: 
                    d_prev = evt_prev.d
                if d_prev != evt_new.d_proc_prev[proc]:
                    has_changed = 1
                    evt_new.d_proc_prev[proc] = d_prev
                if evt_new.d_proc_prev[proc] < entry_time:
                    entry_time = evt_new.d_proc_prev[proc]
                    proc_min = proc

            evt_new.s = evt_new.d - max(evt_new.a, entry_time)
            evt_new.proc = proc_min

            lst.append ( evt_new )
            evt_cur = evt_cur.next_byq
            
        return lst

        
    cdef diffListForArrival (self, Event evt, double a):
        lst = []
        cdef int procs_stabilized = 0
                
        cdef Event this_evt, first_evt, last_evt, evt_new

        # set first_evt == one before first event whose proc might be affected, last_evt to one after the last such evt
        cdef double a_min = min (a, evt.a)
        cdef double a_max = max (a, evt.a)
        
        first_evt = evt
        while first_evt.prev_byq and a_min < first_evt.a:
            first_evt = first_evt.prev_byq
            
        # now set last_evt
        last_evt = first_evt if first_evt else evt
        while last_evt and last_evt.a <= a_max:
            last_evt = last_evt.next_byq

        # N.B. if arrival order unchanged, first_evt = evt.prev_byq and last_evt = evt.next_byq
        
        # now duplicate all the set events
        this_evt = first_evt
        while this_evt and this_evt != last_evt:            
            evt_new = this_evt.dup_with_structure()
            if this_evt == evt: evt_new.a = a
            lst.append (evt_new)
            this_evt = this_evt.next_byq            
        lst.sort (key=get_arrival)
        
        # assertion
        if len(lst) == 0:
            print evt
            print "A_NEW: ", a
            raise Exception ("List is null")

        # Now that it's sorted, reset the next links, etc.
        cdef Event e0, e1
        for i from 0 <= i < len(lst)-1:
            e0 = lst[i]
            e1 = lst[i+1]
            e1.prev_byq = e0
            e0.next_byq = e1
        (<Event>lst[0]).prev_byq = first_evt.prev_byq  # special case, b/c don't want to change next_byq of event in arrivals
        (<Event>lst[-1]).next_byq = last_evt
        
        # recompute the procs
        procs_stable = [1] * self.k
        evt = lst[0]
        cdef int p_old
        for evt in lst:
            p_old = evt.proc
            self.recompute_service (evt)
            if p_old != evt.proc:
                procs_stable[p_old] = 0            
                procs_stable[evt.proc] = 0
            else:
                procs_stable[p_old] = 1

        # extend list if necessary
        cdef Event prev_new = (<Event>lst[-1])
        cdef Event evt_old = prev_new.next_byq
        while evt_old and not all(procs_stable):
            evt_new = evt_old.dup_with_structure()
            evt_new.prev_byq = prev_new
            prev_new.next_byq = evt_new
            p_old = evt_old.proc
            self.recompute_service (evt_new)
            if p_old != evt_new.proc:
                procs_stable[p_old] = 0            
                procs_stable[evt_new.proc] = 0
            else:
                procs_stable[p_old] = 1
            lst.append (evt_new)
            prev_new = evt_new
            evt_old = evt_new.next_byq

        return lst    

    def when_clear (self, Arrivals arrv):
        cdef Event last_evt = arrv.final_arrival_by_queue (self)
        if not last_evt: return 0.0
        self._compute_predecessor (last_evt)
        times = dbl2array (last_evt.d_proc_prev, self.k)
        times[last_evt.proc] = last_evt.d
        return min(times)
        
    def previous_departing (self, e):
        return self.prev_by_proc (e)
        
    def select_job_for_service (self, evt_list):
        return evt_list[0]

    # Computes the departure and processor of Event e,
    #  assuming that the following is already correct: 
    #  (a) arrival and departure for this event, and (b) all info for other events 
    cpdef recompute_service (self, Event e):
        self._compute_predecessor (e)
        e.s = e.d - max(e.a, e.d_prev)

    cdef double departure_of_service  (self, Event e, double s):
        cdef Event pred = self.prev_by_proc (e)
        d_prev = pred.d if pred else 0.0
        return s + max(e.a, d_prev)
        
    cdef Event prev_by_proc (self, Event e):
        cdef Event pred = e.prev_byq
        while pred:
            if pred.proc == e.proc: return pred
            pred = pred.prev_byq
        return None
        
    # Computes the departure and processor of Event e,
    #  assuming that the following is already correct: 
    #  (a) arrival and service for this event, and (b) all info for other events 
    cpdef recompute_departure (self, Arrivals arrv, Event e): 
        cdef Event pred
        self._compute_predecessor (e)
        e.set_departure (e.s + max(e.a, e.d_prev))


    ## utilities
        
    # recomputes the proc-clear-times for each queue, assuming that everything is
    #  correct for the last event
    cdef _compute_predecessor (self, Event e):
        self._recompute_cache (e)
        cdef Event prev = <Event>e.prev_byq        

        if prev is None:
            e.d_prev = 0
            e.proc = 0
        else:
            e.d_prev = famin(e.d_proc_prev, self.k, &e.proc)

    cdef setup_d_proc_prev (self, Event e):
        if e.d_proc_prev == NULL:
            e.d_proc_prev = <double*> PyMem_Malloc (self.k * sizeof(double))
            e.num_proc = self.k

    cdef int _recompute_cache (self, Event e):
        self.setup_d_proc_prev (e)
        cdef double* times = e.d_proc_prev
        cdef double* times_prev

        cdef Event prev = <Event>e.prev_byq        
        cdef int i
        cdef int allEq = 1

        if prev is None:
            for i from 0 <= i < self.k: times[i] = 0
            e.d_prev = 0
            e.proc = 0
            allEq = 0
        else:
            if prev.d_proc_prev == NULL: 
                self._recompute_cache (prev)
            for i from 0 <= i < self.k:
                if i != prev.proc:
                    if times[i] != prev.d_proc_prev[i]: 
                        allEq = 0
                    times[i] = prev.d_proc_prev[i]
            if prev.d != times[prev.proc]: 
                allEq = 0
            times[prev.proc] = prev.d

        return allEq

            
    def update_caches_for (self, Event e):
        self._recompute_cache (e)

        cdef double* times
        cdef int i 
        cdef int alleq = 0

        cdef Event next = e.next_byq
        while (next is not None) and (alleq == 0):
            alleq = self._recompute_cache (next)
            next = next.next_byq

    def validate_service (self, Event e): 
        cdef Event pred = self.prev_by_proc (e)
        if pred:
            assert abs(e.s - (e.d - max(e.a, pred.d))) < 1e-10, "Departure mismatch (was %.15f expected %.15f):\n  %s\n  %s\n" % (e.s, e.d - max(e.a, pred.d), pred.dump(), e.dump())
        else:
            assert abs(e.s - (e.d - e.a)) < 1e-10, "Departure mismatch (no predecssor): %s" % e.dump()
        if self.k == 1:
            prev_byq = e.prev_byq
            if prev_byq:
                assert prev_byq.a <= e.a, "Prev_byq arrival mismatch\n  %s\n  %s\n" % (prev_byq.dump(), e.dump())
                assert prev_byq.d <= e.d, "Prev_byq departure mismatch\n  %s\n  %s\n" % (prev_byq.dump(), e.dump())
            next_byq = e.next_byq
            if next_byq:
                assert e.a <= next_byq.a, "GG1q arrival mismatch\n  %s\n  %s" % (e, next_byq)
                assert e.d <= next_byq.d, "GG1q departure mismatch\n  %s\n  %s" % (e, next_byq)

    def validate_caches (self, Arrivals arrv):
        cdef Event evt = arrv.first_arrival_by_queue (self)
        cdef int i
        times = [ 0.0 ] * self.k

        while evt:
            if evt.d_proc_prev:
                my_times = dbl2array (evt.d_proc_prev, self.k)
                for i from 0 <= i < self.k:
                    assert times[i] == my_times[i], "Bad times for event. Expected %s, was%s\n %s" % (times[i], my_times[i], evt.dump())
            times[evt.proc] = evt.d
            evt = evt.next_byq


# Delay station

cdef class ServiceClosure:
  cdef double t0
  cdef object sfn
  cdef object dsfn
  def __init__ (self, FactorTemplate tmpl, Event evt, double t0):
      self.sfn = tmpl.serviceLikelihood (evt)
      self.dsfn = tmpl.serviceDerivative (evt)
      self.t0 = t0

  def p (self, double t1):
      return self.sfn (abs(t1-self.t0))

  def dp_dr (self, double t1):
      # if this is a function of the departure time, then t1 >= self.t0
      cdef int sgn = 1 if t1 >= self.t0 else -1
      return sgn * self.dsfn (abs(t1-self.t0))


cdef class DelayStation (Queue):

    def __repr__ (self):
        return "<QUEUE TYPE='DELAY' NAME='%s'>\n<SERVICE>%s</SERVICE>\n</QUEUE>" % (self.name, self.service)

    cdef Pwfun arrivalLik (self, Arrivals arrv, Event evt): 
       cdef ServiceClosure sc = ServiceClosure (self.service, evt, evt.d)
       return Pwfun ([0.0, evt.d], [ sc.p ], [ sc.dp_dr ])

    cdef Pwfun departureLik (self, Arrivals arrv, Event evt): 
       cdef ServiceClosure sc = ServiceClosure (self.service, evt, evt.a)
       return Pwfun ([evt.a, numpy.inf], [ sc.p ], [ sc.dp_dr ])

    cdef diffListForArrival (self, Event evt, double a):
       e1 = evt.dup_with_structure ()
       e1.a = a
       e1.s = e1.d - a
       return [e1]

    cdef diffListForDeparture (self, Event evt, double d):
       e1 = evt.dup_with_structure ()
       e1.d = d
       e1.s = d - e1.a
       return [e1]

    cpdef recompute_service (self, Event e):
       e.s = e.d - e.a

    cpdef recompute_departure (self, Arrivals arrv, Event e):
       e.set_departure (e.a + e.s)

    def validate_service (self, e):
        assert e.d == e.a + e.s, "Departure/service mismatch in delay station\n  %s" % e
           
    def select_job_for_service (self, evt_list):
        return evt_list[0]

    def when_clear (self, Arrivals arrv):
        # always room for one more!
        return 0

# GG1 queue that selects which event to enter service randomly

cdef class QueueGG1R (Queue):

    def __repr__ (self):
        return "<QUEUE TYPE='GG1R' NAME='%s'>\n<SERVICE>%s</SERVICE>\n</QUEUE>" % (self.name, self.service)

    def as_yaml (self):
        return "  - { name: %s, type: GG1R, service: %s }" % (self.name, self.service.as_yaml())

    def queue_order (self):
        return arrivals.ARRV_BY_D

    def when_clear (self, Arrivals arrv):
        cdef Event last_evt = arrv.final_arrival_by_queue (self)
        return last_evt.d if last_evt else 0.0
        
    def previous_departing (self, Event e):
        return e.prev_byq # prev_byq actually is departure order
        
    cpdef recompute_service (self, Event e1):
        cdef Event e0, e2
        e0 = e1.prev_byq
        e2 = e1.next_byq
        
        if e0:
            e1.s = e1.d - max(e1.a, e0.d)
        else:
            e1.s = e1.d - e1.a
    
    cdef double departure_of_service  (self, Event e, double s):
        cdef Event byq = e.prev_byq
        cdef double d_prev = byq.d if byq else 0.0
        return s + max (e.a, d_prev)
    
    cpdef recompute_departure (self, Arrivals arrv, Event e1):
        cdef Event e_prev = e1.prev_byq
        cdef double d_prev = e_prev.d if e_prev else 0
        e1.set_departure (e1.s + max(e1.a, d_prev))
    
    def validate_service (self, Event e):
        cdef Ninq ninq
        d_prev = (e.prev_byq.d if e.prev_byq else 0)
        expected = e.d - max(e.a, d_prev)
        assert abs(e.s - expected) < 1e-10, "Invalid service time (expected %s was %s)\n  %s\n  %s" % (expected, e.s, e.prev_byq, e)
        if e.prev_byq:
            assert e.d >= e.prev_byq.d, "Out-of-order departure times for events PREV: %s CURR: %s\n  %s\n  %s" % (e.prev_byq.d, e.d, e.prev_byq, e)
        if e.next_byq:
            assert e.next_byq.d >= e.d, "Out-of-order departure times for events CURR: %s NEXT: %s\n  %s\n  %s" % (e.d, e.next_byq.d, e, e.next_byq)
            if False:  # comment out, for now
             if e.next_byq.a >= e.d:
                # or I must have emptied the queue, and A must be the minimum arriavl
                ninq = e.arrv.ninq_of_queue (self)
                assert 0 == ninq.N(e.d), "Inconsistent arrivals: next departure arrives after I depart (but queue not empty)\n%s\n%s" % (e, e.next_byq)
                arrivals = e.arrv.arrivals_in_range (self, e.d, e.next_byq.a + 1e-15)
                assert 1 == len(arrivals) and arrivals[0] is e.next_byq, "Inconsistent arrivals: queue empty, but first arrival didn't get served\nE0 %s\nE1 %s\nALL %s" % (e, e.next_byq, arrivals)

    def validate_caches (self, Arrivals arrv):
        cdef Event e0
        cdef Ninq ninq
        ninq = arrv.ninq_of_queue (self)
        evts_by_d = arrv.events_of_queue (self)
        for e0 in evts_by_d:
            w = e0.wait()
            assert 0 < ninq.N(e0.a) # should always count yourself
            if w < 1e-10:
                assert 1 == ninq.N(e0.a), "Event got started immediately, but others in q\n%s" % e0                
            else:
                assert 1 < ninq.N(e0.a), "Event entered empty queue but didn't start (w = %.15f)\n%s" % (w, e0)
            
    cdef Pwfun arrivalLik (self, Arrivals arrv, Event evt):
       raise NotImplementedError ("Use slice-only for this guy.")
        
    cdef Pwfun departureLik (self, Arrivals arrv, Event evt):
        # Return a Pwfun f_{service} (d_e - max{a_e, c_e})
        cdef Event pred = evt.prev_byq
        cdef double d_prev = pred.d if pred else 0.0
        return DepartureProposal (self.service, evt, max(evt.a, d_prev)).proposal()
           
    cdef diffListForArrival (self, Event evt, double a):
        cdef Event e_new = evt.dup_with_structure()
        e_new.a = a
        self.recompute_service (e_new)
        e_new.reset_commencement ()
        return [e_new]
        
    cdef diffListForDeparture (self, Event evt, double d):
        cdef Event old_prev = evt.prev_byq
        cdef Event old_next = evt.next_byq

        cdef Event new_next = old_next
        cdef Event new_prev = old_prev

#        print "DLFD ", d, evt 
#        print "  Old prev %s\n  Old next %s " % (old_prev, old_next) 

        while not (get_d (new_prev, 0) <= d):
            new_prev = new_prev.prev_byq
            new_next = new_next.prev_byq if new_next else evt
            if new_next is evt: new_next = new_next.prev_byq

#        print "New prev ", new_prev  

        while not (d <= get_d (new_next, INFINITY)):
            new_prev = new_prev.next_byq if new_prev else evt
            if new_prev is evt: new_prev = new_prev.next_byq
            new_next = new_next.next_byq

#        print "  New prev %s\n  New next %s " % (new_prev, new_next)

        lst = []
        cdef Event e0

        # special case: evt alone in its queue:
        if new_prev is None and new_next is None:
            e0 = evt.dup_with_structure()
            e0.d = d
            self.recompute_service (e0)
            return [e0]

        cdef Event e_first
        if old_prev is None and new_prev is None:
            e_first = evt
        elif old_prev != None and new_prev is None:
            e_first = new_next
        elif old_prev is None and new_prev != None:
            e_first = evt
        else: # both not none. take earliest
            e_first = old_prev if old_prev.d < new_prev.d else new_prev

        cdef Event e_last   # e_last is one past last event to copy
        if old_next is None or new_next is None:
            e_last = None
        else: # both not none. take earliest
            e_last = old_next if old_next.d > new_next.d else new_next
            e_last = e_last.next_byq 

#        print "  E_FIRST %s\n  E_LAST %s\n" % (e_first, e_last)

        # collect all in lst
        cdef Event e_curr = e_first
        while e_curr is not e_last:
#            print "Duping ", e_curr
            e0 = e_curr.dup_with_structure()
            if e_curr == evt: e0.d = d
            e_curr = e_curr.next_byq
            lst.append (e0)
            
        lst = [ (evt.d, evt) for evt in lst ]
        lst.sort()
#        print "Sorted ", lst
        lst = [ tup[1] for tup in lst ]

        # set next pointers
        for i in range(len(lst)-1):
            e0 = lst[i]
            e1 = lst[i+1]
            e0.set_next_by_queue (e1)
        (<Event>lst[0]).prev_byq = e_first.prev_byq
        (<Event>lst[-1]).next_byq = e_last

        # set service
        for i in range(len(lst)):
            self.recompute_service (lst[i])
            lst[i].reset_commencement ()

        # done!
#        print_diff_list (lst)
        return lst

    cpdef double likelihoodDelta (self, Arrivals arrv, diffList) except *:
        cdef Event e_new, e_old
        cdef Ninq ninq0
        cdef Overlay ninq1
        cdef double lik = Queue.likelihoodDelta (self, arrv, diffList)
        if lik == -numpy.inf: return lik
#        print "/////// LD", lik, self.as_yaml(), "\n DL", diffList #GGG
        ninq0 = arrv.ninq_of_queue (self)
        ninq1 = Overlay (ninq0)
        for e_new in diffList:
            if e_new.q == self:
                e_old = arrv.event (e_new.eid)
                if e_new.d != e_old.d: ninq1.move_departure (e_old.d, e_new.d)
                if e_new.a != e_old.a: ninq1.move_arrival (e_old.a, e_new.a)
        # collect all events who could have seen their nq change under them
        # N.B. techically, a depature change can't do this (although I'm not really enforcing the constraint)
        nq_change = ninqueue.EventSet()
        for e_new in diffList:
            if e_new.q == self:
                nq_change.add (e_new)
                e_old = arrv.event (e_new.eid)
                a0 = e_old.a
                a1 = e_new.a
                if a0 != a1:
                    lst = arrv.c2e_range (self, min (a0,a1), max(a0, a1))
                    for x in lst:
                        if not nq_change.contains(x):
                            nq_change.add (x)
#        print "DL LEN = ", len(diffList) 
#        print "EVTSET LEN = ", len(nq_change.items()) 
        # for each event, contribution to log-likelihood is -log(Nq_e)
        for e_new in nq_change.items():
            if e_new.q == self:
                # first check the rss constraint
                if e_new.wait() < 1e-10:                
                    Na = ninq1.N(e_new.a)
                    if Na > 1:
                        print "likelihoodDelta: RSS constraint violated [w %.15f but N(a) = %d]\n %s" % (e_new.wait(), Na, e_new)
                        return -INFINITY                
                e_old = arrv.event (e_new.eid)
                n0 = ninq0.N (e_old.commencement())
                n1 = ninq1.N (e_new.commencement())  
                if e_old.service() != 0:
                    assert n0 > 0, "Huh? N0 == 0 at %s\n%s\nE.A %.15f  E.C %.15f\nNINQ %s\n%s" % (e_old, diffList, e_old.a, e_old.commencement(), ninq0, arrv)
                    lik += numpy.log (n0)
                if e_new.service() != 0:
                    assert n1 > 0, "Huh? N1 == 0 at %s\nE.A %.15f E.C %.15f\n%s\nNINQ %s\n%s" % (e_new, e_new.a, e_new.commencement(), diffList, ninq1, arrv)
#                print "DELTA", n0, n1, e_old
                    lik -= numpy.log (n1)
#        print "\\\\\\\\\\\\\\\\ LD", lik
        return lik

    cdef double allEventLik (self, Arrivals arrv) except *:
       cdef double result = Queue.allEventLik (self, arrv)
       cdef Event evt
       cdef Ninq ninq
       ninq = arrv.ninq_of_queue (self)
       for evt in arrv.events_of_queue (self):
           n = ninq.N (evt.commencement ())
           if evt.s != 0:
               assert n > 0, "Huh? None in queue for %s\nE.A %.15f E.C %.15f\nNINQ %s\n%s" % (evt, evt.a, evt.commencement(), ninq, arrv)
               result -= numpy.log (n)
           # check rss constraints
           if evt.wait() < 1e-10:
               Na = ninq.N(evt.a)
               if Na > 1:
                   print "allEventLik: RSS constraint violated [w %.15f but N(a) = %d]\n %s" % (evt.wait(), Na, evt)
                   return -INFINITY
       return result

    def select_job_for_service (self, evt_list):
        return sampling.random_elt (evt_list)


# Processor sharing queue

cdef class QueuePS (Queue):

    def __repr__ (self):
        return "<QUEUE TYPE='PS' NAME='%s'>\n<SERVICE>%s</SERVICE>\n</QUEUE>" % (self.name, self.service)

    def as_yaml (self):
        return "  - { name: %s, type: PS, service: %s }" % (self.name, self.service.as_yaml())


    def when_clear (self, Arrivals arrv):
        return 0.0
        
    def previous_departing (self, Event e):
        raise NotImplementedError()

    # This can be called with the "real" ninq for the Arrivals,
    #  or with an overlay
    def recompute_service_internal (self, Event e1, Overlay ninq):
        cdef double t0,t1
        cdef int n
        cdef double s = 0
        cdef OverlayIterator iterator
#        print "RSI", e1

        if e1.d - e1.a < 1e-15:
            e1.s = (e1.d - e1.a) / (1+ninq.N(e1.a))
        else:
            iterator = ninq._interval_iterator (e1.a, e1.d)
            while iterator.has_next():
                t0 = iterator.T0()
                t1 = iterator.T1()
                n = iterator.N()
 #               print "T0 %.5f T1 %.5f N %d" % (t0, t1, n)
                if n > 0:
                    s += (t1 - t0) / n 
                else:
                    assert (t1-t0) < 1e-15, "N < 0\n%s\nEVT: %s\n%s" % (n, e1, ninq.knots())
                iterator.advance()
            e1.s = s

    cpdef recompute_service (self, Event e1):
        cdef Ninq ninq = e1.arrv.ninq_of_queue (self)
        self.recompute_service_internal (e1, Overlay(ninq))

    cdef double departure_of_service  (self, Event e, double s):
        cdef double remaining = s
        cdef double t0 = e.a
        cdef double d = e.a
        cdef double diff
        cdef int n0, n_last

        cdef Ninq ninq0 = e.arrv.ninq_of_queue (self)
        cdef Overlay ninq1 = Overlay (ninq0)
        cdef OverlayIterator iterator
        ninq1.move_departure (e.d, e.a)

        n_last = ninq1.N (t0)

 #       print "**D_O_S", e.eid, ninq0, ninq1
            
        while remaining > 0:
            iterator = ninq1.interval_iterator (t0, t0 + (1+n_last) * remaining + 0.01) # adding slack can't hurt
            while iterator.has_next():
                t0 = iterator.T0()
                t1 = iterator.T1()
                n0 = iterator.N()
                if n0 >= 0:
                    diff = min (t1 - t0, remaining * (1+n0))
                    remaining -= diff / (1+n0)
                    d += diff
                else:
                    assert (t1-t0) < 1e-15, "N0 < 0\n%s" % ninq1
                iterator.advance()            
            n_last = n0
            t0 = t1
        
        return d


    # recomputes all departure times for a change in e1's service
    cpdef recompute_departure (self, Arrivals arrv, Event e0):
        cdef double d0,d1
        cdef Event evt, e
        changed = dict()
        changed[e0.eid] = (e0, e0.d)

#        print "REC_D", e0 
#        print arrv.inq[e0.qid]
        
        while changed:
            eid = changed.keys()[0]
            evt,d_old = changed[eid]
            del changed[eid]
#            print "   UPD", evt 

            evt.set_departure (self.departure_of_service (evt, evt.s))
            d0 = min (d_old, evt.d)
            d1 = max (d_old, evt.d)
            if d1-d0 > 1e-15:
                to_update = arrv.evts_in_range (self, d0, d1)
                for e in to_update:
                    if e is not evt:
                        changed[e.eid] = (e, e.d)

    def validate_service (self, Event e):
        expected = self.departure_of_service (e, e.s)
        assert abs (expected - e.d) < 1e-8, "Service mismatch (expected departure %.16f)\n%s" % (expected, e.dump())

    cdef Pwfun arrivalLik (self, Arrivals arrv, Event evt):
       raise NotImplementedError ("Use slice-only for this guy.")
        
    cdef Pwfun departureLik (self, Arrivals arrv, Event evt):
        # Return a Pwfun f_{service} (d_e - max{a_e, c_e})
        cdef Event pred = evt.prev_byq
        cdef double d_prev = pred.d if pred else 0.0
        return DepartureProposal (self.service, evt, max(evt.a, d_prev)).proposal()
           
    cdef diffListForArrival (self, Event evt, double a):
        cdef double a0 = min (a, evt.a) - 1e-15
        cdef double a1 = max (a, evt.a) + 1e-15
        evts = evt.arrv.evts_in_range (self, a0, a1)
        
        ninq0 = evt.arrv.ninq_of_queue (self)
        ninq1 = Overlay (ninq0)
        ninq1.move_arrival (evt.a, a)

#        print "DLFA", evts
        lst = self.adjust_all_service (evts, ninq1, evt)
        cdef Event e0 = evt.dup_with_structure ()
        e0.a = a
        self.recompute_service_internal (e0, ninq1)
        lst.append (e0)

        return lst
        
    cdef diffListForDeparture (self, Event evt, double d):
        cdef double d0 = min (d, evt.d) - 1e-15
        cdef double d1 = max (d, evt.d) + 1e-15
        evts = evt.arrv.evts_in_range (self, d0, d1)

        ninq0 = evt.arrv.ninq_of_queue (self)
        ninq1 = Overlay (ninq0)
        ninq1.move_departure (evt.d, d)

#        print "DLFD", evts
        lst = self.adjust_all_service (evts, ninq1, evt)
        cdef Event e0 = evt.dup_with_structure ()
        e0.d = d
        self.recompute_service_internal (e0, ninq1)
        lst.append (e0)

        return lst
            

    def adjust_all_service (self, evtl, Overlay ninq, exclude):
        cdef Event e0, e1
        result = []
        for e0 in evtl:
            if e0 is exclude: continue
#            print "AAS", e0.arrv.evts_in_range (self, e0.a, e0.d) 
            e1 = e0.dup_with_structure()
            self.recompute_service_internal (e1, ninq)
            result.append (e1)
        return result

    def select_job_for_service (self, evt_list): 
        return evt_list[0]

    cdef double allEventLik (self, Arrivals arrv) except *: 
        cdef double result = Queue.allEventLik (self, arrv)
        return result
        cdef Event evt
        cdef Ninq ninq = arrv.ninq_of_queue (self)
        for evt in arrv.events_of_queue (self):
            result -= log (ninq.N (evt.d) + 1)
#            print evt.eid, ninq.N (evt.d - EPS)
        return result

    cdef Overlay overlayDiffList (self, Arrivals arrv, Ninq ninq, diffList):
        cdef Overlay ninq1 = Overlay(ninq)
        cdef Event e_new
        for e_new in diffList:
            if e_new.q == self:
                e_old = arrv.event (e_new.eid)
                if e_new.d != e_old.d: ninq1.move_departure (e_old.d, e_new.d)
                if e_new.a != e_old.a: ninq1.move_arrival (e_old.a, e_new.a)
        return ninq1
        
    cpdef double likelihoodDelta (self, Arrivals arrv, diffList) except *:
        cdef double lik = Queue.likelihoodDelta (self, arrv, diffList)
        cdef Event e_new, e_old
        cdef Ninq ninq = arrv.ninq_of_queue (self)
        cdef Overlay overlay = self.overlayDiffList (arrv, ninq, diffList)
        for e_new in diffList:
            if e_new.q == self:
                e_old = arrv.event (e_new.eid)
                N_old = ninq.N (e_old.d) + 1
                N_new = overlay.N (e_new.d) + 1
                lik = lik - log(N_new) + log(N_old)
#                print "%d %d %d" % (e_new.eid, N_old, N_new) 
        return lik


    
# MH proposal switching


cdef Pwfun _pair_proposal (Arrivals arrv, Event e0, Event e1): 
    cdef Pwfun d_fn, a_fn, product

    if proposal == TWOS_PROPOSAL:
        d_fn = e0.q.departureLik (arrv, e0)
        a_fn = e1.q.arrivalLik (arrv, e1)
        
        if DEBUG:
            print "a_fn:  ", a_fn
            print "d_fn:  ", d_fn
            a_fn.validate()
            d_fn.validate()
            
        product = pwfun.add (a_fn, d_fn)
        return product

    elif proposal == UNIFORM_PROPOSAL:
        return netutils.uniform_proposal (e0.a, e1.d)
    elif proposal == BAR_PROPOSAL:
        return netutils.bar_pair_proposal (arrv, e0, e1)
    elif proposal == ZOID_PROPOSAL:
        return netutils.zoid_pair_proposal (arrv, e0, e1)
    else:
        print "Could not understand propsal type %d" % proposal
        return None

cdef Pwfun _departure_proposal (Arrivals arrv, Event e0):
    if proposal == TWOS_PROPOSAL:
        return (<Queue> e0.q).departureLik (arrv, e0)
    elif proposal == UNIFORM_PROPOSAL:
        return netutils.uniform_proposal (e0.a, numpy.inf)
    elif proposal == BAR_PROPOSAL:
        return netutils.bar_final_proposal (arrv, e0)
    elif proposal == ZOID_PROPOSAL:
        return netutils.zoid_final_proposal (arrv, e0)
    else:
        print "Could not understand propsal type %d" % proposal
        return None
        raise Exception

def pair_proposal (arrv, e0, e1): 
    return _pair_proposal (arrv, e0, e1)

def departure_proposal (arrv, e0): 
    return _departure_proposal (arrv, e0)

# Factor class

# This is sort of a general graphical model framework, but using
#  the Arrivals datastructures, and making the assumption that all
#  variables are replicated for each event

cdef class FactorTemplate:
    cdef object estimate (self, evtl): raise NotImplementedError()
    cdef object serviceLikelihood (self, Event evt): raise NotImplementedError()
    cdef object serviceDerivative (self, Event evt): raise NotImplementedError()
    cdef void resample_auxiliaries (self, Event evt): raise NotImplementedError()
    cdef void initialize_auxiliaries (self, Event evt): raise NotImplementedError()
    def numParameters (self): raise NotImplementedError()
    def sample_service (self): raise NotImplementedError()

    property parameters:
        def __get__ (self): return self.getParameters ()
        def __set__ (self, v): self.setParameters (v)
    
# A factor object that represents a univariate service distribution
cdef class DistributionTemplate (FactorTemplate):
    cdef Distribution f
    cdef double zero_prob

    def __init__ (self, f): 
        self.f = f
        self.zero_prob = 0.01

    def __repr__(self):
        return "<DISTRIBUTION ZERO_PROB=\"%s\">%s</DISTRIBUTION>" % (self.zero_prob, self.f)

    def sample_parameters (self, evtl):
        s_obs = [ evt.s for evt in evtl if evt.s > 0 ]
        self.parameters = self.f.sample_parameters (s_obs)

    def parameter_kernel (self, theta0, theta1, evtl):
        s_obs = [ evt.s for evt in evtl if evt.s > 0 ]
        return self.f.parameter_kernel (theta0, theta1, s_obs)
 
    cdef object estimate (self, evtl): 
        cdef Event evt
        cdef Arrivals arrn
        
        cdef int N0 = 0
        for evt in evtl:
            if evt.s == 0:
                N0 += 1
        self.zero_prob = (<double>N0) / len(evtl)

        s_obs = [ evt.s for evt in evtl if evt.s > 0 ]
        return self.f.estimate (s_obs)

    cdef void resample_auxiliaries (self, Event evt): pass
    
    cdef void initialize_auxiliaries (self, Event evt): pass

    def sample_service (self): 
        return self.f._sample(1)[0], arrivals.Vars ({})
    
    # not currently used
#    def _service_loglik (self, s):
#        if s == 0:
#            return log(self.zero_prob)
#        else:
#            return log(1-self.zero_prob) + self.f.lpdf(s)

    cdef object serviceLikelihood (self, Event evt):
        return self.f.lpdf
    
    cdef object serviceDerivative (self, Event evt):
        return self.f.dx_lpdf
    
    def numParameters (self): return self.f.numParameters()

    def getParameters (self): 
        v = self.f.parameters[:]
        return v

    def setParameters (self, v): 
        self.f.parameters = v

    def mean (self): return self.f.mean()
    def std (self): return self.f.std()
        
    # maybe hacky
    def as_yaml (self):
        return dist2yaml (self.f)
    
cdef class MixtureTemplate (FactorTemplate):
    
    cdef object var_name
    cdef object mixing
    cdef object dists
    cdef double zero_prob

    def __init__(self, var_name, mixing, dists):
        self.var_name = var_name
        self.mixing = mixing[:]
        self.dists = dists[:]
        self.zero_prob = 0

    def __repr__(self):
        result = "<MIXTURE>\n  probs: %s\n  zero_prob: %s\n" % (self.mixing, self.zero_prob)
        for d in self.dists:
            result += "  "
            result += str(d)
            result += "\n"
        result += "</MIXTURE>"
        return result

    cdef void initialize_auxiliaries (self, Event evt):
       p = [1.0] * len(self.mixing)
       cdef int i = _roll_die_unnorm (p)
       if i >= len(p):
          print "ERROR: initialize_auxiliaries: sampled %d but probs were %s" % (i, p)
       evt.set_variable (self.var_name, i)

    cdef void resample_auxiliaries (self, Event evt):
        p = self.mixing[:]
        cdef int i,j
        for i from 0 <= i < len(p):
            p[i] *= exp(self.dists[i].lpdf (evt.s))
        j = _roll_die_unnorm (p)
        if j >= len(p):
           print "ERROR: initialize_auxiliaries: sampled %d but probs were %s\n %s" % (j, p, evt)
        evt.set_variable (self.var_name, j)
#        print "==== resample_auxiliaries"
#        print evt.dump()
#        print p
        
    # estimates mixture components & likelihoods
    def collect_sobs (self, evtl):
        cmp_counts = [1] * len(self.mixing)  # small dirichlet prior
        cdef int N = len(self.mixing)
        cdef int N0 = 0
        
        s_obs = []
        for i in xrange(len(self.mixing)): s_obs.append([])
        for evt in evtl:
            component = evt.variable (self.var_name)
            cmp_counts[component] += 1
            if evt.s > 0:
                s_obs[component].append (evt.s)
            else: N0 += 1
            N += 1
        return (s_obs, cmp_counts, N, N0)

    cdef object estimate (self, evtl):

        s_obs, cmp_counts, N, N0 = self.collect_sobs (evtl)

        # estimate zero prob
        self.zero_prob = (<double>N0) / N

        # estimate inner distributions
        for i in xrange(len(self.mixing)):
            self.dists[i].estimate (s_obs[i])
    
        # estimate mixing proportions
        for i in xrange(len(self.mixing)): 
            self.mixing[i] = cmp_counts[i] / float(N)
        
        return self.parameters

    def sample_parameters (self, evtl):
        s_obs, cmp_counts, N, N0 = self.collect_sobs (evtl)
        for i in xrange(len(self.mixing)): 
            self.mixing[i] = cmp_counts[i] / float(N)
        for i in xrange(len(self.mixing)):
            params = self.dists[i].sample_parameters (s_obs[i])
            self.dists[i].parameters = params
        return self.parameters

    def parameter_kernel (self, theta0, theta1, evtl):
        raise NotImplementedError()

    cdef object serviceLikelihood (self, Event evt):
        i = evt.variable(self.var_name)
        return self.dists[i].lpdf

    cdef object serviceDerivative (self, Event evt):
        i = evt.variable(self.var_name)
        return self.dists[i].dx_lpdf

    def sample_service (self):
        i = _roll_die_unnorm (self.mixing)
        aux = arrivals.Vars({ self.var_name: i })
        return self.dists[i].sample(1)[0], aux
        
    def numParameters (self):
        n = len(self.mixing)
        for d in self.dists: n += len(d.parameters)
        return n
        
    def getParameters (self):
        l = []
        l.extend (self.mixing)
        for d in self.dists: l.extend (d.parameters)
        return l
        
    def setParameters (self, v):
        j = len(self.mixing)
        self.mixing = v[0:j]
        i = j

        for d in self.dists:
            j = i + len(d.parameters)
            d.parameters = v[i:j]
            i = j

    def as_yaml (self):
        ret = "[ MIX, "
        for i in xrange(len(self.mixing)): 
            ret += "%.5f, %s, " % (self.mixing[i], dist2yaml(self.dists[i]))
        ret += "]"
        return ret
                
    def mean (self):
        sub_means = [ w*f.mean() for w,f in zip(self.mixing,self.dists) ]
        return sum(sub_means)

    def std (self):
        sub_stds = [ f.std() for f in self.dists ]
        weighted = [ w*std*std for w,std in zip(self.mixing, sub_stds) ]
        return math.sqrt (sum (weighted))
    
cdef object dbl2array (double* d, int N):
    cdef int i
    l = []
    for i from 0 <= i < N:
        l.append (d[i])
    return l

cdef double famin (double* d, int len, int *out_argmin):
    cdef int i
    cdef int mi
    cdef double m = INFINITY
    for i from 0 <= i < len:
        if d[i] < m: 
            m = d[i]
            mi = i
    out_argmin[0] = mi
    return m
 
def allequal (l1, l2):
    if len(l1) != len(l2): return False
    for x1, x2 in zip(l1, l2):
        if l1 != l2: return False
    return True

cdef double get_d (Event evt, double default):
    if evt is None:
        return default
    else:
        return evt.d
    
def print_diff_list (dl):
    cdef Event e
    for e in dl:
        print e,
        print (e.prev_byq.eid if e.prev_byq else -1),
        print (e.next_byq.eid if e.next_byq else -1)

def get_arrival (evt): return evt.a

def dist2yaml (dist):
    if isinstance (dist, distributions.Exponential):
        ltr = "M"
    elif isinstance (dist, distributions.Gamma):
        ltr = "G"
    elif isinstance (dist, distributions.LogNormal):
        ltr = "LN"
    else:
        raise Exception ("Can't interpret %s" % dist)
    paramStr = ",".join (map (str, dist.parameters))
    return "[ %s, %s ]" % (ltr, paramStr)
