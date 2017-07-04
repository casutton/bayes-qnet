##  Data structure for managing Arrivals.
##   This needs to be really fast


cdef extern from "math.h":
    enum: INFINITY

cdef enum:
    DEBUG = 0
    DEBUG_RTP = 0
    DEBUG_INITIAL = 0

cdef enum:
    SEC_ARRIVAL = 0
    SEC_CLEAR = 1
    SEC_DEPARTURE = 2

cdef extern from "math.h":
    double exp(double x)
    double log(double x)
    
from heapq import heappush, heappop
import bisect

import numpy
from numpy import random

import yaml

from distributions cimport Distribution
from queues cimport FactorTemplate
from queues cimport Queue, _departure_proposal, _pair_proposal
from pwfun cimport Pwfun
from arrivals cimport Arrivals
from arrivals cimport Event

import queues, distributions, arrivals, qstats

from sampling cimport _slice_pwfun, _slice, rand, uniform, exponential
from misc cimport _integrate

import netutils
import sampling
import pwfun
import stupidlp
import hmm
import warnings
import misc
          
import sys, traceback      
import math # temporary

from scipy import integrate

# configuration: need to expose in a better way
cdef enum SamplingType:
    ARS_PWFUN = 0
    SLICE_PWFUN = 1
    SLICE = 2
    STUPID_REJECT = 3

cdef int doUseRtp = 0
cdef SamplingType sampler = SLICE

cdef int BRUTE_FORCE_RATIO = 0

cdef object initializer

cdef int is_slice_sorted = 0

def set_sampler (str):
    global sampler
    if str == "SLICE": sampler = SLICE
    elif str == "SLICE_PWFUN": sampler = SLICE_PWFUN
    elif str == "ARS_PWFUN": sampler = ARS_PWFUN
    elif str == "STUPID_REJECT": sampler = STUPID_REJECT
    else: raise Exception("Could not understand one_d_sampler: %s" % str)

def set_use_rtp (doitp):
    global doUseRtp
    doUseRtp = doitp

def set_initialization_type (type):
    global initializer
    type = type.upper() if isinstance(type, str) else type
    if callable (type):        
        initializer = type
    elif type == "DFN1":
        initializer = gibbs_initialize_via_dfn1
    elif type == "DFN2":
        initializer = gibbs_initialize_via_dfn2
    elif type == "LP":
        initializer = gibbs_initialize_via_lp
    elif type == "MINILP":
        initializer = gibbs_initialize_via_minilp
    elif type == "PS":
        initializer = gibbs_initialize_for_ps
    elif type == "NONE":
        initializer = null_initializer
    else:
        raise Exception ("Cannot understand type %s" % type)
    print "New QNET initializer = ", initializer

def get_initializer (): return initializer

def set_slice_sorted (val):
    global is_slice_sorted
    is_slice_sorted = val
    "QNET slice_sorted = ", val

# end configuration

# stastics

cdef int nreject = 0
cdef int nevents = 0
cdef int nrtp_reject = 0
cdef int nrtp_calls = 0
cdef int diff_len = 0
cdef int num_diff = 0

def reset_stats():
    global nreject, nevents, nrtp_calls, nrtp_reject, diff_len
    nreject = 0
    nevents = 0
    nrtp_reject = 0
    nrtp_calls = 0
    diff_len = 0
    num_diff = 0
    
def stats():
    return {'nreject':nreject, 'nevents':nevents, 'nrtp_calls':nrtp_calls, 'nrtp_reject':nrtp_reject, 'diff_len':diff_len, 'num_diff':num_diff}
    
reset_stats()

# end stats

INITIAL_QID = 0
INITIAL_STATE = 0

cdef class GGkGibbs
 
cdef class Qnet:
        
    def __init__ (self, queues, templates, universe, qname2id, sname2id, fsm):
        self.queues = queues
        self.templates = templates
        self.universe = universe
        self.fsm = fsm
        self.qname2id = qname2id
        self.sname2id = sname2id
        self.eidx = -1

    def as_yaml (self):
        ret = "states:\n"
        for si in range(self.fsm.num_states()):
            sname = self.sid_to_name (si)
            ret += "  - name: %s\n" % sname
            qnames = [ q.name for q in self.queues_of_state (si) ]
            ret += "    queues: [ %s ]\n" % ",".join (qnames)
            successors = [ self.sid_to_name (si_next) for si_next in self.fsm.successors (si) ]
            if len(successors) > 0:
                ret += "    successors: [ %s ]\n" % ",".join(successors)
        ret += "queues:\n"
        for qi in range(len(self.queues)):
            ret += self.queues[qi].as_yaml ()
            ret += "\n"
        return ret

            

    def as_xml (self):
        ret = "<QNET>\n  <QUEUES>\n"
        for q in self.queues: ret = ret + "  " + str(q) + "\n"
        ret = ret + "\n  </QUEUES>\n</QNET>"
        return ret


    # forward sampling code. complex but extremely general
    #  SEC_* constants defined at top of file

    # samples assuming that an infinite number of customers queue up according
    #  to the initial state.
    def sample (self, n):
        """Generates n jobs from this qnet."""
        arrv = Arrivals (self)
        
        heap = []
        for i in xrange(n):
            s_initial = self.fsm.initial_state()
            heappush (heap, (0.0, i, s_initial, SEC_ARRIVAL, None))

        return self._sample_discrete_events (heap, arrv, True)

    def sample_by_time (self, T):
        """Simulates the system from time 0...T.  Arrivals must be poisson."""
        q0 = self.queues[0]
        delta_bar = q0.service.parameters[0]
        if not q0.is_markov():
            raise Exception ("ERROR: sample_by_time works only for M/G/xxx queues")
        if not len(q0.service.parameters) == 1:
            raise Exception ("Huh?")

        arrv = Arrivals (self)
        n = random.poisson (T / delta_bar)        
        heap = []

        allD = random.uniform(0, T, n)
        allD.sort()

        sid = self.fsm.initial_state()
        qid = self.fsm.sample_one_obs (sid)
        q = self.queues[qid]
        
        d_prev = 0
        for i in xrange(n):
            d = allD[i]
            s = d - d_prev
            tid = eid = i
            evt = Event (eid, tid, qid, q, 0.0, s, d, sid, obs_a=1, obs_d=1)
            arrv.append (evt)
            heappush (heap, (evt.d, tid, None, SEC_DEPARTURE, evt))
            d_prev = d

        return self._sample_discrete_events (heap, arrv, True)

    # discrete event simulation, private
    def _sample_discrete_events (self, heap, Arrivals arrv, obs):
        job_lists = [ list() for i in range(len(self.queues)) ]

        cdef Queue q
        cdef Event evt

        eid = arrv.max_eidx() + 1

        while heap:
            tup = heappop(heap)
            event_type = tup[3]
#            print "TYPE %d: @ %.5f  TID: %d QID: %s STATE: %s" % (tup[3], tup[0], tup[1], tup[4], tup[2])

            if event_type == SEC_ARRIVAL:
                a_i, tid, sid, event_type, arg = tup
                qid = self.fsm.sample_one_obs (sid)
                q = self.queues[qid]
                evt = Event (eid, tid, qid, q, a_i, 0.0, a_i, sid, obs_a=obs, obs_d=obs)
                eid += 1
                bisect.insort (job_lists[qid], evt)
                if len(job_lists[qid]) == 1: # wake up queue, if necessary
                    t = q.when_clear (arrv)
#                    print "Pushing ", (t, tid, 0, SEC_CLEAR, qid)
                    heappush (heap, (t, tid, 0, SEC_CLEAR, qid))

            elif event_type == SEC_CLEAR:
                qid = tup[4]
                if job_lists[qid]:
                    q = self.queues[qid]

                    # select job to serve, remove from queue
                    evt = q.select_job_for_service (job_lists[qid])
                    tid = evt.tid
                    job_lists[qid].remove (evt)

                    #  add served event to Arrivals
                    s_i, aux = q.sample_service()
                    evt.s = s_i
                    evt.auxiliaries = aux
                    arrv.append (evt)
                    q.recompute_departure (arrv, evt)
#                    print "   Added event: ", evt
#                    arrv.validate()
                    
                    # create arrival event
                    heappush (heap, (evt.d, tid, None, SEC_DEPARTURE, evt))

                    # create next clear event
                    t = q.when_clear (arrv)
                    heappush (heap, (t, 0, 0, SEC_CLEAR, qid))


            elif event_type == SEC_DEPARTURE:
                t, tid, foo1, foo2, evt = tup
                # Below check necessary b/c of PS queues
                if evt.d == t:
                    s_next = self.fsm.sample_state (evt.state)
                    if s_next: 
                        heappush (heap, (t, tid, s_next, SEC_ARRIVAL, None))
                else:
                    # oops, job isn't really finished yet
                    #  (other concurrent jobs came into the PS queue)
                    heappush (heap, (evt.d, tid, None, SEC_DEPARTURE, evt))


            else: # event_type == something else
                raise Exception ("Could not understand event type %d" % event_type)

        arrv.regenerate_all_caches()

        return arrv
        
    def num_queues (self): return len(self.queues)
    def num_states (self): return self.fsm.num_states()
    def get_fsm (self): return self.fsm
    
    def sid_by_name (self, sname):
        return self.sname2id.get (sname, None)
        
    def sid_to_name (self, sid):
        for k,v in self.sname2id.iteritems():
            if v == sid:
                return k
        raise KeyError ("No state with id %s" % sid)

    def queues_of_state (self, sid):
        qids = self.fsm.possible_observations (sid)
        return [ self.queues[qi] for qi in qids ]

    # inefficient
    def qid_of_queue (self, queue):
        if queue in self.queues:
            return self.queues.index (queue)
        else:
            raise Exception ("Can't find %s in net %s" % (queue, self))
        
    def queue_by_name (self, qname):
        qid = self.qname2id.get (qname, None)
        if qid != None:
            return self.queues[qid]
        else:
            return None
    
    def all_queues (self): return self.queues[:]

    def all_templates (self): return self.templates[:]

    # Gibbs sampling
    
    def gibbs_resample (self, arrv, burn, num_iter, report_fn=None, return_arrv=True, sample_final=True):

        cdef int gi, i
        cdef Event evt, e_next
        cdef FactorTemplate factor
        
        cdef int nt05 = len(arrv.all_tasks()) / 20
        if nt05 == 0: nt05 = 1

        all_arrv = []

        for gi from 0 <= gi < num_iter + burn:
            
            # resample latent variables for all factors
            for evt in arrv:
                evt.q.resample_auxiliaries (evt)
                    
            for evt in arrv:
#                if not arrv.is_final(evt): continue
#                if not evt.qid != 1: continue

                if not evt.obs_d:
                    if arrv.is_final (evt):
                        if sample_final: 
                            if DEBUG: print "Resampling final event:\n  %s\n  %s\n  %s" % (evt.prev_byq, evt, evt.next_byq)
                            self.gibbs_resample_final (arrv, evt)
                    else:
                        e_next = evt.next_byt
                        if not e_next.obs_a:
                            if DEBUG: print "Resampling:\n", format_six (evt)
                            self.gibbs_resample_pair (arrv, evt, e_next)                                            
                if DEBUG: arrv.validate() 
            
            if doUseRtp:
                for i from 0 <= i < nt05:
                    arrv = self.remove_task_proposal (arrv)

            if burn <= gi:
                if report_fn:
                    report_fn (self, arrv, gi)
                if return_arrv:
                    a2 = arrv.duplicate()
                    all_arrv.append (a2)

        if not return_arrv:
            all_arrv = [arrv]

        return all_arrv                



    cdef int gibbs_resample_final (self, Arrivals arrv, Event evt) except -1:
        cdef Pwfun fn = _departure_proposal (arrv, evt)
        
        if DEBUG: print "Final d_fn: ", fn
        if fn.U() - fn.L() < 1e-10: return 0  # don't sample small intervals
        
        cdef double new_s, d_new

        cdef int tries = 0
        cdef double x_initial = evt.d
        while (<double> fn(x_initial)) == -INFINITY:
            new_s = exponential (1.0)
            x_initial = (<Queue> evt.q).departure_of_service (evt, new_s)
            tries += 1
            if tries > 25: 
                if DEBUG: print "Can't find valid initialization"
                return 0 # can't find a valid initialize
        
        try:
            d_new = call_sampler (fn, x_initial)            
        except:
            print "Error resampling event %s" % evt
            print format_six (evt)
            raise

        cdef Event e_new, e_old
        cdef double ratio, u
        global nreject, nevents 
        nevents += 1

        if DEBUG: print "was", evt.d, "sampled", d_new
        diff_list = (<Queue> evt.q).diffListForDeparture (evt, d_new)

        global diff_len, num_diff
        num_diff += 1
        diff_len += len (diff_list)

        ratio = math.exp(fn(evt.d) - fn(d_new))
        for e_new in diff_list:
            if e_new.s < 0:
                ratio = 0
                break
            e_old = arrv.event (e_new.eid)
            ratio *= exp(e_new.q.service_lpdf (e_new) - e_old.q.service_lpdf (e_old))

        u = rand()
        if DEBUG: print "Ratio: ", ratio

#        print "diff_list\n", display_diff_list (arrv, evt, diff_list) 
        
        if u < ratio:
            arrv.applyDiffList (diff_list)
        else:
            if DEBUG: print "Rejected:\n", display_diff_list (arrv, evt, diff_list)
            nreject += 1

        return 0
        
        
    # Resamples departure of event E0 == arrival of event E1
    cdef int gibbs_resample_pair (self, Arrivals arrv, Event e0, Event e1) except -1:
        cdef Pwfun product = _pair_proposal (arrv, e0, e1)

        if product.U() - product.L() < 1e-10 or e1.d - e0.a < 1e-10: return 0  # don't sample small intervals

        cdef int tries = 0
        cdef double x_initial = e0.d
        while (<double>product(x_initial)) == -INFINITY:
            if DEBUG: print "Using strange initialization."
            x_initial = uniform (e0.a, e1.d)
            tries += 1            
            if tries > 25:
                if DEBUG: print "Can't find valid initialization"                
                return 0
        
        try:        
            d_new = call_sampler (product, x_initial)
        except:
            print "Error resampling event %s" % e0
            print format_six (e0)
            L = product.L()
            U = product.U() - 1e-8
            mid = (L+U)/2
            raise
            
        if DEBUG: print "product: ", product
        if DEBUG: print "was", e0.d, "sampled", d_new

        cdef Event e_new, e_old
        cdef Queue q0 = e0.q
        cdef Queue q1 = e1.q
        cdef double ratio, u
        global nreject, nevents
        nevents += 1
        
        diff_list = q0.diffListForDeparture (e0, d_new)
        diff_list.extend (q1.diffListForArrival (e1, d_new))
                
        global diff_len, num_diff
        num_diff += 1
        diff_len += len (diff_list)

        ratio = math.exp(product(e0.d) - product(d_new))
        for e_new in diff_list:
            if e_new.s < 0: 
                ratio = 0
                break
            e_old = arrv.event (e_new.eid)
            ratio_this = exp(e_new.q.service_lpdf (e_new) - e_old.q.service_lpdf (e_old))
            ratio *= ratio_this
#            print "%s\n%s\n%.5f %.5f\n\n" % (e_new, e_old, ratio_this, ratio) 

        u = rand()
        if DEBUG: print "Ratio", ratio

        if u < ratio:
#            print "diff list: ", display_diff_list (arrv, e0, diff_list) 
            arrv.applyDiffList (diff_list)
#             lpOld = self.log_prob (arrv)
#             self.applyDiffList (arrv, diff_list)
#             lpNew = self.log_prob (arrv)
#             print "Ratio real: ", exp(lpNew - lpOld)*ratio0
        else:
            if DEBUG: print "Rejected:\n", display_diff_list (arrv, e0, diff_list)
            nreject += 1
        
        return 0
                
    # For debugging / unit tests
    def gibbs_for_departure (self, a, evt):
        cdef Arrivals arrv = a
        cdef Event e0 = evt
        cdef Event e1 = e0.next_byt
        if e1:
            a_fn = e0.q.departureLik (arrv, e0)
            d_fn = e1.q.arrivalLik (arrv, e1)
            product = pwfun.add (a_fn, d_fn)
            return product
        else:
            raise NotImplementedError()
        
    NUM_TRIES = 1  # Number of tries for random GGk initializations
        
    def gibbs_initialize (self, a, samplePaths=False):
#        return self.gibbs_initialize_via_mmk (a)
        return initializer (self, a) 
#        return self.gibbs_initialize_via_dfn3 (a) 

    def gibbs_initialize_via_mmk (self, a):
        cdef Arrivals arrv2 = self.gibbs_initialize_via_rejection (a.duplicate())
        cdef Queue q
   
        old_tmpls = self.markovize (arrv2)
         
        arrvls = self.gibbs_resample (arrv2, burn=9, num_iter=1)

        self.demarkovize (old_tmpls)

        cdef Arrivals result = arrvls[-1]
        cdef Event evt
        for evt in result:
            evt.q.initialize_auxiliaries (evt)

        return result

    def markovize (self):
        cdef Queue q
        old_tmpls = []
        for q in self.queues:
            old_tmpls.append (q.service)
            q.service = queues.DistributionTemplate (distributions.Exponential ())
        theta = [ f.mean() for f in old_tmpls ]
        self.parameters = theta
        print "gibbs_initialize: mini_params %s" % self.parameters
        return old_tmpls
        
    def markovize_from_arrivals (self, a):
        cdef Queue q
        old_tmpls = []
        for q in self.queues:
            old_tmpls.append (q.service)
            q.service = queues.DistributionTemplate (distributions.Exponential ())
        self.estimate (a)
        print "gibbs_initialize: mini_params %s" % self.parameters
        return old_tmpls

    def demarkovize (self, old_tmpls):
        for i in xrange(len(old_tmpls)):
            (<Queue>self.queues[i]).service = old_tmpls[i]

    def gibbs_initialize_via_rejection (self, Arrivals arrv):
        # Initialize by sampling response times from initial distribution
        #  Assumes that data is sampled by task

        cdef Queue q
        cdef Event evt, evt_new
        
        max_backtracks = 10000

        result = None
        
        s_initial = self.fsm.initial_state()
        
        result_snapshot = Arrivals (self)
        heap_snapshot = self.heap_from_arrivals (arrv)
        eid_snapshot = 0
        s_max = numpy.inf
        
        for rep in xrange (max_backtracks):
            eid = eid_snapshot
            heap = heap_snapshot[:]
            result = result_snapshot.duplicate()
        
            while heap:
                a_i, tid, sid, evt = heappop (heap)
            
                if evt:                
                    # it was an observed event
                    d_min = evt.q.when_clear (result)
                    if d_min > evt.d:
                        # oops. backtrack
#                        print "BACKTRACK: Wanted to add %s to [s_max %f]\n%s" % (evt, s_max, result)
                        prev_d = evt.q.when_clear (result_snapshot)
                        s_max = evt.d - prev_d
                        result = None
                        break
                    else:                    
                        evt_new = evt.duplicate()
                        evt_new.eid = eid
                        result.append (evt_new)                
                        evt_new.q.recompute_service (evt_new)
                        
                        # create new snapshot
                        eid_snapshot = eid+1
                        heap_snapshot = heap[:]
                        result_snapshot = result.duplicate()
#                        print "SNAPSHOT: %s" % result_snapshot
                        s_max = numpy.inf
                else:    
                    # create new event
                    qid = self.fsm.sample_one_obs (sid)
                    q = self.queues[qid]
                    
                    if sid == s_initial: # sample initial arrivals
                        s_i = numpy.inf
                        while s_i >= s_max:
                            s_i, aux = q.sample_service()
                    else:
                        s_i = 0

                    evt = Event (eid, tid, qid, q, a_i, s_i, 0.0, sid, obs_a=0, obs_d=0, proc=0, aux=aux)
                    result.append (evt)
                    q.recompute_departure (result, evt)

                    s_next = self.fsm.sample_state (sid)
                    if s_next: heappush (heap, (evt.d, tid, s_next, None))

                eid = eid + 1

            if result: break  # we made it!
            
        # max backtracks reached, no initialization
        print "Gibbs initialization: Num backtracks = %d" % rep
        if not result:
            raise Exception ("Could not initialize\n%s" % arrv)
            
        return result


    # Event-wise slice sampling
    
    def slice_resample (self, Arrivals arrv, int num_iter, report_fn=None, return_arrv=False, reverse=False, return_lp=False):

        cdef int ne = 0
        cdef int gi
        cdef Event evt, evt_next
        cdef double upper, d_new
        cdef double lik

        global nevents 

        all_arrv = []

        lik = self.log_prob (arrv)
        Tmax = 2 * max([ e.d for e in arrv ])
        if DEBUG: print "slice_resample: initial log_prob ", lik
        if DEBUG: print "slice_resample: Tmax ", Tmax

        for gi in range(num_iter):
            
            # resample latent variables for all factors
            for evt in arrv:
                evt.q.resample_auxiliaries (evt)
                    
            allEvents = arrv.all_events()
            if is_slice_sorted: allEvents.sort(reverse=True, key=get_departure)
            if reverse: allEvents.reverse()

            for evt in allEvents:
                evt_next = evt.next_byt
                if not evt.obs_d and (not evt_next or not evt_next.obs_a):
                    if DEBUG: print "=========\nslice_resampling EID: %d\n%s" % (evt.eid, format_six(evt))

                    if lik == -INFINITY:
                        print "Before update: lik is infinity"
                        print "    EID: %d\n%s" % (evt.eid, format_six(evt))
                        arrv.validate()

                    ne += 1
                    dfn = GGkGibbs (self, arrv, evt, lik).dfn()
                    upper = self.determine_upper (evt)
                    upper = min(upper, Tmax)

                    d_new = _slice (evt.d, dfn, evt.a, upper)
                    lik = apply_departure (arrv, evt, d_new, lik)

                    if lik == -INFINITY:
                        print "After update: lik is infinity"
                        print "    EID: %d\n%s" % (evt.eid, format_six(evt))
                        arrv.validate()

                    nevents += 1

                    if DEBUG:
                        print "Sampled: ", d_new
                        arrv.validate()

            if report_fn:
                report_fn (self, arrv, gi, lik)
            if return_arrv:
                a2 = arrv.duplicate()
                all_arrv.append (a2)

        # sanity check
#        arrv.validate()
#        assert abs(lik - self.log_prob (arrv)) < 1e-5, "Log prob mismatch! Cached %.5f true %.5f" % (lik, self.log_prob(arrv))

        if not return_arrv: all_arrv = [arrv]

        if return_lp:
            return all_arrv, lik
        else:
            return all_arrv

    cdef double determine_upper (self, Event evt):
        next = evt.next_by_task()
        return next.d if next else numpy.inf

    cdef object dfn_range (self, Arrivals arrv, Event evt):
        return (evt.a, self.determine_upper(evt))
                            

    ## Computing the exact Gibbs kernel
    def gibbs_kernel (self, Arrivals arrv_from, Arrivals arrv_to, double lik=-INFINITY):

      cdef Event evt, evt_next
      cdef double upper, lower, d_new
      cdef double prob = 0, Z
       
      if lik == -INFINITY:
           lik = self.log_prob (arrv_from)

      arrv_new = arrv_from.duplicate()
      evtlst = arrv_new.all_events()
#        arrivals.sort_by_departure (evtlst)

        # logic to really try to find a good order

      deferred = arrv_new.all_events()
      while deferred:
        evtlst = deferred
        deferred = []
        ne = len (evtlst)
        for evt in evtlst:
            evt_next = evt.next_byt
            if not evt.obs_d and (not evt_next or not evt_next.obs_a):
                gibbs = GGkGibbs (self, arrv_new, evt, lik)
                log_fn = gibbs.dfn()

                lower = evt.a
                upper = self.determine_upper (evt)
                evt_to = arrv_to.event (evt.eid)
                d_new = evt_to.d
                lf_new = log_fn(d_new)

                # Not reachable
                if d_new <= lower or upper <= d_new:
                    print "gibbs_kernel: Kernel 0: %.5f not in %.5f .. %.5f\n  %s" % (d_new, lower, upper, evt)               
                    deferred.append (evt)
                    continue
#                    return -numpy.inf
                if lf_new == -numpy.inf: 
                    print "gibbs_kernel -inf: Kernel 0\nFROM  %s\nTO    %s" % (evt, evt_to)
                    deferred.append (evt)
                    continue
#                    return -numpy.inf

                Z = gibbs.integrate(lf_new) # pre-normalize by lf_new

                if Z == 0:
                    print "Bad event:", evt
                    print "Initial log prob", self.log_prob (arrv_from)
                    print "Log prob (from scratch)", self.log_prob (arrv_new)
                    if upper == numpy.inf: upper = lower + 5
                    x = lower
                    for i in range(1000):
                        print x, log_fn(x)
                        x += (upper - lower) / 1000
                    print "arrv_new"
                    arrv_new.write_csv(sys.stdout)
                    raise Exception ("Huh? Z=0. Lik ", lik)

#                print "d: %.5f  lp: %.5f  log(Z): %.5f  lower: %.5f  upper: %.5f\n   %s\n   %s" % (d_new, lf_new, log(Z), lower, upper, evt, evt_to)

                prob +=  - log(Z) # lf_new cancels out prenormalization

                lik = apply_departure (arrv_new, evt, d_new, lik)
        if len(deferred) == ne:  # no progress
            print "gibbs_kernel=0 : Could not process events"
            return -numpy.inf
        
      return prob

    def parameter_kernel (self, theta_from, theta_to, arrvl):
        """Returns the value of the transition kernel used by sample_parameters
           when it is given arrv."""
        if isinstance (arrvl, Arrivals):
            arrvl = [ arrvl ]

        cdef Queue q
        cdef Arrivals a
        cdef double ret = 0
        cdef FactorTemplate factor

        for q in self.queues:
            evtl = []
            for a in arrvl:
                evtl.extend (a.events_of_queue (q))

            n = q.service.numParameters()
            mini_from = theta_from[0:n]
            mini_to = theta_to[0:n]

            factor = q.service
            ret += factor.parameter_kernel (mini_from, mini_to, evtl)

            theta_from = theta_from[n:]
            theta_to = theta_to[n:]

        return ret

    ########################
    ## Parameters support ##
    ########################
    
    # todo: assumes exactly one factor per queue.  Not sure what the best way to generalize this is
    def estimate (self, arrvl):
        """Computes ML estimates for the parameters of all queues from the given arrivals list
           Returns the paramters, and also updates the parameters of this net."""
        if isinstance (arrvl, Arrivals):
            arrvl = [ arrvl ]

        cdef Queue q
        cdef Arrivals a

        for q in self.queues:
            evtl = []
            for a in arrvl:
                evtl.extend (a.events_of_queue (q))
            (<FactorTemplate>q.service).estimate (evtl)
        return self.parameters
        
    def estimate_from_obs (self, Arrivals arrv):
        """Computes ML estimates for the parameters of all queues from the given arrivals list.
           Returns the paramters, and also updates the parameters of this net."""
        cdef FactorTemplate f
        cdef Queue q
        for q in self.queues:
            f = q.service
            evts = [ e for e in arrv.events_of_queue (q) if (e.obs_a and e.obs_d) ]
            f.estimate (evts)
        return self.parameters

    def estimate_from_events (self, evts):
        """Computes ML estimates for the parameters of all queues from the given list of events.
           Returns the paramters, and also updates the parameters of this net."""
        cdef FactorTemplate f
        cdef Queue q
        for q in self.queues:
            f = q.service
            evts_of_q = [ e for e in evts if e.queue() == q ]
            f.estimate (evts_of_q)
        return self.parameters
        
    def sample_parameters (self, arrvl):
        """Samples the parameters from the posterior given the full obseravtions in ARRVL.
           Returns the paramters, and also updates the parameters of this net."""
        if isinstance (arrvl, Arrivals):
            arrvl = [ arrvl ]

        cdef Queue q
        cdef Arrivals a

        for q in self.queues:
            evtl = []
            for a in arrvl:
                evtl.extend (a.events_of_queue (q))
            (<FactorTemplate>q.service).sample_parameters (evtl)
        return self.parameters

    property parameters:
    
        def __get__ (self):
            r = []
            cdef Queue q
            cdef FactorTemplate t
            for q in self.queues:  
                t = q.service
                r.extend (t.parameters)
            return r
            
        def __set__ (self, v):
            r = v[:]
            cdef Queue q
            for q in self.queues:
                n = q.service.numParameters()
                q.service.setParameters (r[0:n])
                r = r[n:]
            


    # order proposal 
    
    # N.B. The acceptance ratio for this is certainly wrong
    def remove_task_proposal (self, Arrivals arrv_old):
        cdef Arrivals arrv = arrv_old.duplicate()
        cdef Event last_evt = arrv.final_arrival_by_qid (0)
        cdef Event evt
        cdef Queue q

        if self.eidx < 0:
            self.eidx = arrv_old.max_eidx() + 1

        if DEBUG_RTP:
            print "RTP"
            print arrv
            arrv.validate()

        tasks = arrv.all_tasks()
        k = last_evt.tid
        while k == last_evt.tid or arrv.events_of_task(k)[0].obs_d:
            k = sampling.random_elt (tasks)
        task_old = arrv.events_of_task(k)

        if DEBUG_RTP: print "Deleting task %s\n  %s" % (k, task_old)

        # Step 1: delete old task
        for evt in task_old:
            q = evt.q
            d_new = evt.wait() + evt.a
            dl = q.diffListForDeparture (evt, d_new)
            arrv.applyDiffList (dl)
        arrv.delete_task (k)

        if DEBUG_RTP: arrv.validate()

        # Step 2: compute random path for task
        path = self.sample_path ()
            
        cdef Event prev_evt
        cdef double log_trans1 = 0, log_trans0 = 0

        # Step 3: add events back, one_by_one, using departure proposal
        evt = task_old[0]
        global nrtp_reject, nrtp_calls

        while True:
            q = evt.q
            arrv.insert (evt)
            q.recompute_departure (arrv, evt)
            
            if DEBUG_RTP:
                print "++++++++++++++++\nRTP: Resampling departure"
                print evt.prev_byq
                print evt
                print evt.next_byq

            dfn = q.departureLik (arrv, evt)
            d_new, diff_list = sample_d_until_acceptance (self, arrv, evt)

            if DEBUG_RTP:
                print "+++++++++++++++++"
                print "RTP: Adding event"
                print diff_list

            diff_list = q.diffListForDeparture (evt, d_new)
            arrv.applyDiffList (diff_list)
            if evt.s == 0.0:
                if DEBUG_RTP: print "Rejecting: zero service on %s" % evt
                nrtp_calls += 1
                nrtp_reject += 1 
                return arrv_old
            # todo: log_trans0, log_trans1

            if  len(path) == 0: break

            # compute new evt
            sid, qid = path.pop(0)
            evt = Event (self.eidx, k, qid, self.queues[qid], d_new, 0.0, d_new, sid, 0)
            self.eidx += 1

        if DEBUG_RTP: 
            print "New task ", arrv.events_of_task(k)

        # step 4: acceptance ratio
        pi0 = self.log_prob (arrv_old) + log_trans1
        pi1 = self.log_prob (arrv) + log_trans0

        if DEBUG_RTP: print "pi0: %.10f  pi1: %.10f" % (pi0, pi1)

        nrtp_calls += 1

        u = rand()
    
        if u < math.exp(pi1 - pi0):
            # accept
            if DEBUG_RTP: arrv.validate()
            return arrv
        else:
            # reject
            nrtp_reject += 1
            return arrv_old
    
    # samples a path from the FSM.  Does not include initial state in
    # the returned path
    def sample_path (self):
        path = []
        sid = self.fsm.sample_state (self.fsm.initial_state())
        while sid:
            qid = self.fsm.sample_one_obs (sid)
            path.append ( (sid,qid) )
            sid = self.fsm.sample_state (sid)
        return path

    cpdef double log_prob (self, Arrivals arrv) except *:
        cdef double result = 0
        cdef Queue q
        for q in self.queues:
            result += q.allEventLik (arrv)
        return result
                    

# end (Qnet)
###################################################################


###################################################################
# Utilities

def is_sampled_by_task (arrv):
    tasks_in = {}
    tasks_out = {}

    for evt in arrv:            
        if evt.obs_d:
            if not evt.obs_a: return False
            if evt.tid in tasks_out: return False
            tasks_in[evt.tid] = True
        else:
            if evt.obs_a: return False
            if evt.tid in tasks_in: return False
            tasks_out[evt.tid] = True
            
    return True

def has_tau (Event evt):
    cdef Event byq = evt.prev_byq
#    if evt.prev_byt != None: return False
    if byq:
        return not evt.obs_a or not byq.obs_a
    else:
        return not evt.obs_a
        
def create_qp (a, reg):
    cdef Arrivals arrv = a
    ne = arrv.max_eidx()
    nt = arrv.num_tasks()

    qp = stupidlp.QP()
        
    cdef Event evt
    for evt in arrv:
        eid = evt.eid
        if has_tau(evt):
            qp.add_var_ge ("TAU%d"%eid, 1e-9)

    # objective
    for evt in arrv:
        byq = evt.prev_byq
        if has_tau (evt):
            qp.add_objsq (1.0, "TAU%d" % evt.eid)            
        
    for evt in arrv:
        evt.q.add_initialization_constraints (qp, evt)

    return qp
            
cdef setup_initialization_lp (Arrivals arrv, reg):
        lp = stupidlp.LP()
        r_means = compute_obs_means (arrv)   
        tau_mean = compute_obs_mean_ia (arrv)
        if tau_mean < 1e-5: tau_mean = 1e-5
        
        ne = arrv.max_eidx()
        nt = arrv.num_tasks()

        lp.add_vars ("set E := 0 .. %d;\n" % ne)

        lp.add_vars ("var d{i in E} >= 0;")
        lp.add_vars ("var r{i in E} >= 0;")
        lp.add_vars ("var tau{i in E} >= 0;")
        lp.add_vars ("var a{i in E} >= 0;")
        lp.add_vars ("var slack1{i in E} >= 0;")
        lp.add_vars ("var slack2{i in E} >= 0;")
        lp.add_vars ("var slack3{i in E} >= 0;")
        lp.add_vars ("var slack4{i in E} >= 0;")

        lp.set_objective ("minimize cost: sum{i in E} slack1[i] + sum{i in E} slack2[i] + sum{i in E} slack3[i] + sum{i in E} slack4[i];\n")

        # Constraint 1. d_e = a_e + r_e
        lp.add_constraint ("s.t. allr{i in E}: d[i] = a[i] + r[i];")
        
        cdef Event evt, byq, byt
        
        # Constraint 1b. d_e > a_e
        for evt in arrv:
            if (not evt.obs_d) or (not evt.obs_a):
                lp.add_constraint ("s.t. d_ge_a%d: d[%d] - a[%d] >= %s;\n" % (evt.eid, evt.eid, evt.eid, reg))
            
        # Constraint 2. d_pi(e) = a_e
        for evt in arrv:
            byt = evt.prev_byt
            if not byt:
                lp.add_constraint ("s.t. a_pi%d: a[%d] = 0;\n" % (evt.eid, evt.eid))
            else:
                lp.add_constraint ("s.t. a_pi%d: a[%d] = d[%d];\n" % (evt.eid, evt.eid, byt.eid))
                
        # Constraint 3. d_e >= d_(rho(e))
        for evt in arrv:
            byq = evt.q.previous_departing (evt)
            if byq:
                if (not evt.obs_d) or (not byq.obs_d):
                    lp.add_constraint ("s.t. d_rho%d: d[%d] - d[%d] >= %s;\n" % (evt.eid, evt.eid, byq.eid, reg))

        # Constraint 3b. a_e >= a_(rho(e))
        for evt in arrv:
            byq = evt.prev_byq #  evt.q.previous_departing (evt)
            if byq and evt.prev_byt:  # don't do it for initial events
                if (not evt.obs_a) or (not byq.obs_a):
                    lp.add_constraint ("s.t. a_rho%d: a[%d] - a[%d] >= %s;\n" % (evt.eid, evt.eid, byq.eid, reg))
                 
        # Constraint 4. Define tau
        for evt in arrv:
            if not evt.prev_byt:
                byq = evt.q.previous_departing (evt)
                if byq:
                    lp.add_constraint ("s.t. tau%d: d[%d] + tau[%d] = d[%d];\n" % (evt.eid, byq.eid, evt.eid, evt.eid))
                else:
                    lp.add_constraint ("s.t. tau%d: d[%d] = tau[%d];\n" % (evt.eid, evt.eid, evt.eid))

        # Constraint 5. slacks on R
        for evt in arrv:
            val = r_means [evt.qid]
            lp.add_constraint ("s.t. slack1%d: r[%d] + slack1[%d] - slack2[%d] = %s;\n" % (evt.eid, evt.eid, evt.eid, evt.eid, val))

        # Constraint 5. slacks on tau
        for evt in arrv:
            if not evt.prev_byt:
                lp.add_constraint ("s.t. slack_b%d: tau[%d] + slack3[%d] - slack4[%d] = %s;\n" % (evt.eid, evt.eid, evt.eid, evt.eid, tau_mean))
                
        # Constraint 4. Equality if event observed
        for evt in arrv:
            if evt.obs_d:
                lp.add_constraint ("s.t. dobs_%d: d[%d] = %.17g;\n" % (evt.eid, evt.eid, evt.d))
            if evt.obs_a:
                lp.add_constraint ("s.t. aobs_%d: a[%d] = %.17g;\n" % (evt.eid, evt.eid, evt.a))

        return lp
    
        
def compute_obs_means (arrv):
    cdef Event evt
    tots = [1e-5] * arrv.num_queues()
    Ns = [1] * arrv.num_queues()
    for evt in arrv:
        if evt.obs_d and evt.obs_a:
            tots[evt.qid] += (evt.d - evt.a)
            Ns[evt.qid] += 1
    for i in xrange(len(tots)):
        tots[i] = tots[i] / Ns[i]
    return tots

def compute_obs_mean_ia (arrv):
    cdef Event evt
    last_a = 0.
    N = 0
    for evt in arrv:
        if not evt.prev_byt:
            if evt.obs_a:
                last_a = max (last_a, evt.d)
                N = N + 1
    if N > 0:
        return last_a / N
    else: return 0
    
def compute_mean_ia (arrv):
    cdef Event evt
    times = []
    for evt in arrv.initial_events():
        if evt.next_byq:
            time_next = evt.next_byq.d - evt.d
            if time_next >= 0:   # zeroing out unobserved event can do this
                times.append (time_next)
    return numpy.mean (times)
    
def sorted_uniforms (n, L, U):
    l = [ uniform(L,U) for i in xrange (n) ]
    l.sort()
    return l

def setup_minilp (Arrivals arrv, reg):
    cdef Event evt

    lp = stupidlp.LP()
        
    ne = arrv.max_eidx()
    nt = arrv.num_tasks()

    lp.add_vars ("set E := 0 .. %d;\n" % ne)

    lp.add_vars ("var d{i in E} >= 0;")
    lp.add_vars ("var a{i in E} >= 0;")

    lp.set_objective ("minimize cost: 0;\n")
        
    # Constraint 1. d_e > a_e
    for evt in arrv:
        if (not evt.obs_a) or (not evt.obs_d):
            lp.add_constraint ("s.t. d_ge_a%d: %s - %s >= %s;\n" % (evt.eid, dstr(evt), astr(evt), reg))
            print "s.t. d_ge_a%d: %s - %s >= %s;" % (evt.eid, dstr(evt), astr(evt), reg)
            
    # Constraint 2. d_pi(e) = a_e
    for evt in arrv:
        byt = evt.prev_byt
        if byt:
            if (not evt.prev_byt.obs_d) or (not evt.obs_a):
                lp.add_constraint ("s.t. a_pi%d: %s = %s;\n" % (evt.eid, dstr(byt), astr(evt)))

    # Constraint 3. Initial Queue
    for evt in arrv:
        byt = evt.prev_byt
        if byt is None and not evt.obs_a:
            lp.add_constraint ("s.t. a_ini%d: a[%d] = 0;\n" % (evt.eid, evt.eid))

    # Constraint 4: ORDER
    for evt in arrv:
        byq = evt.prev_byq
        if byq is not None:
            order = evt.q.queue_order()
            if order == arrivals.ARRV_BY_A:
                if (not evt.obs_a) or (not byq.obs_a):
                    lp.add_constraint ("s.t. order%d: %s >= %s;\n" % (evt.eid, astr(evt), astr(byq)))
            elif order == arrivals.ARRV_BY_D:
                if (not evt.obs_d) or (not byq.obs_d):
                    lp.add_constraint ("s.t. order%d: %s >= %s;\n" % (evt.eid, dstr(evt), dstr(byq)))
            else:
                raise Exception ("Huh? Unknown queue order %d" % order)

    return lp

def astr (evt):
    return "a[%d]" % evt.eid if not evt.obs_a else "%.10f" % evt.a
def dstr (evt):
    return "d[%d]" % evt.eid if not evt.obs_d else "%.10f" % evt.d
def compactify_soln (arrv, soln):
    for evt in arrv:
        if evt.obs_a: soln["a[%d]" % evt.eid] = evt.a
        if evt.obs_d: soln["d[%d]" % evt.eid] = evt.d

cdef format_six (Event evt):
    cdef Event byt, byq, byt_minus1, byt_plus1, next_byq
    byq = evt.prev_byq
    byt = evt.next_byt
    byt_minus1 = byt.prev_byq if byt else None
    byt_plus1 = byt.next_byq if byt else None
    next_byq = evt.next_byq
    return "%-75s%s\n%-75s%s\n%-75s%s" % (byq, byt_minus1, evt, byt, next_byq, byt_plus1)
    
def _compute_task_order (Event evt):
    i = 0
    while evt.prev_byt:
        evt = evt.prev_byt
        i = i + 1
    return i
    
cdef double call_sampler (Pwfun fn, double x_initial) except *:
    if sampler == ARS_PWFUN:
        return sampling.arms_pwfun (fn, x_initial)
    elif sampler == SLICE_PWFUN:
        return _slice_pwfun (fn, x_initial)
    elif sampler == SLICE:
        return sampling.slice (fn, x_initial, lower=fn.L(), upper=fn.U(), thin=2)
    elif sampler == STUPID_REJECT:
        return sampling.rejection (fn)
        


def display_diff_list (arrv, Event e, diff_list):
    result = "EID: %d\n" % e.eid
    cdef Event evt
    for evt in diff_list:
        evt_old = arrv.event (evt.eid)
        result += "OLD %s\nNEW %s\n\n" % (evt_old, evt)

    return result

def any_zero_service (diff_list):
    for evt in diff_list:
        if evt.s < 0: 
            return True
    return False



# closure emulation

NUM_INTEGRATION_POINTS = 3

# true gibbs distribution
cdef class GGkGibbs:
    cdef double lik0
    cdef Event evt, evt_next
    cdef Arrivals arrv
    cdef Qnet net

    def __init__(self, Qnet net, Arrivals arrv, Event evt, double lik0):
        self.lik0 = lik0
        self.net = net
        self.evt = evt
        self.evt_next = evt.next_byt
        self.arrv = arrv

    def inner_dfun (self, double d):
        cdef double lik = self.lik0

        cdef Event e_new, e_old
        cdef Queue q0 = self.evt.q
        cdef Queue q1 = self.evt_next.q if self.evt_next else None
        dl0 = q0.diffListForDeparture (self.evt, d)
        dl1 = q1.diffListForArrival(self.evt_next, d) if q1 else []

        global diff_len, num_diff
        diff_len += len(dl0)
        diff_len += len(dl1)
        num_diff += 1

        lik += q0.likelihoodDelta (self.arrv, dl0)
        lik += q1.likelihoodDelta (self.arrv, dl1) if q1 else 0

        return lik

    def dfn (self): return self.inner_dfun

    def integrate (self, C=0):
        fn = netutils.expify(self.dfn(), -C)
        lower = self.evt.a
        upper = self.net.determine_upper (self.evt)
        
        # heuristic: add integration points for past and next 5
        points = [lower, upper]
        self.add_points (points, self.evt, lower, upper, get_departure)
        if self.evt_next: self.add_points (points, self.evt_next, lower, upper, get_arrival)

        points.sort()
        return misc.integrate_by_points (fn, points)

    def add_points (self, points, Event e0, double lower, double upper, fn):
        cdef Event e_curr = e0
        for i in range(NUM_INTEGRATION_POINTS):            
            if not e_curr: break
            val = fn(e_curr)
            if lower <= val and val <= upper:
                points.append (val)
            e_curr = e_curr.prev_byq
        e_curr = e0.next_byq
        for i in range(NUM_INTEGRATION_POINTS):            
            if not e_curr: break
            val = fn(e_curr)
            if lower <= val and val <= upper:
                points.append (val)
            e_curr = e_curr.next_byq


# Initializers

def gibbs_initialize_via_dfn1 (net, Arrivals arrv):
    result = gibbs_initialize_via_dfn (net, arrv, 1.0)
    cdef Event evt
    for evt in result:
        evt.q.initialize_auxiliaries (evt)
    result.regenerate_all_caches()
    return result

def gibbs_initialize_via_dfn2 (net, Arrivals arrv):
    inner = gibbs_initialize_via_dfn (net, arrv, 1.0)
    result = gibbs_initialize_via_dfn (net, inner, 1.0)
    cdef Event evt
    for evt in result:
        evt.q.initialize_auxiliaries (evt)
    result.regenerate_all_caches()
    return result

def gibbs_initialize_via_dfn3 (net, Arrivals arrv):
    result = gibbs_initialize_via_dfn (net, arrv, 1.1)
    cdef Event evt
    for evt in result:
        evt.q.initialize_auxiliaries (evt)
    result.regenerate_all_caches()
    return result

def insert_implied_tasks (Qnet net, Arrivals arrv):
    cdef Event evt, evt0, evt1

    eid = arrv.max_eidx()+1
    d0 = 0
    tid0 = -1
    evt2 = arrv.first_arrival_by_qid (INITIAL_QID)
    initial_q = net.queues[INITIAL_QID]

    eps = 1e-10

    while evt2:
        tid2 = evt2.tid
        if DEBUG_INITIAL: print "INIT  %s\n" % evt2
        if DEBUG_INITIAL: print "TIDS:",tid0, tid2

        for tid1 in range(tid0+1, tid2):            
            evt1 = Event (eid, tid1, INITIAL_QID, initial_q, 0.0, 0.0, d0 + eps, INITIAL_STATE, obs_a=0, obs_d=0)
            if DEBUG_INITIAL: print tid1, evt1
            eid += 1
            arrv.insert (evt1)
            evt1.q.recompute_service (evt1)
            # increment last event
            d0 = evt1.d

        evt0 = evt2
        d0 = evt0.d
        tid0 = evt0.tid
        evt2 = evt0.next_byq

    return arrv

def initialize_inserting_implied_tasks (Qnet net, Arrivals arrv):
    insert_implied_tasks (net, arrv)
    return initializer (net, arrv)

def gibbs_initialize_via_dfn (Qnet net, Arrivals arrv, scale):
    # Initialize by sampling response times from initial distribution
    #  Assumes that data is sampled by task

    cdef Queue q
    cdef Event evt, evt_new, evt_prev, evt1
    cdef Arrivals result = Arrivals(net)

    s_initial = net.fsm.initial_state()
    eid = arrv.max_eidx()+1

    dists = net.markovize_from_arrivals (arrv)

    if DEBUG_INITIAL: 
        print "Parameters: ", net.parameters
        print "Adding observed events"

    # first add observed events
    for tid in arrv.all_tasks():
        evts = arrv.events_of_task(tid)
        if evts[0].obs_d:
            for evt in evts:
                evt_new = evt.duplicate()
                evt_new.s = 0.0
                evt_new.d = evt_new.a + evt_new.wait()   
                result.insert (evt_new)
                diff_list = evt_new.q.diffListForDeparture (evt_new, evt.d)
                result.applyDiffList (diff_list)
                evt_new.q.initialize_auxiliaries (evt_new)

    result.validate()

    if DEBUG_INITIAL: 
        print "Adding unobserved events"
        print result

    # now add unobserved events
    tasks = arrv.all_tasks()
    for i in xrange(len(tasks)):
        tid = tasks[i]
        evts = arrv.events_of_task(tid)
        if not evts[0].obs_d:

            path = net.sample_path()
            evt = (<Event> evts[0]).duplicate()

            # special case: initial arrival (maintain initial arrival order)
            if i > 0:
                eid_prev = (<Event> arrv.events_of_task(tasks[i-1])[0]).eid
                evt_prev = result.event (eid_prev)
                evt.s = 1e-10
                evt.d = evt_prev.d + 1e-10  
                if DEBUG_INITIAL:
                    print "Pre add:\n  prev %s\n   evt %s" % (evt_prev.dump(), evt.dump())
                
            result.insert (evt)

            q = evt.q
            q.recompute_service (evt)
            if evt.next_byq: (<Event> evt.next_byq).q.recompute_service (evt.next_byq)

            while True:
                d_new, diff_list = sample_d_until_acceptance (net, result, evt, service_only=True, scale=scale)

                if DEBUG_INITIAL:
                    print "================="
                    print evt
                    print "D_NEW ", d_new
                    print "::::::::::::::::::::::::::"
                    print evt.prev_byq.dump() if evt.prev_byq else None
                    print evt.dump()
                    print evt.next_byq.dump() if evt.next_byq else None
                    print "::::::::::::::::::::::::::"
                    print "diff_list: ", diff_list
                    #print result
                    print ":::::::  Event %d done :::::::::::::::::::\n"  % evt.eid

                result.applyDiffList (diff_list)

                if DEBUG_INITIAL:  
                    print "diff_list ", diff_list
                    result.validate()

                if len(path) == 0: break

                # compute new evt
                sid, qid = path.pop (0)
                if DEBUG_INITIAL: print "TID: %s  SID: %s  QID: %s" % (tid, sid, qid)

                q = net.queues[qid]
                evt = Event (eid, tid, qid, q, d_new, 0.0, d_new, sid, 0)
                eid += 1

                result.insert (evt)
                q.recompute_departure (result, evt)
                q.update_caches_for (evt)
                q.initialize_auxiliaries (evt)

    net.demarkovize (dists)

    return result


def init_insert_observed_tasks (Arrivals arrv, Arrivals result, add_initial=True):
    cdef Event evt_new, evt
    for tid in arrv.all_tasks():
        evts = arrv.events_of_task(tid)
        if evts[0].obs_d:
            for evt in evts:
                if add_initial or evt.qid > 0:
                    evt_new = evt.duplicate()
                    evt_new.s = 0.0
                    evt_new.d = evt_new.a + evt_new.wait()   
                    if DEBUG_INITIAL: print "Adding observed event ", evt_new
                    result.insert (evt_new)
                    diff_list = evt_new.q.diffListForDeparture (evt_new, evt.d)
                    result.applyDiffList (diff_list)
                    evt_new.q.initialize_auxiliaries (evt_new)

    # now add unobserved initial events
def init_add_initial_events (Qnet net, Arrivals arrv, Arrivals result):
    cdef Event evt
    cdef Queue q
    cdef double d_prev
    heap = []
    tasks = arrv.all_tasks()
    for i in xrange(len(tasks)):
        tid = tasks[i]
        evts = arrv.events_of_task(tid)
        if not evts[0].obs_d:
            path = net.sample_path()
            evt = (<Event> evts[0]).duplicate()

            # special case: initial arrival (maintain initial arrival order)
            if i > 0:
                eid_prev = (<Event> arrv.events_of_task(tasks[i-1])[0]).eid
                evt_prev = result.event (eid_prev)
                d_prev = evt_prev.d
            else:
                d_prev = 0

            evt.s = 1e-10
            evt.d = d_prev + 1e-10  
            if DEBUG_INITIAL: print "Adding initial event ", evt

            result.insert (evt)
            q = evt.q
            q.recompute_service (evt)
            if evt.next_byq: (<Event> evt.next_byq).q.recompute_service (evt.next_byq)

            d_new, diff_list = sample_d_until_acceptance (net, result, evt, service_only=True, scale=1.0)
            result.applyDiffList (diff_list)

            heappush (heap, (d_new, tid, path))
    return heap

def gibbs_initialize_via_order (Qnet net, Arrivals arrv):
    # Initialize by sampling response times from initial distribution
    #  Assumes that data is sampled by task

    cdef Queue q
    cdef Event evt, evt_new, evt_prev, evt1
    cdef Arrivals result = Arrivals(net)

    s_initial = net.fsm.initial_state()
    eid = arrv.max_eidx()+1

    dists = net.markovize_from_arrivals (arrv)

    if DEBUG_INITIAL: 
        print "Parameters: ", net.parameters
        print "Adding observed events"

    # first add observed events
    init_insert_observed_tasks (arrv, result)

#   This actually need not hold for GG1R queues
#    result.validate()

    if DEBUG_INITIAL: 
        print "Adding unobserved events"
        print result

    # now add unobserved initial events
    heap = init_add_initial_events (net, arrv, result)

    while heap:
        a, tid, path = heappop(heap)
        sid, qid = path.pop (0)
        if DEBUG_INITIAL: print "TID: %s  SID: %s  QID: %s" % (tid, sid, qid)

        q = net.queues[qid]
        evt = Event (eid, tid, qid, q, a, 0.0, a, sid, 0)
        eid += 1

        result.insert (evt)
        q.recompute_departure (result, evt)
        q.update_caches_for (evt)
        q.initialize_auxiliaries (evt)

        d_new, diff_list = sample_d_until_acceptance (net, result, evt, service_only=True, scale=1.0)
        if diff_list is None:
            # abort thin task
            print "initialize_by_order: aborting task ", tid
            result.delete_task (tid)
            continue

        if DEBUG_INITIAL: print "DIFF LIST", "\n".join(map(str, diff_list))
        
#         for e0 in diff_list:
#             eold = result.event (e0.eid)
#             if eold.s > e0.s:
#                 print evt
#                 print "---> reducing service time\nOLD  %s\nNEW  %s" % (eold, e0)
#                 break

        result.applyDiffList (diff_list)

        if len(path) > 0:
            heappush (heap, (d_new, tid, path))
        
    net.demarkovize (dists)

    for evt in result:
        evt.q.initialize_auxiliaries (evt)
    result.regenerate_all_caches()

    result.recompute_all_service()

    print "Initialization done.  Service: ", qstats.mean_service (result)
    result.report ()

    return result

def gibbs_initialize_for_gg1r (Qnet net, Arrivals arrv):
    # Initialize by sampling response times from initial distribution
    #  Assumes that data is sampled by task

    cdef Queue q
    cdef Event evt, evt_new, evt_prev, evt1
    cdef Arrivals result = Arrivals(net)

    s_initial = net.fsm.initial_state()
    eid = arrv.max_eidx()+1

    dists = net.markovize_from_arrivals (arrv)

    if DEBUG_INITIAL: 
        print "Parameters: ", net.parameters
        print "Adding observed events"

    # first add observed events
    init_insert_observed_tasks (arrv, result)

    result.validate()

    # now add unobserved initial events
    heap = init_add_initial_events (net, arrv, result)

    # now add all other events with zero service
    unobs = []
    for elem in heap:
        a, tid, path = elem
        for sid, qid in path:
            sid, qid = path.pop (0)
            q = net.queues[qid]
            s = 0.01
            d = a + s
            evt = Event (eid, tid, qid, q, a, s, d, sid, 0)
            result.insert (evt)
            unobs.append (evt)
            eid += 1
            assert evt.next_byt is not evt, "Huh? %s" % evt
            a = d

    # sample service times
    for evt in unobs: 
        s_new, aux = evt.q.sample_service ()
        evt.s = s_new
        
    # iterate: 
    N=0
    to_fix = set ([ e.eid for e in unobs])
    while to_fix:
        eid = to_fix.pop()    
        evt = result.event (eid)
        d_old = evt.d
        evt.q.recompute_departure (result, evt)
        if evt.d != d_old:
            # add old next_byq
            if evt.next_byq:
                to_fix.add (evt.next_byq.eid)
            # resort in queue
            result.resort_byq (evt)
            if evt.next_byt:
                evt.next_byt.a = evt.d
                to_fix.add (evt.next_byt.eid)
            # add new next_byq
            if evt.next_byq:
                to_fix.add (evt.next_byq.eid)

        N += 1
        if N % 1000 == 0:
            print "gg1r_initialize ", len(to_fix)

    for evt in result:
        evt.q.initialize_auxiliaries (evt)
    result.recompute_all_service()
    result.regenerate_all_caches()
    
    net.demarkovize (dists)
    
    print "GG1R initialization done.  Mean service: ", qstats.mean_service (result)

    return result

def diff_list_definitely_bad (net, arrv, q, diff_list):
    for evt1 in diff_list:
        if evt1.s < 0: return True
    if numpy.isinf (q.likelihoodDelta(arrv, diff_list)):
        return True
    return False

def sample_d_until_acceptance (net, Arrivals arrv, Event evt, service_only=False, scale=1.0):
    cdef Queue q = evt.q

    dfn = q.departureLik (arrv, evt)

    d_initial = evt.d

    if evt.s == 0:
        delta = 1.0
        while delta > 0:
            diff_list = q.diffListForDeparture (evt, evt.d + delta)
            if not diff_list_definitely_bad (net, arrv, q, diff_list): break
            delta = delta/2.0
        d_initial += delta

    found_a_good_one = False
    cdef double s_true, wait = evt.wait()
    for j in xrange(100):
        if service_only:
            s_new, aux = q.sample_service()
            s_new = s_new / scale
            d_new = q.departure_of_service (evt, s_new)
        else:
            d_new = call_sampler (dfn, d_initial)

        diff_list = q.diffListForDeparture (evt, d_new)
        if not diff_list_definitely_bad (net, arrv, q, diff_list):
            found_a_good_one = True
            break

        # kind of a hack.  try to fix service being underinitialized, exp in GG1RSS
        s_true = -1
        for e0 in diff_list:
            if e0.eid == evt.eid:
                s_true == e0.s
                wait = e0.wait()
        if abs (s_new - s_true) < 0.1:
            break

    cdef Event evt_to_test
    if not found_a_good_one:
        # one more try
        evt_to_test = evt.next_byq
        while evt_to_test:
            d_new = evt_to_test.d + 1e-5
            diff_list = q.diffListForDeparture (evt, d_new)
            if not diff_list_definitely_bad (net, arrv, q, diff_list):
                found_a_good_one = True
                break
            evt_to_test = evt_to_test.next_byq

    if not found_a_good_one:        
        print "s_d_u_a: Error: Could not find good initial time for EVT %d\n%s" % (evt.eid, format_six(evt))
        d_new = d_initial
        diff_list = q.diffListForDeparture (evt, d_new)
        if diff_list_definitely_bad (net, arrv, q, diff_list):
            d_new = 0
            diff_list = None 

    return d_new, diff_list

def heap_from_arrivals (net, Arrivals arrv):
    heap = []
    cdef Event evt
    for evt in arrv:
        if evt.obs_d:
            heappush (heap, (evt.a, evt.tid, evt.state, evt))
        # add initial event if necessary
        if not evt.prev_byt and not evt.obs_d:
            s_initial = net.fsm.initial_state()
            heappush (heap, (0.0, evt.tid, s_initial, None))
    return heap

def gibbs_initialize_via_lp (net, Arrivals arrv_in):
    cdef Arrivals arrv
    cdef Event evt

    soln = dict()

    for reg in [ 1e-6, 1e-10, 0 ]:
        try:

            arrv = arrv_in.duplicate()

            lp = setup_initialization_lp (arrv, reg)
            soln = lp.solve()
            print "***********"

            if DEBUG_INITIAL:
                for k,v in soln.iteritems():
                    print "%s %.17g" % (k,v)

            for evt in arrv:
                evt.update_from_solution (soln)

            for evt in arrv:
                evt.q.recompute_service (evt)

            arrv.validate()        
            lp.cleanup()

            return arrv

        except Exception, e:
            if reg == 0:  # give up
                for k,v in soln.iteritems():
                    print "%s %s" % (k,v)
                print arrv
                raise e
            else:
                print "Initialization at %s failed: %s" % (reg, e)

def visualize_initialization_lp (net, Arrivals arrv):
    cdef Event evt, byq, byt

    f = open ("model.dot", "w")

    f.write ("digraph G {\n")

    for q in net.queues:

        f.write ("subgraph cluster_%s_D {\nlabel = \"%s\";\ncolor=blue;\n" % (q.name, q.name))

        for evt in arrv.events_of_queue (q):                
            eid = evt.eid

            # Constraint 3. d_e >= d_(rho(e))
            byq = evt.q.previous_departing (evt)
            if byq:
                if (not evt.obs_d) or (not byq.obs_d):
                    f.write ("d%d -> d%d;\n" % (byq.eid, eid))

        f.write("}\n")
        f.write ("subgraph cluster_%s_A {\nlabel = \"%s\";\ncolor=blue;\n" % (q.name, q.name))

        for evt in arrv.events_of_queue (q):
            eid = evt.eid

            # Constraint 3b. a_e >= a_(rho(e))
            byq = evt.q.previous_departing (evt)
            if byq and evt.prev_byt:  # don't do it for initial events
                if (not evt.obs_a) or (not byq.obs_a):
                    f.write ("a%d -> a%d;\n" % (byq.eid, eid))

        f.write("}\n")

    for evt in arrv:
        byt = evt.prev_byt
        if byt:
            f.write ("a%d -> d%d [color=gray60];\n" % (evt.eid, evt.eid))
            f.write ("d%d -> a%d [color=lawngreen];\n" % (byt.eid, evt.eid))
        if evt.obs_a:
            f.write ("a%d [style=filled, color=grey80];\n" % evt.eid)
        if evt.obs_d:
            f.write ("d%d [style=filled, color=grey80];\n" % evt.eid)

    f.write ("}\n")
    f.close ()

def gibbs_initialize_via_minilp (net, Arrivals arrv_in):
    cdef Arrivals arrv
    cdef Event evt

    arrv_in.validate()
    soln = dict()

    for reg in [ 0.1, 0.01, 1e-4, 1e-6, 1e-10, 0 ]:
        try:

            arrv = arrv_in.duplicate()
            arrv.recompute_all_service() # hack for cases in which multiple events have same departure time at same queeu
            arrv.validate()

            prev_dept_old = [ e.queue().previous_departing(e) for e in arrv ] 
            prev_arrv_old = [ e.queue().previous_arriving(e) for e in arrv ] 

            print "********** RUNNING MINI-LP ", reg
            lp = setup_minilp (arrv, reg)
            soln = lp.solve()
            compactify_soln (arrv, soln)
            print "***********"
            
            if DEBUG_INITIAL:
                for k,v in soln.iteritems():
                    print "%s %.17g" % (k,v)

            for evt in arrv:
                evt.update_from_solution (soln)

            for evt in arrv:
                evt.q.recompute_service (evt)

            arrv.regenerate_all_caches()
            arrv.validate()        
            lp.cleanup()

            old_tmpls = net.markovize_from_arrivals (arrv)      
            estimate_minilp_params (net, arrv)
            arrvls = net.slice_resample (arrv, 10)
            arrvls[-1].validate()
            net.demarkovize (old_tmpls)

            return arrvls[-1]

        except Exception, e:
            if reg == 0:  # give up
                for k,v in soln.iteritems():
                    print "%s %s" % (k,v)
                arrv.write_csv (sys.stdout)
                traceback.print_exc()
                raise e
            else:
                print "Initialization at %s failed: %s" % (reg, e)
                traceback.print_exc()


# requires Markovization
def estimate_minilp_params (net, arrv):
    mu_s = qstats.mean_task_service (arrv)
    N = qstats.mean_task_length (arrv)
    net.parameters = [mu_s / N] * len(net.parameters)
    print "minilp mini params", net.parameters


def gibbs_initialize_for_ps (net, arrv):
    cdef Event e, evt, dup
    cdef Queue q

    cdef Arrivals result = Arrivals(net)
    cdef int eid = arrv.max_eidx() + 1
    cdef double d_last = 0

    ievts = arrv.initial_events()
    ievts.sort()
    for evt in arrv.initial_events():
        if evt.obs_d:
            s = evt.d - d_last
            dup = Event (evt.eid, evt.tid, evt.qid, evt.q, 0.0, s, evt.d, evt.state, obs_a=1, obs_d=1)
            d_last = dup.d
            result.append (dup)
        else:
            dup = Event (eid, evt.tid, evt.qid, evt.q, 0, 0, d_last, evt.state, 0)
            result.append (dup)
            eid += 1
        (<Queue>dup.q).recompute_service(dup)

    result.regenerate_all_caches()
#    print "STEP0 ", result
#    result.validate()
    
    # spread out initial events
    ievts = result.events_of_queue (net.queues[0])
    next_obs = [0] * len(ievts)
    n_hidden = [0] * len(ievts)
    i = len(ievts)-1
    Z = numpy.inf
    nh = 0
    d_max = 0
    while i >= 0:
        evt = ievts[i]
        if not evt.obs_d:
            nh += 1
        else:
            next_obs[i] = Z
            n_hidden[i] = nh
            Z = evt.d
            d_max = max (d_max, evt.d)
            nh = 0
        i -= 1
    # estimate of mean service
    rate = len(ievts) / d_max
    # now set them up
    d_last = ievts[0].a
    delta = Z / (1+nh)
    for i in range(len(ievts)):
        evt = ievts[i]
        if not evt.obs_d:
            d_last += delta
            apply_departure (result, evt, d_last, 0)
        else:
            if numpy.isinf(next_obs[i]):
                delta = 1/rate
            else:
                delta = (next_obs[i] - evt.d) / (1+n_hidden[i])
            d_last = evt.d
    # finally reset service
    q0 = net.queues[0]
    for e in ievts:
        q0.recompute_service (e)
#    print "STEP1 ", result
#    result.validate()

    # sample extensions of tasks
    heap = []
    for evt in ievts:        
        if not evt.obs_d:
            s_next = net.fsm.sample_state (evt.state)
            if s_next: heappush (heap, (evt.d, evt.tid, s_next, SEC_ARRIVAL, None))
    net._sample_discrete_events (heap, result, False)

#    print "STEP2 ", result
#    result.validate()
    
    # add observed last so above processing doesn't screw with their departure times
    init_insert_observed_tasks (arrv, result, add_initial=False)
    result.recompute_all_service ()

#    print "STEP3 ", result
#    result.validate()
    
    return result

## scrap
    
# NUM_SWEETEN_ITER = 5
# def gibbs_initialize_via_spread (net, arrv):
#     old_tmpls = net.markovize ()
    
#     # set all service times to 0
#     cdef Event e
#     for e in arrv:
#         if not e.obs_d:
#             e.s = 0
#             (<Queue>e.q).recompute_departure (arrv, e)

#     # some of the observed events may need their service times udpated
#     arrv.recompute_all_service()

#     print arrv #GGG

#     # spread the intial events out evenly
#     q0 = net.queues[0]
#     ievts = arrv.events_of_queue (q0)
#     delta = -1
#     d_current = 0
#     obs_i = 0

#     for ei in range(len(ievts)):
#         e = ievts[ei]
#         if e.obs_d:
#             delta = -1
#             obs_i = ei
#             d_current = e.d
#         else:
#             if delta < 0:
#                 ei2 = find_next_obs (ievts, ei+1)
#                 if ei2 is None:
#                     delta = ievts[obs_i].d / obs_i
#                 else:
#                     delta = (ievts[ei2].d - ievts[obs_i].d) / (ei2 - obs_i)
#                     print "Delta ", delta, " ei2", ei2, "obs_i", obs_i, "ei", ei
                    
#             d_current += delta
#             apply_departure (arrv, e, d_current, 0)
                   
#     # sweeten using the Markovized service distributions
#     arrvl = net.slice_resample (arrv, NUM_SWEETEN_ITER, return_arrv=False, return_lp=False)

#     net.demarkovize (old_tmpls)

#     return arrvl[-1]

# def find_next_obs (evtl, idx):
#     while idx < len(evtl):
#         if evtl[idx].obs_d: return idx
#         idx += 1
#     return None

# def diff (l):
#     return [ (i1-i0) for i0,i1 in zip(l, l[1:]) ]

# def cumsum(l):
#     result = []
#     cum = 0
#     for i in l:
#         cum += i
#         result.append (cum)
#     return result

def null_initializer (net, arrv): return arrv

#cas?? initializer = gibbs_initialize_via_order
initializer = gibbs_initialize_via_lp


def get_departure(evt): return evt.d
def get_arrival(evt): return evt.a


cpdef double apply_departure (Arrivals arrv, Event evt, double d, double lik) except *:
    cdef Event evt_next

    dl0 = evt.q.diffListForDeparture (evt, d)

    evt_next = evt.next_by_task()
    if evt_next:
        dl0.extend (evt_next.q.diffListForArrival(evt_next, d))

    if DEBUG:
        print "diff_list was\n", display_diff_list (arrv, evt, dl0)

    return arrv.applyDiffList (dl0, lik)

