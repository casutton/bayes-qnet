# Functions for adding and removing queues and states to models,
#  keeping an associated Arrivals object up to date

import qnet
import arrivals
import qnetu
from numpy import linalg
import numpy
import hmm
import estimation
import sys
import mytime

# TODO: Arrivals conversion (should be easy)
def split_state (net, sname_old, sname_new1, sname_new2, qnames0, qnames1):
    """Splits a state in a Qnet.
         NET is a qnet object, and ARRV an associated arrivals object.
         SNAME_OLD is the state to split, into SNAME_NEW1 and SNAME_NEW2.
         (SNAME_OLD, SNAME_NEW1, and SNAME_NEW2 are string state names).
         Q0L and Q1L are a partition of the queues that SNAME_OLD
         can generate.  These are lists of names.
       Returns tuple (NEW_NET, FUNC).  FUNC is a function
         that will convert arrivals to the new NETWORK.
      Non-destructive."""
    
    q0l = map(lambda n: net.qid_of_queue (net.queue_by_name (n)), qnames0)
    q1l = map(lambda n: net.qid_of_queue (net.queue_by_name (n)), qnames1)

    sname_map_new = dict()

    # update the FSM.  This is the tricky part
    fsm_old = net.fsm
    A,O,sid2 = _split_state (fsm_old.a, fsm_old.o, net.sname2id[sname_old], q0l, q1l)
    fsm_new = hmm.HMM (A,O)
    
    sname_map_new = dict(net.sname2id)
    sid_old = sname_map_new [sname_old]
    del sname_map_new[sname_old]
    sname_map_new[sname_new1] = sid_old
    sname_map_new[sname_new2] = sid2

    net_new = qnet.Qnet (net.all_queues(), net.all_templates(), net.universe, dict(net.qname2id), sname_map_new, fsm_new)

    converter_fn = gen_arrv_converter (sid_old, sid_old, sid2, q0l, q1l)

    return net_new, converter_fn

def _split_state (a_old, o_old, state, q0l, q1l):
# use least-squares to find new parameters
# add new state to end
    ns = a_old.shape[0]
    no = o_old.shape[1]
    A = numpy.zeros((ns+1,ns+1))
    O = numpy.zeros((ns+1,no))
    A[0:ns,0:ns] = a_old
    O[0:ns,:] = o_old
    # fill in p ( * | snew2)
    A[ns,:] = A[state,:]
    # fill in p (q | snew2)
    O[ns,:] = O[state,:]  # order of these three lines matters
    O[ns,q0l] = 0
    w2 = numpy.sum (O[ns,:])
    O[ns,:] = O[ns,:] / w2
    # fill in p (q | snew1)
    O[state,q1l] = 0
    w1 = numpy.sum(O[state,:])
    O[state,:] = O[state,:] / w1
    # finally, fill in p( snew1 | *)  p ( snew2 | *)
    A[:,ns] = A[:,state]
    A[:,state] = w1 / (w1+w2) * A[:,state]
    A[:,ns] = w2 / (w1++w2) * A[:,ns]
    # done
    return A,O,ns

    # tricky part is the incoming state reallocate by solving linear system
    #  (we'll call the coefficients M)
#     M = numpy.zeros(( ns + no*ns,  2*(ns+1) ))
#     b = numpy.zeros( M.shape[0] )
#     ri = 0
#     def idx (s0, s1): return s0 if s1 == state else ns + s0
#     for si in range(ns):
#         M[ri,idx(si,state)] = M[ri,idx(si,ns)] = 1
#         b[ri] = a_old[si,state]
#         ri += 1
#     # q1 obs constraint
#     print M.shape
#     for si in range(ns):
#         for qi in range(no):
#             print "O[state,qi]", O[state,qi]
#             M[ri,idx(si,state)] = O[state,qi]
#             M[ri,idx(si,ns)] = O[si,qi]
#             print ri, "o_old", state, qi, o_old[state,qi]
#             print ri, b[ri]
#             b[ri] = o_old[state,qi]
#             ri += 1
#     # solve        
#     x,resids,rank,s = linalg.lstsq(M,b)
#     print "X", x
#     print "RESIDS", resids
#     A[:,state] = x[0:ns+1]
#     A[:,state] = A[:,state] / numpy.sum(A[:,state])
#     A[:,ns] = x[ns+1:]
#     A[:,ns] = A[:,ns]  / numpy.sum(A[:,ns])
#     # finally
#     return A,O,ns

def gen_arrv_converter (state_old, state1, state2, q0l, q1l):
    def converter (arrv):
        arrv_new = arrivals.Arrivals (arrv.qnet())
        for evt in arrv:
            e_new = evt.duplicate()
            if e_new.state == state_old:
                if e_new.qid in q0l:
                    e_new.set_state (state1)
                elif e_new.qid in q1l:
                    e_new.set_state (state2)
                else:
                    raise Exception ("Can't find %s in %s or %s" % evt, q0l, q1l)
            arrv_new.insert(e_new)
        # below is necessary for things like, events that have
        #  equal depature times whose orders get flipped.
        arrv_new.recompute_all_service() 
        return arrv_new
    return converter

def add_hidden_queue (net, arrv, sname_from, sname_to, qtext):
    # add state
    s1 = net.sname2id[sname_from]
    s2 = net.num_states()
    sname_map_new = dict(net.sname2id)
    sname_map_new[sname_to] = s2

    # add queue
    q2 = qnetu.queue_from_text (qtext, net.universe)
    qname_map_new = dict(net.qname2id)
    qname_map_new[q2.name] = len(net.qname2id)
    qlist = net.all_queues()[:]
    qlist.append (q2)
    tmpllist = net.all_templates()[:]
    tmpllist.append (q2.service)

    # update the FSM.  Only slightly tricky
    fsm_old = net.fsm
    ns = fsm_old.a.shape[0]
    no = fsm_old.o.shape[1]
    A = numpy.zeros((ns+1, ns+1))
    O = numpy.zeros((ns+1, no+1))
    A[0:ns,0:ns] = fsm_old.a
    A[s2,:] = A[s1,:]  # hidden state goes to same places as old state
    A[s1,:] = 0; A[s1,s2] = 1   # old state goes only to hidden
    O[0:ns,0:no] = fsm_old.o
    O[s2,no] = 1  # hidden state must output hidden queue
    fsm_new = hmm.HMM (A,O)

    net_new = qnet.Qnet (qlist, tmpllist, net.universe, qname_map_new, sname_map_new, fsm_new)

    # now... update the arrivals
    maxeid = arrv.max_eidx()+1
    eps = 1e-10 # hack; initialization will pick a better value
    arrv_new = arrivals.Arrivals (net_new)
    for evt in arrv:
        qnew = net_new.queue_by_name (evt.queue().name)  # currently unnecessary
        if evt.state == s1:
            first_evt = arrivals.Event (evt.eid, evt.tid, evt.qid, qnew, evt.a, evt.s, evt.d - eps, evt.state, obs_a=evt.obs_a, obs_d=0, proc=evt.proc, aux=evt.auxiliaries.duplicate())
            arrv_new.insert (first_evt)
            # insert hidden state
            qid_new = net_new.qid_of_queue(q2)
            arrv_new.insert (arrivals.Event (maxeid, evt.tid, qid_new, q2, evt.d - eps, eps, evt.d, s2, obs_a=0, obs_d=evt.obs_d))
            maxeid += 1
        else:
            arrv_new.insert (arrivals.Event (evt.eid, evt.tid, evt.qid, qnew, evt.a, evt.s, evt.d, evt.state, evt.obs_a, evt.obs_d, evt.proc, evt.auxiliaries.duplicate()))
    # recompute services b/c of substracting eps
    arrv_new.recompute_all_service()
            
    # done!
    return net_new, arrv_new


import pdb

#  Use the above for a unified model selection procedure.
#   IDEA: Report
def bottleneck_model_selection (net0, sname, arrv, qtext, outf=sys.stdout, gibbs_iter = 100, cluster_size=2, do_init=True, mdlIx=-1):
    tmr = mytime.timeit()

    if cluster_size != 2: raise NotImplementedError()  # need to add a choose() function here
    arrv0 = arrv.duplicate()
    if do_init: 
        arrv.recompute_all_service()
        net0.gibbs_initialize (arrv)
        arrv.validate()
    sid = net0.sid_by_name (sname)

    ix = 0


    qs = net0.queues_of_state (sid)
    for ix0 in range(len(qs)):
        for ix1 in range(ix0+1, len(qs)):
            ix += 1
            if mdlIx >= 0 and mdlIx != ix:
               print "MODEL_SELECTION :: SKIPPING %d (running %d)" % (ix, mdlIx)
               continue

            q0 = qs[ix0]
            q1 = qs[ix1]
            
            qlist = map(lambda q: q.name, qs)
            qlist.remove(q0.name)
            qlist.remove(q1.name)

            print "MODEL_SELECTION :: SPLITTING ", q0.name, q1.name

            sbad = "_%s_ABNORMAL" % sname
            sgood = "_%s_OK" % sname

            net1, converter = split_state (net0, sname, sbad, sgood, [q0.name, q1.name], qlist)
            print "INTERMEDIATE NET ", net1.as_yaml()
            arrv1 = converter(arrv)
            arrv1.validate()

            net2, arrv2 = add_hidden_queue (net1, arrv1, sbad, "HIDDEN", qtext)
            initial = net2.gibbs_initialize (arrv2)
            tmr.tick ("Initialization time [%d %d] " % (ix0, ix1))
            arrv2.validate()

            print "NET ", ix, net2.as_yaml()

            allmu = []
            def reporter (net, arrv, iter, lp):
                hdnq = net.queue_by_name ("HIDDEN")
                evts = arrv.events_of_queue (hdnq)
                allmu.append (numpy.mean ([ e.s for e in evts ]))

            estimation.bayes (net2, initial, gibbs_iter, report_fn = reporter)

            # output statistics
            mu_this = numpy.mean(allmu)
            std_this = numpy.std (allmu)
            print "HDN_Z %s__%s %.15f %.15f %.15f" % (q0.name, q1.name, mu_this, std_this, mu_this/std_this)

            f = open ("mmgmt_%s_%s.txt" % (q0.name, q1.name), "w")
            f.write ("\n".join(map(str, allmu)))
            f.close ()
            

