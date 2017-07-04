import unittest
import qnet
import qnetu
import numpy
import mytime
import netutils
import estimation
import yaml

import pwfun
import sampling
import arrivals
import qstats
import queues
import test_qnet

from scipy import integrate
from numpy import random

import gc
import arrivals

import sys
from math import sqrt

class TestMHGGk (unittest.TestCase):  

    def test_sem_via_gg1 (self):
        sampling.set_seed (342)

        nrep = 1
        nt = 250
        pct = 0.25
        
        neta = self.ln1a
        netb = self.ln1b
        
        for q in netb.queues:
            print "%s %f %f" % (q.name, q.service.mean(), q.service.std())
            
        all_mus = []
        for i in xrange(nrep):            
            sample_a = neta.sample (nt)
            observed = sample_a.subset_by_task (pct)
            arrv_a = neta.gibbs_initialize (observed.duplicate())
            arrv_b = qnetu.load_arrivals (netb, yaml.dump (arrv_a))
            arrv_a.validate()
            arrv_b.validate()
            
            print arrv_a
            print "NET A:\n", neta
            print "NET B:\n", netb
            
            # Sampling from net A (QueueGGk objects)
            mus_a,arrv_out = estimation.sem (neta, arrv_a, burn=50, num_iter=25, debug=False)
            avg_a = map(numpy.mean, zip(*mus_a))
            
            # Sampling from net B (Queue GG1 objects)
            mus_b,arrv_out = estimation.sem (netb, arrv_b, burn=50, num_iter=25, debug=False)
            avg_b = map(numpy.mean, zip(*mus_b))

            for q in netb.queues:
                print "%s %f %f" % (q.name, q.service.mean(), q.service.std())
            
            for thetas in zip(avg_a, avg_b): print "%.5f %.5f" % thetas
            for theta_a, theta_b in zip(avg_a, avg_b):
                self.assertAlmostEquals (theta_a, theta_b, 2)


    def test_sem_ggk (self):
        sampling.set_seed (342)
        nrep = 5
        nt = 100
        pct = 0.25        
        net = self.lnk
        self.do_test_sem (net, nrep, nt, pct, 1)

    def test_sem_mm3 (self):
        sampling.set_seed (342)
        nrep = 10
        nt = 250
        pct = 0.1
        net = self.mm3
        gibbs_iter = 1
#        qnet.set_sampler("ARS_PWFUN")
        self.do_test_sem (net, nrep, nt, pct, gibbs_iter, num_iter=250, do_init=True)

    def test_sem_mm020 (self):
        sampling.set_seed (342)
        nrep = 1
        nt = 100
        pct = 0.1
        net = self.mm020
        gibbs_iter = 1
#        qnet.set_sampler("ARS_PWFUN")
        self.do_test_sem (net, nrep, nt, pct, gibbs_iter=10, num_iter=250, do_init=True)

    def test_momsem_mm020 (self):
        sampling.set_seed (342)
        nrep = 1
        nt = 100
        pct = 0.1
        net = self.mm020
        gibbs_iter = 1
#        qnet.set_sampler("ARS_PWFUN")
        self.do_test_sem (net, nrep, nt, pct, gibbs_iter=1, num_iter=250, do_init=True, momentum=0.1)

    def test_gibbs_mm020 (self):
        sampling.set_seed (342)
        nrep = 1
        nt = 250
        pct = 0.1
        net = self.mm020

        arrv = net.sample (nt)
        sampled = arrv.subset_by_task (pct, adapt_fn=test_qnet.copy_evt)
        init = net.gibbs_initialize (sampled.duplicate())
    
        print "Initialization"
        print init

        def print_means (net, arrv, gi):
            print "GIBBS ITER", gi
            print "MEAN RESPONSE", " ".join (map(str,qstats.mean_response_time (arrv)))
            print "MEAN SERVICE", " ".join (map(str,qstats.mean_service (arrv)))
            print "MEAN SERVICE (OBS)", " ".join(map(str,qstats.mean_obs_service (arrv)))

        net.gibbs_resample (init, 0, 250, report_fn=print_means)


    def test_sem_gc (self):
        sampling.set_seed (342)
        gc.set_debug(gc.DEBUG_LEAK)
        nrep = 1
        nt = 250
        pct = 0.1
        net = self.mm020
        gibbs_iter = 1
        self.do_test_sem (net, nrep, nt, pct, gibbs_iter, num_iter=50, do_init=True)
#        self.do_test_bayes (net, nrep, nt, pct, gibbs_iter, num_iter=25, do_init=True)
        dump_garbage()


    def test_mm1_gc (self):
        sampling.set_seed (342)
        gc.set_debug(gc.DEBUG_LEAK)
        nrep = 1
        nt = 250
        pct = 0.1
        net = self.mm1
        gibbs_iter = 1
        self.do_test_sem (net, nrep, nt, pct, gibbs_iter, num_iter=50, do_init=True)
#        self.do_test_bayes (net, nrep, nt, pct, gibbs_iter, num_iter=25, do_init=True)
        dump_garbage()

    def test_gc (self):
        sampling.set_seed (342)
        gc.set_debug(gc.DEBUG_LEAK)
        arrv = self.mm020.sample (10)
        arrv = arrv.duplicate()
        dump_garbage()

    def test_gibbs_gc (self):
        sampling.set_seed (342)
        gc.set_debug(gc.DEBUG_LEAK)
        arrv = self.mm020.sample (10)
        smp = arrv.subset_by_task (0.5)
        arrvl = self.mm020.gibbs_resample(arrv, 0, num_iter=1000000, return_arrv = False)
        print arrvl[-1]
        arrvl = None
        dump_garbage()

    def do_test_sem (self, net, nrep, nt, pct, gibbs_iter, num_iter=100, do_init=True, momentum=0):
        mu0 = net.parameters[:]
        all_mus = []
        tot_theta = numpy.zeros(len(mu0))

        for i in xrange(nrep):         
            net.parameters = mu0[:]
            sample = net.sample (nt)
            observed = sample.subset_by_task (pct, adapt_fn=test_qnet.copy_evt)
            if do_init:
                arrv = net.gibbs_initialize (observed.duplicate())
            else:
                arrv = observed.duplicate()
            arrv.validate()
            
            print "Initialization:"
            print arrv
            print "Mean service (initialized): ", qstats.mean_service (arrv)
            print "Mean obs service (initialized): ", qstats.mean_obs_service (arrv)
            print "Mean service (true): ", qstats.mean_service (sample)
            mus,arrv_out = estimation.sem (net, arrv, burn=0, num_iter=num_iter, gibbs_iter=gibbs_iter, debug=False, return_arrv = False, momentum=momentum)
            avg = map(numpy.mean, zip(*mus))

#             for ai in range(len(arrv_out)):
#                 a = arrv_out[ai]
#                 f = open ("sem%d_%d.txt" % (i, ai), "w")
#                 qstats.write_arrv (f, a)
#                 f.close()

            for q in net.queues:
                print "%s %f %f" % (q.name, q.service.mean(), q.service.std())
            
            print "WAIT"
            w_gold = qstats.mean_wait (sample)
            w_last = qstats.mean_wait (arrv_out[-1])            
            for ws in zip (w_gold, w_last): print "%.5f %.5f" % ws

            print "RHO"
            rho_gold = qstats.utilization (sample)
            rho_last = qstats.utilization (arrv_out[-1])            
            for ws in zip (rho_gold, rho_last): print "%.5f %.5f" % ws


            print "THETA"
            for thetas in zip(avg, mu0): print "%.5f %.5f" % thetas

            tot_theta += numpy.array (avg)
#            for theta_a, theta_b in zip(avg, mu0):
#                self.assertAlmostEquals (theta_a, theta_b, 2)

        print "AVG_THETA"
        tot_theta /= nrep
        for thetas in zip(tot_theta, mu0):
            print "%.5f %.5f" % thetas


    def test_bayes_mm020 (self):
        sampling.set_seed (342)
        nrep = 1
        nt = 250
        pct = 0.1
        net = self.mm020
        gibbs_iter = 1
#        qnet.set_sampler("ARS_PWFUN")
        self.do_test_bayes (net, nrep, nt, pct, gibbs_iter, num_iter=5000, do_init=True)

    def test_bayes_mm020b (self):
        sampling.set_seed (342)
        nrep = 1
        nt = 250
        pct = 0.1
        net = self.mm020b
        gibbs_iter = 1
#        qnet.set_sampler("ARS_PWFUN")
        queues.set_proposal(queues.UNIFORM_PROPOSAL)
        self.do_test_bayes (net, nrep, nt, pct, gibbs_iter, num_iter=10000, do_init=True, do_output=0)

    def do_test_bayes (self, net, nrep, nt, pct, gibbs_iter, num_iter=100, do_init=True, do_output=0):
        mu0 = net.parameters[:]
        all_mus = []
        tot_theta = numpy.zeros(len(mu0))

        def do_report (net, arrv, i):
            if do_output > 0 and (i % do_output) == 0:
                f = open ("arrv_bayes_%d.txt" % i, "w")
                qstats.write_arrv (f, arrv)
                f.close()

        for i in xrange(nrep):         
            net.parameters = mu0[:]
            sample = net.sample (nt)
            observed = sample.subset_by_task (pct, adapt_fn=test_qnet.copy_evt)
            if do_init:
                arrv = net.gibbs_initialize (observed.duplicate())
            else:
                arrv = observed.duplicate()
            arrv.validate()
            
            print "Initialization:"
            print arrv
            print "Mean service (initialized): ", qstats.mean_service (arrv)
            print "Mean obs service (initialized): ", qstats.mean_obs_service (arrv)
            print "Mean service (true): ", qstats.mean_service (sample)
            mus,arrvl = estimation.bayes (net, arrv, num_iter=num_iter, return_arrv=False, report_fn=do_report)
            arrv_out = arrvl[-1]
            avg = map(numpy.mean, zip(*mus))

            for q in net.queues:
                print "%s %f %f" % (q.name, q.service.mean(), q.service.std())
            
            print "WAIT"
            w_gold = qstats.mean_wait (sample)
            w_last = qstats.mean_wait (arrv_out)
            for ws in zip (w_gold, w_last): print "%.5f %.5f" % ws

            print "RHO"
            rho_gold = qstats.utilization (sample)
            rho_last = qstats.utilization (arrv_out)
            for ws in zip (rho_gold, rho_last): print "%.5f %.5f" % ws


            print "THETA"
            for thetas in zip(avg, mu0): print "%.5f %.5f" % thetas

            tot_theta += numpy.array (avg)
#            for theta_a, theta_b in zip(avg, mu0):
#                self.assertAlmostEquals (theta_a, theta_b, 2)

        print "AVG_THETA"
        tot_theta /= nrep
        for thetas in zip(tot_theta, mu0):
            print "%.5f %.5f" % thetas

    def test_mm3_wait (self):
        net = self.mm3
        net.parameters = [ 30.0, 30.0, 100.0 ]
        arrv = self.mm3.sample (2500)
        w = qstats.mean_wait (arrv)
        print w
        self.assertTrue (w[1] < 2.0, "Weird waiting times: %s" % w)

    def test_mm1_stationary (self):
        net = self.mm1
        print net.parameters
        net.parameters = [ 30.0, 40.0, 40.0, 20.0 ]
        nt = 100
        nreps = 10
        giter = 250
        pct = 0.2
        self.check_param_stationary (net, nt, nreps, giter, pct)

    def test_mm3_stationary (self):
        sampling.set_seed (321223)
        net = self.mm3
#        net.parameters = [ 10.0, 20.0, 30.0 ]
        nt = 5
        nreps = 500
        giter = 50
        pct = 0.0
        qnet.set_sampler ("ARS_PWFUN")
#        qnet.set_sampler ("SLICE_PWFUN")
#        qnet.set_sampler ("SLICE")
        self.check_param_stationary (net, nt, nreps, giter, pct)

    def test_mm020_stationary (self):
        sampling.set_seed (321223)
        net = self.mm020
#        net.parameters = [ 10.0, 20.0, 30.0 ]
        nt = 100
        nreps = 50
        giter = 50
        pct = 0.0
#        qnet.set_sampler ("ARS_PWFUN")
#        qnet.set_sampler ("SLICE_PWFUN")
#        qnet.set_sampler ("SLICE")
        self.check_param_stationary (net, nt, nreps, giter, pct)

    def test_mm020b_stationary (self):
        sampling.set_seed (321223)
        net = self.mm020b
#        net.parameters = [ 10.0, 20.0, 30.0 ]
        nt = 250
        nreps = 50
        giter = 10
        pct = 0.0
#        qnet.set_sampler ("ARS_PWFUN")
#        qnet.set_sampler ("SLICE_PWFUN")
#        qnet.set_sampler ("SLICE")
        self.check_param_stationary (net, nt, nreps, giter, pct)

    def check_param_stationary (self, net, nt, nreps, giter, pct, do_output=0):
        theta0 = net.parameters[:]
        tot = numpy.zeros(len(net.parameters))
        tot_gibbs = [ numpy.zeros(len(net.parameters)) ]

        Ng = [0]

        tmr = mytime.timeit()

        def do_report (net, arrv, i):
            net.estimate (arrv)
            print "MU ", i, " ".join (map(str, net.parameters))
            Ng[0] += 1
            tot_gibbs[0] += numpy.array (net.parameters)
            net.parameters = theta0[:]
            print "THETA", net.parameters
            print "UT ", " ".join (map(str, qstats.utilization (arrv)))
            if do_output > 0 and (i % do_output) == 0:
                f = open ("cps_%d.txt" % i, "w")
                qstats.write_arrv (f, arrv)
                f.close()

        for rep in range(nreps):
            net.parameters = theta0[:]

            arrv = net.sample (nt)
            net.estimate (arrv)
            tot += numpy.array(net.parameters)
            net.parameters = theta0
            print "GOLD"
            print arrv

            obs = arrv.subset_by_task (pct, adapt_fn=test_qnet.copy_evt)
            resampled = net.gibbs_resample (obs, burn=0, num_iter=giter, return_arrv=False, report_fn=do_report)


        tot = tot / nreps
        tot_gibbs = tot_gibbs[0] / Ng[0]

        st = qnet.stats()
        print st
        secs = tmr.total ("Gibbs time")
        print "Events per sec = ", (st['nevents'] / secs)

        print
        print "TOT_MC  ",
        for x in tot: print " %.10f" % x,
        print

        print "TOT_GIB ", 
        for x in tot_gibbs: print " %.10f" % x,
        print

    def test_mm1_response_stationary (self):
        sampling.set_seed(2334)

        net = self.mm1
        nt = 250
        pct = 0.1
        nreps = 10
        giter = 25

        self.test_response_stationarity (net, nt, pct, nreps, giter)

    def test_mm3_response_stationary (self):
        sampling.set_seed(2334)

        net = self.mm3
        nt = 250
        pct = 0.1
        nreps = 1
        giter = 25

        self.test_response_stationarity (net, nt, pct, nreps, giter)


    def test_mm3_multiinit_response_stationary (self):
        sampling.set_seed(2334)

        net = self.mm3
        nt = 50
        pct = 0.2
        nreps = 10
        giter = 100

        self.test_response_stationarity (net, nt, pct, nreps, giter, resample_gold=False)


    def test_ln1b_response_stationary (self):
        sampling.set_seed(2334)

        net = self.ln1b
        nt = 50
        pct = 0.2
        nreps = 10
        giter = 25

        self.test_response_stationarity (net, nt, pct, nreps, giter)


    def test_lnk_response_stationary (self):
        sampling.set_seed(2334)

        net = self.lnk
        nt = 250
        pct = 0.2
        nreps = 10
        giter = 25

#        qnet.set_use_rtp(1)
        self.test_response_stationarity (net, nt, pct, nreps, giter)

    def test_ln3tier_response_stationary (self):
        sampling.set_seed(2334)

        net = self.ln3tier
        nt = 1
        pct = 0.0
        nreps = 10000
        giter = 10

#        qnet.set_use_rtp(1)

        # BAR
#        queues.set_proposal(queues.BAR_PROPOSAL)
#        qnet.set_sampler ("ARS_PWFUN")
        # END BAR

        # ZOID
        queues.set_proposal(queues.ZOID_PROPOSAL)
        qnet.set_sampler ("ARS_PWFUN")
        # END ZOID
#        queues.set_proposal(queues.TWOS_PROPOSAL)

        self.test_response_stationarity (net, nt, pct, nreps, giter, final_only=True)


    def test_response_stationarity (self, net, nt, pct, nreps, giter, resample_gold=True, do_init=False, do_report=False, final_only=False):
        params = net.parameters[:]

        tot_gold_mrt = numpy.zeros (net.num_queues())
        tot_sampled_mrt = numpy.zeros (net.num_queues())
        tot_gold_mst = numpy.zeros (net.num_queues())
        tot_sampled_mst = numpy.zeros (net.num_queues())

        arrv = None

        for rep in range(nreps):
            net.parameters = params[:]
            if not arrv or resample_gold:
                arrv = net.sample(nt)
                print "GOLD ARRIVALS"
                print arrv

                obs = arrv.subset_by_task (pct, adapt_fn=test_qnet.copy_evt)
                
                if final_only:
                    for e in obs:
                        if not obs.is_final(e):
                            e.obs = True

                if do_report:
                    f = open ("tmh_gold_%d.txt" % rep, "w")
                    qstats.write_arrv (f, obs)
                    f.close()

            tot_gold_mrt += qstats.mean_response_time (arrv)
            tot_gold_mst += qstats.mean_service (arrv)

            if do_init:
                initial = net.gibbs_initialize (obs.duplicate())
                print "INITIALIZATION"
                print initial
                initial.validate()
            else:
                initial = obs.duplicate()

            def do_write_arrv_text (net, arrv, iter):
                if do_report:
                    f = open ("tmh_arrv_%d_%d.txt" % (rep, iter), "w")
                    qstats.write_arrv (f, arrv)
                    f.close ()
                
            resampled = net.gibbs_resample(initial, burn=0, num_iter=giter, return_arrv=True, report_fn=do_write_arrv_text)
            arrvl = resampled

#            resampled = net.gibbs_resample (initial, burn=0, num_iter=giter, return_arrv=True)
#            mus, arrvl = estimation.sem (net, initial, burn=0, num_iter=giter)
#            for i in xrange(len(arrvl)):
#                print "ARRIVALS %d" % i
#                print arrvl[i]

            print "FINAL"
            print arrvl[-1]

            print "PARAMS: ", net.parameters

            print "SERVICE:"
            print "GOLD (OBS) ", qstats.mean_obs_service (obs)
            print "GOLD       ", qstats.mean_service (arrv)
            print "RESAMPLED  ", qstats.mean_service (arrvl[-1])
            print "RESPONSE:"
            print "GOLD (OBS) ", qstats.mean_obs_response (obs)
            print "GOLD       ", qstats.mean_response_time (arrv)
            print "RESAMPLED  ", qstats.mean_response_time (arrvl[-1])

            print "MU"
            for i in xrange(len(arrvl)):
                print i, qstats.mean_service (arrvl[i])

#            resampled = net.gibbs_resample (arrvl[-1], burn=0, num_iter=giter, return_arrv=True)
            tot_sampled_mrt += qstats.mean_response_time (resampled[-1])
            tot_sampled_mst += qstats.mean_service (resampled[-1])
            print resampled[-1]
            
            print qnet.stats()

        print "AVG_MRT_GOLD    ", tot_gold_mrt / nreps
        print "AVG_MRT_SAMPLED ", tot_sampled_mrt / nreps
        print "AVG_MST_GOLD    ", tot_gold_mst / nreps
        print "AVG_MST_SAMPLED ", tot_sampled_mst / nreps
    


    def test_ars_vs_slice (self):
        sampling.set_seed (178)

        net = self.mm3
        nt = 5
        nreps = 1000
        arrv = net.sample (nt)

        eid = 6
        eid_next = 13

        e = arrv.event (eid)
        e_next = arrv.event (eid_next)
        e.obs = 0
        e_next.obs = 0
        print arrv

        s_ars = []
        s_slice = []

        qnet.set_sampler ("ARS_PWFUN")
        for rep in range(nreps): 
            arrvl = net.gibbs_resample (arrv.duplicate(), burn=0, num_iter=1, return_arrv=True)
            s_ars.append (arrvl[-1].event(eid).d)
        print "STATS (ARS)"
        print qnet.stats()
        qnet.reset_stats()

        qnet.set_sampler ("SLICE_PWFUN")
        for rep in range(nreps): 
            arrvl = net.gibbs_resample (arrv.duplicate(), burn=0, num_iter=1, return_arrv=True)
            s_slice.append (arrvl[-1].event(eid).d)
        print "STATS (SLICE)"
        print qnet.stats()

        # for visualization
        d_fn = e.queue().departure_proposal (arrv, e)
        a_fn = e_next.queue().arrival_proposal (arrv, e_next)
        fn = pwfun.add (a_fn, d_fn)
        print e
        print e_next
        print d_fn.dump_table()
        print a_fn.dump_table()
        print fn.dump_table()
 
        limits = fn.range()
        rng = limits[1] - limits[0]
        if numpy.isinf (rng): rng = 100.0
        x = limits[0]
        for i in xrange(1000):
            print "FN %.5f %.5f" % (x, fn(x))
            x += rng/1000

        s_ars.sort()
        s_slice.sort()
        for d_ars, d_sl in zip(s_ars, s_slice): 
            print "D", d_ars, d_sl

        netutils.check_quantiles (self, s_ars, s_slice, nreps)
        
    def test_single_departure (self):
        sampling.set_seed (178)

        net = self.mm3
        nt = 5
        nreps = 100

        d_ex = []
        d_gibbs = []

        giter = 1

        for rep in xrange(nreps):
            arrv = net.sample (nt)
            evts = arrv.events_of_task (2)
            e = evts[1]
            e_next = evts[2]
            d_ex.append (e.d)
            e.obs = 0 
            e_next.obs = 0
            
            for i in xrange(giter):
                arrvl = net.gibbs_resample (arrv, burn=0, num_iter=1)
                arrv = arrvl[-1]
                evts = arrv.events_of_task (2)
                e = evts[1]
                e_next = evts[2]
                d_gibbs.append (e.d)

        netutils.check_quantiles (self, d_ex, d_gibbs, nreps)

        print qnet.stats()

    def test_ln3tier_via_slice (self):
        sampling.set_seed(2334)
        net = self.ln3tier
        nt = 250
        self.do_test_via_slice(net, nt, lambda arrv: arrv.events_of_qid(1)[0])

    def test_ln3tier_zoid_slice (self):
        sampling.set_seed(2334)
        net = self.ln3tier
        nt = 250

        queues.set_proposal(queues.ZOID_PROPOSAL)
        qnet.set_sampler ("ARS_PWFUN")

        self.do_test_via_slice(net, nt, lambda arrv: arrv.events_of_qid(1)[0])

    def test_mm020_slice (self):
        sampling.set_seed(2334)
        net = self.mm020
        nt = 250
        self.do_test_via_slice(net, nt, lambda arrv: arrv.events_of_qid(1)[0])

    def do_test_via_slice (self, net, nt, grab_evt):
        arrv = net.sample (nt)
        evt5 = grab_evt(arrv)
        evt5n = evt5.next_by_task()

        evt5.obs = False
        evt5n.obs = False
        
        print arrv
        print "Events selected:\n %s\n %s" % (evt5, evt5n)
        
        N = 100

        alld = []
        for i in xrange (N):
            arrvl = net.gibbs_resample (arrv, burn=0, num_iter=1)
            evt_rs = arrvl[-1].event (evt5.eid)
            alld.append (evt_rs.d)

        d_mean = numpy.mean(alld)

        prp = queues.pair_proposal (arrv, evt5, evt5n)
        a,b = prp.range()

        xs = list(random.uniform (a, b, size=N))
        ws = []
        for x in xs:
            dl = evt5.queue().pyDiffListForDeparture (evt5, x)
            dl.extend (evt5n.queue().pyDiffListForArrival (evt5n, x))
            arrv.applyDiffList (dl)
            evt5.d = x
            evt5n.a = x
            ws.append(net.log_prob (arrv) - numpy.log(b-a))
            
        Z = logsumexp(ws)
        Ws = [ numpy.exp(w - Z) for w in ws ]

        print "Z = ", Z

        sum = 0
        for w,x in zip(Ws,xs):
            sum += w*x
        mc_mean = sum

        print qnet.stats()
        print "Sample mean = ", d_mean
        print "MCINT mean = ", mc_mean

        f = open("prp-zoid.txt", "w")
        L,U = prp.range()
        eps = (U-L)/1000
        x = L
        while x < U:
            f.write ("%.5f %.5f\n" % (x, prp(x)))
            x += eps
        f.close()

        dfn = evt5.queue().pyDepartureLik (arrv, evt5)
        afn = evt5n.queue().pyArrivalLik (arrv, evt5n)
        product = pwfun.add(dfn, afn)

        f = open ("prp-product.txt", "w")
        L,U = product.range()
        eps = (U-L)/1000
        x = L
        while x < U:
            f.write ("%.5f %.5f\n" % (x, product(x)))
            x += eps
        f.close()

        # ZOID DEBUGGING
        for x in prp.knots():
            print "%.5f %.5f %.5f %.5f" % (x, prp(x), product(x), dfn(x)+afn(x)) 

        self.assertTrue (abs(d_mean - mc_mean) < 0.1, "Not close enough")

    def test_logprob (self):
        arrvd = self.delay.sample (50)
        qnetu.write_multif_to_prefix ("tlp.", arrvd)
        arrvk = qnetu.read_multif_of_prefix ("tlp.", self.delay)
        lpd = self.delay.log_prob (arrvd)
        lpk = self.k100.log_prob (arrvk)
        self.assertTrue ( abs(lpd-lpk) < 1e-5 )
        
    mm3_text = """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: TIER1
  queues: [ TIER1_0 , TIER1_1 ]
  successors: [TIER2]
- name: TIER2
  queues: [ TIER2_0 , TIER2_1 ]
  successors: [TIER3]
- name: TIER3
  queues: [ TIER3_0 ]
queues:
- { name: INITIAL, service: [M, 1.0 ] }
- { name: TIER1_0, processors: 5, service: [M, 6] }
- { name: TIER1_1, processors: 5, service: [M, 6] }
- { name: TIER2_0, processors: 2, service: [M, 3] }
- { name: TIER2_1, processors: 2, service: [M, 3] }
- { name: TIER3_0, processors: 3, service: [M, 2] }
"""

    mm020_text = """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: TIER1
  queues: [ TIER1_0 ]
  successors: [TIER2]
- name: TIER2
  queues: [ TIER2_0 ]
queues:
- { name: INITIAL, service: [M, 0.5 ] }
- { name: TIER1_0, processors: 100, service: [M, 1.0] }
- { name: TIER2_0, processors: 100, service: [M, 1.0] }
"""
#- { name: TIER2_0, processors: 5, service: [M, 1.0] }

#- { name: INITIAL, service: [M, 2.0 ] }
#- { name: TIER1_0, processors: 5, service: [M, 2] }
#- { name: TIER2_0, processors: 2, service: [M, 0.8] }

    mm020b_text = """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: TIER1
  queues: [ TIER1_0 , TIER1_1 ]
  successors: [TIER2]
- name: TIER2
  queues: [ TIER2_0 , TIER2_1 ]
  successors: [TIER3]
- name: TIER3
  queues: [ TIER3_0 ]
queues:
- { name: INITIAL, service: [M, 1.0 ] }
- { name: TIER1_0, processors: 5, service: [M, 2] }
- { name: TIER1_1, processors: 5, service: [M, 2] }
- { name: TIER2_0, processors: 2, service: [M, 0.8] }
- { name: TIER2_1, processors: 2, service: [M, 0.8] }
- { name: TIER3_0, processors: 3, service: [M, 0.6] }
"""
    mm3_text_foo = """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: TIER1
  queues: [ TIER1_0  ]
  successors: [TIER2]
- name: TIER2
  queues: [ TIER2_0 ]
queues:
- { name: INITIAL, service: [M, 2.0 ] }
- { name: TIER1_0, processors: 5, service: [M, 6] }
- { name: TIER2_0, processors: 2, service: [M, 3] }
"""
#- { name: INITIAL, service: [M, 2.0 ] }
#- { name: INITIAL, type: GGk, processors: 1, service: [M, 2.0 ] }
    
    
    ln1a_text = """
    states:
      - name: INITIAL
        queues: [ INITIAL ]
        successors: [ TIER1 ]
        initial: TRUE 
      - name: TIER1
        queues: [ WEB1 ]
        successors: [ TIER2 ]
      - name: TIER2
        queues: [ APP1 ]
    queues:
      - { name: INITIAL, service: [M, 10.0]  }
      - { name: WEB1, processors: 1, service: [LN, 1.75, 1.0] }
      - { name: APP1, processors: 1, service: [LN, 0.1, 2.0] }
    """
    
    ln1a_text = """
    states:
      - name: INITIAL
        queues: [ INITIAL ]
        successors: [ TIER1 ]
        initial: TRUE 
      - name: TIER1
        queues: [ WEB1 ]
        successors: [ TIER2 ]
      - name: TIER2
        queues: [ APP1 ]
    queues:
      - { name: INITIAL, service: [M, 10.0]  }
      - { name: WEB1, processors: 1, type: GGk, service: [LN, 1.75, 1.0] }
      - { name: APP1, processors: 1, type: GGk, service: [LN, 0.1, 2.0] }
    """
    
    ln1b_text = """
    states:
      - name: INITIAL
        queues: [ INITIAL ]
        successors: [ TIER1 ]
        initial: TRUE 
      - name: TIER1
        queues: [ WEB1 ]
        successors: [ TIER2 ]
      - name: TIER2
        queues: [ APP1 ]
    queues:
      - { name: INITIAL, service: [M, 15.0]  }
      - { name: WEB1, processors: 1, type: GG1, service: [LN, 1.75, 1.0] }
      - { name: APP1, processors: 1, type: GG1, service: [LN, 0.1, 2.0] }
    """
    
    lnk_text = """
    states:
      - name: INITIAL
        queues: [ INITIAL ]
        successors: [ TIER1 ]
        initial: TRUE 
      - name: TIER1
        queues: [ WEB1 ]
        successors: [ TIER2 ]
      - name: TIER2
        queues: [ APP1 ]
    queues:
      - { name: INITIAL, service: [M, 5.0]  }
      - { name: WEB1, processors: 3, type: GGk, service: [LN, 1.75, 1.0] }
      - { name: APP1, processors: 3, type: GGk, service: [LN, 0.1, 2.0] }
    """

    ln3tier_text = """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: TIER1
  queues: [ TIER1_0 ]
  successors: [TIER2]
- name: TIER2
  queues: [ TIER2_0 ]
queues:
- { name: INITIAL, service: [G, 2, 10.0] }
- { name: TIER1_0, service: [LN, 0.9, 0.5 ] }
- { name: TIER2_0, service: [LN, 0.9, 0.5 ] }

"""

    foo_text = """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: TIER1
  queues: [ TIER1_0 , TIER1_1 ]
  successors: [TIER2]
- name: TIER2
  queues: [ TIER2_0 , TIER2_1 ]
  successors: [TIER3]
- name: TIER3
  queues: [ TIER3_0 , TIER3_1 ]
queues:
- { name: INITIAL, service: [G, 2, 1.0] }
- { name: TIER1_0, service: [LN, 0.9, 0.5 ] }
- { name: TIER1_1, service: [LN, 0.9, 0.5 ] }
- { name: TIER2_0, service: [LN, 0.9, 0.5 ] }
- { name: TIER2_1, service: [LN, 0.9, 0.5 ] }
- { name: TIER3_0, service: [LN, 0.9, 0.5 ] }
- { name: TIER3_1, service: [LN, 0.9, 0.5 ] }

"""

    mm1_text ="""
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: TIER1
  queues: [ TIER1_0, TIER1_1] 
queues:
- { name: INITIAL, service: [M, 0.02] }
- { name: TIER1_0, service: [M, 0.038] }
- { name: TIER1_1, service: [M, 0.038] }
"""

    mm1_text = \
    """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: TIER1
  queues: [ TIER1_0, TIER1_1] 
  successors: [TIER2]
- name: TIER2
  queues: [TIER2_0]
queues:
- { name: INITIAL, service: [M, 0.02] }
- { name: TIER1_0, service: [M, 0.035] }
- { name: TIER1_1, service: [M, 0.035] }
- { name: TIER2_0, service: [M, 0.015] }
"""


    k100_text = \
    """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: TIER1
  queues: [ TIER1_0 ]
queues:
- { name: INITIAL, service: [M, 1.0] }
- { name: TIER1_0, service: [M, 0.5], processors: 100 }
    """

    delay_text = \
    """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: TIER1
  queues: [ TIER1_0 ]
queues:
- { name: INITIAL, service: [M, 1.0] }
- { name: TIER1_0, service: [M, 0.5], type: DELAY }
"""

    def setUp (self):
        self.mm1 = qnetu.qnet_from_text (TestMHGGk.mm1_text)
        self.mm3 = qnetu.qnet_from_text (TestMHGGk.mm3_text)
        self.mm020 = qnetu.qnet_from_text (TestMHGGk.mm020_text)
        self.mm020b = qnetu.qnet_from_text (TestMHGGk.mm020b_text)
        self.ln1a = qnetu.qnet_from_text (TestMHGGk.ln1a_text)
        self.ln1b = qnetu.qnet_from_text (TestMHGGk.ln1b_text)    
        self.ln3tier = qnetu.qnet_from_text (TestMHGGk.ln3tier_text)    
        self.lnk = qnetu.qnet_from_text (TestMHGGk.lnk_text)
        self.k100 = qnetu.qnet_from_text (TestMHGGk.k100_text)
        self.delay = qnetu.qnet_from_text (TestMHGGk.delay_text)


def main():
    if len(sys.argv) > 1:
        for test_name in sys.argv[1:]:
            suite = unittest.TestLoader().loadTestsFromName("test_mh_ggk.TestMHGGk.%s" % (test_name,))
            unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()

def print_with_type(x):
    s = str(x)
    if len(s) > 80: s = s[:120]
    print type(x),"\n  ", s

def dump_garbage():
    print "\nGARBAGE:"
    gc.collect()

    print "\nGARBAGE OBJECTS:"
    for x in gc.garbage:
        print_with_type (x)
        if isinstance(x, arrivals.Event):
            print "REFERS:"
            for y in gc.get_referrers (x):
                print_with_type (y)
            print "End refers"

def logsumexp (args):
    x_max = max (args)
    if x_max == -numpy.inf: return x_max
    return x_max + numpy.log (sum ([ numpy.exp(x - x_max) for x in args ]))

