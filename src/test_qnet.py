import unittest
import qnet
import qnetu
import numpy
import numpy.random
import mytime
import netutils
import estimation
import yaml
import arrivals
import qstats
import pwfun
import distributions

import sys
import sampling
from math import sqrt

class TestQnet (unittest.TestCase):  

    def test_read_from_text (self):
        self.assertEquals (7, self.qnet3.num_queues ())
        self.assertEquals (4, self.qnet3.num_states ())
        
        initial = self.qnet3.sid_by_name ("INITIAL")
        tier1 = self.qnet3.sid_by_name ("TIER1")
        tier2 = self.qnet3.sid_by_name ("TIER2")
        tier3 = self.qnet3.sid_by_name ("TIER3")
        
        fsm = self.qnet3.get_fsm()
        self.assertEquals (set([ tier1 ]), set(fsm.successors (initial)))
        self.assertEquals (set([ tier2 ]), set(fsm.successors (tier1)))
        self.assertEquals (set([ tier3 ]), set(fsm.successors (tier2)))
        self.assertEquals (set([]), set(fsm.successors (tier3)))

    def test_sample (self):
        sampling.set_seed (1897239)

        N = 50
        arrvls = self.qnet3.sample (N)
        arrvls.validate ()
        
        s_tot = 0
        s2_tot = 0
        
        web2_tot = 0
        app1_tot = 0
        db_tot = 0
        
        web2 = self.qnet3.queue_by_name ("WEB2")
        app1 = self.qnet3.queue_by_name ("APP1")
        db = self.qnet3.queue_by_name ("DB")
        
        for e in arrvls.all_events():
            s_tot += e.service()
            s2_tot += e.service()*e.service()
            q = e.queue()
            
            if q == web2:
                web2_tot += 1
            elif q == app1:
                app1_tot += 1
            elif q == db:
                db_tot += 1
                
#        print "%d %d %d" % (web2_tot, app1_tot, db_tot)
#        print "%f %f" % (s_tot, s2_tot)
        
        self.assertEquals (web2_tot, len(arrvls.events_of_queue (web2)))
        self.assertEquals (app1_tot, len(arrvls.events_of_queue (app1)))
        self.assertEquals (N, len(arrvls.events_of_queue (db)))
        
        ne = len(arrvls.all_events())
        
        self.assertTrue (abs(s_tot/ne - 1.0) < 2.0/sqrt(ne), "Average S wrong, was %.10f expected 1.0 (STD %.10f)" % (s_tot/ne, 1/sqrt(ne)))
        self.assertTrue (abs(web2_tot - N/3) < 2*N*0.333*0.6666)  # 2 SD
        self.assertTrue (abs(app1_tot - N/2) < 2*N*0.25)  # 2 SD
        self.assertEquals (db_tot, N)
            

    def test_mm3_sample (self):
        N = 25
        arrvls = self.mm3.sample (N)
        print arrvls
        arrvls.validate ()
        
        evts = arrvls.all_events()
        evts[55].d += 1
        
        was_except = False
        try:
            arrvls.validate()
        except AssertionError:
            was_except = True # exception expected
            
        self.assertTrue (was_except)

    def test_sample_quantile (self):
        sampling.set_seed (1897239)
        theta = self.qnet3.parameters
        theta[0] = 2.0
        self.qnet3.parameters = theta

        Asum = []
        nreps = 1000
        for i in range(nreps):
            arrvls = self.qnet3.sample (10)
            maxA = max(e.d for e in arrvls if e.state == 0)
            Asum.append (maxA)

        Asum.sort()

        self.assertTrue (abs (25.0 - Asum[800]) < 0.1)
        
        
    def test_sample_by_time (self):
        sampling.set_seed (28922980)
        theta = self.qnet3.parameters
        theta[0] = 2.0
        self.qnet3.parameters = theta
        
        nreps = 1000
        T = 30.0

        d1 = []
        d5 = []
        N = []
        for i in range(nreps):
            arrvls = self.qnet3.sample_by_time (T)
            arrvls.validate ()
            N.append (arrvls.num_tasks())
            d1.append (arrvls.event(0).d)
            d5.append (arrvls.event(4).d)
            maxA = max(e.d for e in arrvls if e.state == 0) 
            self.assertTrue (maxA < T, "maxA too big, was %.5f should have been < %.4f (N was %s)\n%s" % (maxA, T, N, arrvls))

        Nbar = numpy.mean(N)
        self.assertTrue (abs (Nbar - 15) < 3*sqrt(2),
                         "Mean ntasks %.4f expected 15" % Nbar)

        d1.sort()
        d1mean = numpy.mean(d1)
        self.assertTrue (abs(d1mean - 2.0) < 0.1, "Huh? %.5f" % d1mean)
        print d1mean
        
        d5.sort()
        d5mean = numpy.mean(d5)
        d5q7 = d5[int(0.7 * nreps)]

        print d5mean, numpy.var (d5)
        self.assertTrue (abs(d5mean - 10.0) < 0.25)
        self.assertTrue (abs(d5q7 - 11.78) < 0.25,
                         "Bad 0.7 quantile %.5f" % d5q7)  # qgamma (0.7, shape=5, rate=0.5)


    # TODO: Segfault, appears to be d_proc_prev not being initialized            
    def fixme_test_gibbs_final (self):
        sampling.set_seed (12312)
        filter_fn = lambda arrv,evt: not arrv.is_final (evt)
        print self.oneq.parameters
        self.check_gibbs_stationarity (ns=250, nt=5, num_giter=25, net=self.oneq, filter_fn=filter_fn, corrupt_fn=zero_departure)


    # TODO: Segfault, appears to be d_proc_prev not being initialized
    def fixme_test_gibbs_final2 (self):
        net = self.oneq
        self.assertEquals (2, len(net.parameters))
        net.parameters = [ 1.0, 0.1 ]
        filter_fn = lambda arrv,evt: not arrv.is_final (evt)
        self.check_gibbs_stationarity (ns=500, nt=5, num_giter=3, net=net, filter_fn=filter_fn, corrupt_fn=copy_evt)


    def check_gibbs_stationarity (self, ns=1, nt=1, num_giter=1, net=None, filter_fn=None, corrupt_fn = None):
        sampled = [ net.sample (nt) for i in xrange(ns) ]
        avg_times = self.mean_event_times (sampled)
        
        tmr = mytime.timeit()
        
        all_resampled = []
        for i in xrange(ns):
            arrvls = sampled[i]
            arrv_obs = arrvls.subset (filter_fn, corrupt_fn) 
            print arrv_obs
            resampled = net.gibbs_resample (arrv_obs, 0, num_giter)
            all_resampled.append (resampled[-1])
        elapsed = tmr.total()


#        for i in xrange(len(sampled)):
#            print "ORIG  %s\nGIBBS %s\n=====" % (sampled[i].event(3), all_resampled[i].event(3))
        
        gibbs_times = self.mean_event_times (all_resampled)
        
        smpl0 = all_resampled[0]
        final_eidx = [ i for i in xrange(smpl0.num_events()) if not filter_fn (smpl0, smpl0.event(i)) ]
        
        print "Events resampled per sec = ", (len(final_eidx) * ns * num_giter) / elapsed
        
        self.assertEquals (len(avg_times), len(gibbs_times))
        print avg_times
        print gibbs_times
        for i in xrange(len(avg_times)): print "%s  %s" % (avg_times[i], gibbs_times[i])
        for i in final_eidx:
            self.assertTrue (abs(avg_times[i] - gibbs_times[i]) > 1e-20)
            self.assertTrue (abs(avg_times[i] - gibbs_times[i]) / avg_times[i] < 0.05, \
                            "Mismatch event %d  Sampled: %.4f  Gibbs: %.4f  NormDiff: %.4f" % (i, avg_times[i], gibbs_times[i], abs(avg_times[i] - gibbs_times[i]) / avg_times[i]))


    def mean_event_times (self, arrv_list):
        ne = arrv_list[0].num_events()
        tot = numpy.zeros ((ne,))
        
        for arrv in arrv_list:
            times = numpy.array ([ evt.d for evt in arrv ])
            tot = tot + times
        
        tot = tot / len(arrv_list)
        
        return tot


    def test_m2_via_gamma (self):
        """Runs the gibbs sampler, checks that the sum of exponential service times is gamma."""
        nt = 5
        ns = 100
        num_giter = 10
        net = self.m2
        
        tmr = mytime.timeit()
        
        times = []
        for i in xrange(ns):
            arrv = net.sample (nt)
            obs = arrv.subset (lambda a,e: a.is_initial (e), copy_evt)
            gsmp = net.gibbs_resample (obs, 0, num_giter)
            resampled = gsmp[-1]
            for tid in xrange(nt):
                byt = resampled.events_of_task (tid)
                self.assertEquals (3, len(byt))
                times.append (byt[1].s + byt[2].s)
        
        elapsed = tmr.total()        
        print "Events resampled per sec = ", (nt * 2 * ns * num_giter) / elapsed
        
        s1 = numpy.array ([ numpy.random.exponential (0.5) for i in xrange(ns) ])
        s2 = numpy.array ([ numpy.random.exponential (0.25) for i in xrange(ns) ])
        s = s1+s2
        s.sort()
        times.sort()
        print "FROM QNET:\n", summarize(times)
        print "EXACT SAMPLE:\n", summarize(s)
        
        netutils.check_quantiles (self, s, times, ns)
        
        
    def test_m2b_via_uniform (self):
        """Runs the gibbs sampler, checks that the sum of exponential service times is gamma."""
        nt = 5
        ns = 1
        num_giter = 100
        net = self.m2b

        tmr = mytime.timeit()

        # For this test, each sample is tested independently rather than aggregated
        for i in xrange(ns):
            arrv = net.sample (nt)
            print arrv
            obs = arrv.subset (lambda a,e: a.is_initial (e), copy_evt)
            gsmp = net.gibbs_resample (obs, 0, num_giter, sample_final=False)
            for tid in xrange(nt):
                # For each task, check that the Gibbs distribution is correctly uniform
                times = []
                for smp_id in xrange(1,len(gsmp)):
                    byt = gsmp[smp_id].events_of_task (tid)
                    self.assertEquals (3, len(byt))
                    times.append (byt[1].d)
                
                # examine gibbs function
                e0 = arrv.events_of_task (tid)[1]
                e1 = arrv.events_of_task (tid)[2]
                L = e0.a
                U = e1.d
                cdist = net.gibbs_for_departure (obs, e0)
                xs = [ L+ i*(U-L)/10 for i in xrange(10) ]
                for x in xs: print "  x %.4f  p(d = x | A) %.4f" % (x, cdist(x))
                
                # generate true sample
                s = [ numpy.random.uniform (L, U) for i in xrange(num_giter) ]            

                # now check the cdfs
                s.sort()
                times.sort()
                print summarize (times)
                netutils.check_quantiles (self, s, times, num_giter)

        elapsed = tmr.total()        
        print "Events resampled per sec = ", (nt * 2 * ns * num_giter) / elapsed

    def test_poisson(self):
        """Tests that Gibbs sampling the initial process yields a Poisson process."""
        nt = 50
        ns = 1000
        num_giter = 5
        net = self.poisson

        times = []
        for i in range(ns):
            arrv = net.sample (nt)
            obs = arrv.subset (lambda a,e: a.is_last_in_queue(e), copy_evt)
            gsmp = net.gibbs_resample (arrv, 0, num_giter)
            resampled = gsmp[-1]
            evts = resampled.events_of_task (2)
            times.append (evts[0].d)
            
        exact_sample = [ numpy.random.gamma (shape=3, scale=0.5) for i in xrange (ns) ]
        times.sort()
        exact_sample.sort()
        
        print summarize(times)
        print summarize(exact_sample)
        
        netutils.check_quantiles (self, exact_sample, times, ns)

    def check_obs_equal (self, arrv, initialized):
        for evt in arrv:
            if evt.obs_a or evt.obs_d:
                new_evt = initialized.event (evt.eid)
                self.assertTrue (new_evt.obs_a == evt.obs_a, "No obs match\n  %s\n  %s" % (evt, new_evt))
                self.assertTrue (new_evt.obs_d == evt.obs_d, "No obs match\n  %s\n  %s" % (evt, new_evt))
                if evt.obs_d:
                    self.assertAlmostEquals (evt.d, new_evt.d, 8)
                if evt.obs_a:
                    self.assertAlmostEquals (evt.a, new_evt.a, 8)

    def test_initialize (self):
        print self.oneq
        sampling.set_seed (438)
        nt = 500
        for i in xrange(1):
            arrv = self.oneq.sample (nt)
            obs = arrv.subset_by_task (0.5)

#             print "ORIG"
#             print arrv
#             print "OBS"
#             print obs

            arrv.validate()
            obs.validate()

#            print obs
#            print qnet.create_qp (obs, 10)
            
            result = self.oneq.gibbs_initialize (obs)
            print "RESULT"
            print result
            
            result.validate()
            self.check_obs_equal (obs, result)
            
            self.assertEquals (len(obs.all_tasks()), len(result.all_tasks()))

    def test_initialize2 (self):
        nt = 250
        self.oneq.parameters = [ 0.1, 0.5 ]
        sampling.set_seed (339)
        
        arrv = self.oneq.sample (nt)
        obs = arrv.subset_by_task (0.5)
        result = self.oneq.gibbs_initialize (obs)
        result.validate()
    
    def test_initialize3 (self):
        
        N = 250
        pct = 0.5
        net = self.qnet3
        net.parameters = [1.0, 0.2, 0.2, 0.2, 0.05, 0.05, 0.3]

        sampling.set_seed (342)        
        for i in xrange(1):
            print "#",
            sample = net.sample (N)
            observed = sample.subset_by_task (pct)

            print "ORIG ", sample

#            print observed
#            print qnet.create_qp (observed, 10)

            try:

                arrv = net.gibbs_initialize (observed.duplicate())
                print "RESULT ", arrv
            
                arrv.validate()
                
            except AssertionError, e:
                print observed
                raise e
                
    def test_initialize4 (self):
        sampling.set_seed (1335)
        
        net = self.rails3
        arrv = self.rails3_arrv
        
        obs = arrv.subset_by_task (0.2)
        print obs
        
        result = net.gibbs_initialize (obs)

        result.validate()
        
    def test_mix_initialize (self):
        sampling.set_seed (484)
        nt = 100
        arrv = self.mix1.sample (nt)
        obs = arrv.subset_by_task (0.5)

        result = self.mix1.gibbs_initialize (obs)
        result.validate()
        
        N = [0, 0]
        Ntot = 0
        for evt in result:
            if evt.queue().name == "THE_Q":
                component = evt.variable("COMPONENT_THE_Q")
                N[component] += 1
                Ntot += 1

        print "RESULT"
        print result

        print "PROPORTIONS: ", N
        self.assertTrue (0.5 * abs(N[0] - N[1]) < 2*sqrt(0.5*0.5*Ntot), "Difference too big: was %d, expected no more than %s" % (0.5*abs(N[0]-N[1]), 2*sqrt(0.5*0.5*Ntot)))
        
        self.check_obs_equal (obs, result)
        

    def test_mm3_initialize (self):
         sampling.set_seed (1338)

         N = 30
         net = self.mm3
         arrv = net.sample (N)
         print "ORIGINAL"
         print arrv
         
         obs = arrv.subset_by_task (0.2)
         arrv = net.gibbs_initialize (obs)
     
         print "INITIALIZED"
         print arrv
         arrv.validate()

    def test_minilp_initialize (self):
        sampling.set_seed (1335)
        
        net = self.rails3
        arrv = net.sample(1000)
        print "GOLD ", arrv

        q = net.queue_by_name ("DB")
        evts = arrv.events_of_queue (q)

        for e in evts:
            e.obs_a = 0
            e.previous_by_task().obs_d = 0

#        print arrv
        
        result = qnet.gibbs_initialize_via_minilp (net, arrv.duplicate())

        result.validate()
#        print arrv
#        print result

        self.check_obs_equal (arrv, result)

    def test_initialize_spread (self):
        net = self.oneq
        init_old = qnet.get_initializer()
        qnet.set_initialization_type ("SPREAD")
        sampling.set_seed (438)
        nt = 500
        nt = 10
        
        for i in xrange(10):
            arrv = net.sample (nt)
            obs = arrv.subset_by_task (0.5)

#             print "ORIG"
#             print arrv
#             print "OBS"
#             print obs

            arrv.validate()
            obs.validate()

            result = self.oneq.gibbs_initialize (obs)
            print "RESULT"
            print result
            
            result.validate()
            self.check_obs_equal (obs, result)
            
            self.assertEquals (len(obs.all_tasks()), len(result.all_tasks()))

            print qstats.mean_response_time (arrv)
            print qstats.mean_response_time (result)
            
        # don't screw up other tests
        qnet.set_initialization_type (init_old)
        
    def test_inserting (self):
        sampling.set_seed (1492)
        N = 10
        reps = 10
        
        for ri in range(reps):
            arrv = self.qnet3.sample (N)
            obs = arrv.subset_by_task (0.5)
            nt_expected = 1+max(range(N), key=lambda tid: tid if obs.events_of_task(tid)[0].obs_d else -1)

            for tid in range(N):
                tevts = obs.events_of_task (tid)
                if not tevts[0].obs_d:
                    obs.delete_task (tid)
            result = qnet.initialize_inserting_implied_tasks (self.qnet3, obs)

            self.assertEquals (nt_expected, len(result.initial_events()))        
            self.assertEquals (4*nt_expected, len(result.all_events()))
         
    def test_subset_by_task (self):
        nt = 1000
        arrv = self.qnet3.sample (nt)
        obs = arrv.subset_by_task (0.4)
        
        tasks_in = {}
        tasks_out = {}
        
        for evt in obs:            
            if evt.obs_a or evt.obs_d:
                self.assertTrue (evt.tid not in tasks_out)
                tasks_in[evt.tid] = True
            else:
                self.assertTrue (evt.tid not in tasks_in)
                tasks_out[evt.tid] = True
            
        num_obs = sum([ 1 for d in tasks_in ])
        
        self.assertTrue (abs(num_obs/float(nt) - 0.4) <= 2*nt*sqrt(0.4*0.6))  # 2 sd


    def test_sampled_by_task (self):
        nt = 100
        
        for i in xrange(10):
            arrv = self.qnet3.sample (nt)
            obs = arrv.subset_by_task (0.4)
            self.assertTrue (qnet.is_sampled_by_task (obs))

        for i in xrange(10):
            arrv = self.qnet3.sample (nt)
            obs = arrv.subset (lambda a,e: numpy.random.rand() < 0.4)
            self.assertTrue (not qnet.is_sampled_by_task (obs))


    def test_parameters (self):
        theta = self.qnet3.parameters
#        expected = [0.0, 1.0] * 7 # for zero-prob
        expected = [ 1.0 ] * 7
        self.assertEquals (expected, theta)
        
        theta2 = [ 0.25, 0.5, 0.5, 0.5, 1.2, 1.2, 0.05 ]
        self.qnet3.parameters = theta2
        self.assertEquals (theta2, self.qnet3.parameters)
        
        
    def test_sem (self):
        sampling.set_seed (345)
        
        nouter = 25
        nt = 50
        pct = 0.8
        net = self.qnet3
        param_orig = [1.0, 0.2, 0.2, 0.2, 0.05, 0.05, 0.3]
 
        all_mus = []
        for i in xrange(nouter):
#            print numpy.random.get_state()
            net.parameters = param_orig[:]
            sample = net.sample (nt)
            sample.validate()

            observed = sample.subset_by_task (pct)
            print observed

            arrv = net.gibbs_initialize (observed.duplicate())
            arrv.validate()        
#            print arrv
            
            mus,arrv = estimation.sem (net, arrv, burn=25, num_iter=10, return_arrv=False, debug=False)
            avg = map(numpy.mean, zip(*mus))
            all_mus.append (avg)
            
        allavg = map(numpy.mean, zip(*all_mus))
        print zip(allavg, param_orig)
        for smp, true in zip(allavg, param_orig):
            self.assertTrue (abs(smp - true) < 0.1)

    # TODO: Not sure why this fails
    def fixme_test_sem_reporting (self):
        sampling.set_seed (345)

        nt = 50
        pct = 0.5
        net = self.qnet3

        sample = net.sample (nt)
        observed = sample.subset_by_task (pct)
        arrv = net.gibbs_initialize (observed.duplicate())
        arrv.validate()        
            
        burn_mrt = [] 
        mrt = []        
        def gather_mrt (net, arrv, iter):
            mrt.append (qstats.mean_response_time (arrv))        
        def gather_burn_mrt (net, arrv, iter):
            burn_mrt.append (qstats.mean_response_time (arrv))

        mus,arrv = estimation.sem (net, arrv, burn=25, num_iter=10, return_arrv=True, report_fn=gather_mrt, burn_report_fn=gather_burn_mrt)

        self.assertEquals (25, len(burn_mrt))
        self.assertEquals (10, len(mrt))
        
        for i in xrange(10):
            self.assertEquals (mrt[i], qstats.mean_response_time(arrv[i]))
    
    def test_sem_small (self):
        """Don't break if a queue has all unobserved tasks."""
        sampling.set_seed (343)

        nt = 20
        pct = 0.8
        net = self.qnet3

        all_mus = []

        sample = net.sample (nt)
        observed = sample.subset_by_task (pct)
        arrv = net.gibbs_initialize (observed.duplicate())
        arrv.validate()        
        print arrv

        mus,arrv = estimation.sem (net, arrv, burn=25, num_iter=10, return_arrv=False, debug=False)

    # TODO: Compare mean and variance of distribution rather than parameters
    def fixme_test_gamma_sem (self):
        sampling.set_seed (342)

        nouter = 1
        nt = 250
        pct = 0.5
        net = self.gamma3
        param_orig = net.parameters[:]
        print net

        all_mus = []
        for i in xrange(nouter):
#            print numpy.random.get_state()
            net.parameters = param_orig[:]
            sample = net.sample (nt)
            
#            s = []
#            app1 = net.queue_by_name ("APP1")            
#            for evt in sample.events_of_queue (app1):
#                s.append(evt.s)
#            print "MU %s SD %s N %s" % (numpy.mean(s), numpy.std(s), len(s))
            
            observed = sample.subset_by_task (pct)
            arrv = net.gibbs_initialize (observed.duplicate())
            arrv.validate()        
            print arrv

            mus,arrv = estimation.sem (net, arrv, burn=25, num_iter=10, return_arrv=False, debug=False)
            avg = map(numpy.mean, zip(*mus))
            all_mus.append (avg)

        allavg = map(numpy.mean, zip(*all_mus))
        print zip(allavg, param_orig)
        for smp, true in zip(allavg, param_orig):
            self.assertTrue (abs(smp - true) < 0.1)
            
    def test_insert (self):
        sampling.set_seed(765)
        net = self.oneq
        arrv = net.sample(5)
        print arrv
        self.assertEquals (10, len(arrv.all_events()))

        q = net.queue_by_name("THE_Q")

        e = arrivals.Event(90, 90, 1, q, 0.45, 0.0, 2.6, 1)
        arrv.insert (e)

        self.assertEquals (11, len(arrv.all_events()))
        prev = e.previous_by_queue()
        self.assertTrue (prev is not None)
        self.assertEquals (7, prev.eid)
        self.assertTrue (prev.a < e.a)
        self.assertTrue (prev.d < e.d)
        next = e.next_by_queue()
        self.assertTrue (next is not None)
        self.assertEquals (8, next.eid)
        self.assertTrue (e.a < next.a)
        self.assertTrue (e.d < next.d)

        e = arrivals.Event(91, 91, 1, q, 0.75, 0.0, 3.5, 1)
        arrv.insert (e)

        self.assertEquals (12, len(arrv.all_events()))
        prev = e.previous_by_queue()
        self.assertTrue (prev is not None)
        self.assertEquals (9, prev.eid)
        self.assertTrue (prev.a < e.a)
        self.assertTrue (prev.d < e.d)
        next = e.next_by_queue()
        self.assertTrue (next is None)

        e = arrivals.Event(92, 92, 1, q, 0.025, 0.0, 1.5, 1)
        arrv.insert (e)

        self.assertEquals (13, len(arrv.all_events()))
        prev = e.previous_by_queue()
        self.assertTrue (prev is None)
        next = e.next_by_queue()
        self.assertTrue (next is not None)
        self.assertEquals (5, next.eid)
        self.assertTrue (e.a < next.a)
        self.assertTrue (e.d < next.d)

        print arrv

        
    def test_delete (self):
        sampling.set_seed(765)
        net = self.gamma3
        
        arrv = net.sample (100)
        for i in xrange(100):
            if numpy.random.rand() < 0.5:
                arrv.delete_task (i)
                
        arrv.validate()
        
    def test_subset (self):
        sampling.set_seed (443)
        net = self.gamma3
        arrv = net.sample (1000)
        obs = arrv.subset_by_task (0.2)
        
        arrv.validate()
        obs.validate()

        nhidden = obs.num_hidden_events ()
        num_events = obs.num_events ()
        evts_pert = 4
        
        mean = 0.8 * 1000 * evts_pert
        sd = evts_pert * sqrt (1000 * 0.2 * 0.8)
        print "Events hidden = %d (expected %d)" % (nhidden, mean)
        
        self.assertEquals (num_events, arrv.num_events())
        self.assertTrue (abs(nhidden - mean) < 2*sd, "ERROR: %d events hidden, expected 200 pm %s" % (nhidden, 2*sd))

    def test_subset2 (self):
        sampling.set_seed (446)
        net = self.gamma3
        arrv = qnetu.load_arrivals (net, TestQnet.weird_arrv_text)
        obs = arrv.subset_by_task (0.2)

        nhidden = obs.num_hidden_events ()
        num_events = obs.num_events ()

        self.assertEquals (num_events, arrv.num_events())
        self.assertTrue (nhidden > 12, "ERROR: %d events hidden, expected more" % nhidden)
        
    # TODO: Causes bus error
    def fixme_test_near_zero (self):
        sampling.set_seed (789)
        net = self.oneq
        net.parameters = [ 1.0, 1e-10 ]
        
        arrv = net.sample(500)
        obs = arrv.subset_by_task(0.5)
        
        net.gibbs_initialize (obs)
        estimation.sem (net, obs, burn=25, num_iter=25)
        
    def test_mix1_sample (self):
        sampling.set_seed (883)
        
        N = 1000
        arrv = self.mix1.sample (N)
        
        S = 0
        Sv = [0, 0]
        Nv = [0, 0]
        for evt in arrv.events_of_qid (1):
            S += evt.s
            component = evt.variable("COMPONENT_THE_Q")
            Sv[component] += evt.s
            Nv[component] += 1
        S = S / N
        for i in xrange(len(Sv)): Sv[i] = Sv[i] / Nv[i]
        
        print Sv
        
        self.assertTrue (abs(S - 0.830) < 0.01, "Mean mismatch: expected %.4f  actual %.4f" % (0.245, S))
        self.assertTrue (abs(Sv[0] - 0.05) < 0.01, "Mean mismatch: expected %.4f  actual %.4f" % (0.05, Sv[0]))
        self.assertTrue (abs(Sv[1] - 2.0) < 2.0/sqrt(Nv[1]), "Mean mismatch: expected %.4f  actual %.4f tol %.4f" % (2.0, Sv[1], 2*2.0/sqrt(Nv[1])))

    # TODO: Errors
    # Exception exceptions.AttributeError: "'NoneType' object has no attribute 'get'" in 'queues.Queue.service_lpdf' ignored
    def fixme_test_mix1_sem (self):
        nouter = 1
        nt = 1000
        pct = 0.5 #0.75
        net = self.mix1
        sampling.set_seed (342)

        self.do_test_sem (net, nt, pct, nouter)

    # TODO: Exception exceptions.AttributeError: "'NoneType' object has no attribute 'get'" in 'queues.Queue.service_lpdf' ignored
    def fixme_test_lnmix1_sem_small (self):
        sampling.set_seed (342)
        nouter = 1
        nt = 100
        pct = 0.25
        net = self.lnmix1
        self.do_test_sem (net, nt, pct, nouter)

    def fixme_test_lnmix1_sem (self):
        sampling.set_seed (342)
        self.do_test_sem (net=self.lnmix1, nt=1000, pct=0.25, nouter=1)

    def do_test_sem (self, net, nt, pct, nouter):
        param_orig = net.parameters[:]
        print param_orig
        
        
        all_mus = []
        for i in xrange(nouter):
            net.parameters = param_orig[:]
            
            sample = net.sample (nt)

            net.estimate (sample)
            print "GOLD ", net.parameters

            f = open("arrv_true.txt", "w")
            f.write (sample.as_csv())
            f.close()

#            for evt in sample:
#                if evt.queue().name == "THE_Q":
#                    if evt.variable("COMPONENT_THE_Q") == 1:
#                        print evt.dump()
#            print "***********************"

            observed = sample.subset_by_task (pct)
            arrv = net.gibbs_initialize (observed.duplicate())
            arrv.validate()        

            mus,arrv_out = estimation.sem (net, arrv, burn=50, num_iter=50, return_arrv=True, debug=False)
            avg = map(numpy.mean, zip(*mus))
            all_mus.append (avg)
            
            # debug
#            for i in xrange(len(arrv_out)):
#                f = open("arrv%d.txt" % i, "w")
#                f.write (arrv_out[i].as_csv())
#                f.close()

        allavg = map(numpy.mean, zip(*all_mus))
        for tup in zip(allavg, param_orig):
            print "%15.10f%15.10f" % tup
        
        
    def test_mix1_parameters (self):
        net = self.mix1
        param_orig = net.parameters[:]
        net.parameters = param_orig[:]
        param_new = net.parameters[:]
        
        self.assertEquals (len(param_orig), len(param_new))
        for x1,x2 in zip (param_orig, param_new):
            self.assertEquals (x1,x2)

    # TODO: Debug exception
    def fixme_test_mix1_param_stationary (self):
        qnetu.set_seed (3732)
        net = self.mix1
        nt = 100
        nreps = 50
        giter = 0
        pct = 0.5
        self.check_param_stationary (net, nt, nreps, giter, pct)
        print qnet.stats()

    def test_net3_param_stationary (self):
        net = self.qnet3
        nt = 100
        nreps = 10
        giter = 0
        pct = 0.25
        self.check_param_stationary (net, nt, nreps, giter, pct)

    def test_ln_param_stationary (self):
        net = self.ln1
        nt = 500
        nreps = 10
        giter = 10
        pct = 0.9
        self.check_param_stationary (net, nt, nreps, giter, pct)

    def check_param_stationary (self, net, nt, nreps, giter, pct):
        theta0 = net.parameters[:]
        tot = numpy.zeros(len(net.parameters))
        tot_gibbs = numpy.zeros(len(net.parameters))

        for rep in range(nreps):
            net.parameters = theta0[:]

            arrv = net.sample (nt)
            net.estimate (arrv)
            tot += numpy.array(net.parameters)

            obs = arrv.subset_by_task (pct, adapt_fn=copy_evt)
#            initial = net.gibbs_initialize (obs)
            resampled = net.gibbs_resample (obs, burn=giter, num_iter=1, return_arrv=True)
            net.estimate (resampled[-1])

            tot_gibbs += numpy.array (net.parameters)

        tot = tot / nreps
        tot_gibbs = tot_gibbs / nreps

        print
        print "TOT_MC  ",
        for x in tot: print " %.10f" % x,
        print

        print "TOT_GIB ", 
        for x in tot_gibbs: print " %.10f" % x,
        print

    def test_ln1_sem (self):
        qnetu.set_seed (342)

        nouter = 1
        nt = 500
        pct = 0.5 #0.75
        net = self.ln1
        param_orig = net.parameters[:]
        print param_orig
        
        all_mus = []
        for i in xrange(nouter):
            net.parameters = param_orig[:]
            
            sample = net.sample (nt)

            net.estimate (sample)
            print "GOLD ", net.parameters

            f = open("arrv_true.txt", "w")
            f.write (sample.as_csv())
            f.close()

            observed = sample.subset_by_task (pct)
            arrv = net.gibbs_initialize (observed.duplicate())
#            observed = sample.subset_by_task (pct, adapt_fn=copy_evt)
#            arrv =observed.duplicate()
            arrv.validate()        

            mus,arrv_out = estimation.sem (net, arrv, burn=0, num_iter=50, return_arrv=True, debug=False)

            # debug
#            for i in range(len(arrv_out)):
#                f = open("arrv%d.txt" % i, "w")
#                f.write (arrv_out[i].as_csv())
#                f.close()

            avg = map(numpy.mean, zip(*mus))
            all_mus.append (avg)

        allavg = map(numpy.mean, zip(*all_mus))
        print zip(allavg, param_orig)

    def test_gibbs_kernel (self):
      sampling.set_seed (23412832)

      net = self.ini_only
      nt = 5
      N = 1000

      arrv_star = net.sample (nt)
      for evt in arrv_star: 
         if evt.next_by_queue():
            evt.obs_a = evt.obs_d = 0
      print arrv_star

      # OK, so only the last task is observed.  True p(v) is gamma shape=nt

      allLp = []
      arrv0 = arrv_star.duplicate()
      for i in range(N):
         arrv1 = net.slice_resample (arrv0.duplicate(), 1)[-1]
#         print "======"
#         print arrv0
#         print arrv1 #GGG
         kp = net.gibbs_kernel (arrv1, arrv_star)
         print "KP %.5f" % kp
         allLp.append (kp)
         arrv0 = arrv1

      print allLp
      gamma = distributions.Gamma (nt, net.parameters[0])
      ievts = arrv_star.initial_events()
      expected = gamma.lpdf (ievts[-1].d)

      actual = net.log_prob (arrv_star) - (pwfun.logsumexp0 (numpy.array (allLp)) - numpy.log(N))

      self.assertTrue (abs(expected - actual) < 0.01, "Mismatch: log prob %.10f  R-B %.10f" % (expected, actual))

    
    def test_big_initialize (self):
        net = self.qnet3
        all_nt = [ 100, 1000, 5000, 10000, 100000 ]

        for nt in all_nt:
            all = net.sample (nt)
            sub = all.subset_by_task (0.1)

            tmr = mytime.timeit()
            arrv = net.gibbs_initialize (sub)
            tmr.tick ("Gibbs initialize @ %d" % nt)
            
            arrv.validate()
        
        
    qnet3_text = """
    states:
      - name: INITIAL
        queues: [ INITIAL ]
        successors: [ TIER1 ]
        initial: TRUE 
      - name: TIER1
        queues: [ WEB1, WEB2, WEB3 ]
        successors: [ TIER2 ]
      - name: TIER2
        queues: [ APP1, APP2 ]
        successors: [ TIER3 ]
      - name: TIER3
        queues: [ DB ]
    queues:
      - { name: INITIAL, service: M }
      - { name: WEB1, service: M }
      - { name: WEB2, service: M }
      - { name: WEB3, service: M }
      - { name: APP1, service: M }
      - { name: APP2, service: M }
      - { name: DB, service: M }
    """

    gamma3_text = """
    states:
      - name: INITIAL
        queues: [ INITIAL ]
        successors: [ TIER1 ]
        initial: TRUE 
      - name: TIER1
        queues: [ WEB1, WEB2, WEB3 ]
        successors: [ TIER2 ]
      - name: TIER2
        queues: [ APP1, APP2 ]
        successors: [ TIER3 ]
      - name: TIER3
        queues: [ DB ]
    queues:
      - { name: INITIAL, service: [M, 0.04] }
      - { name: WEB1, service: [G, 2, 0.025 ] }
      - { name: WEB2, service: [G, 2, 0.025 ] }
      - { name: WEB3, service: [G, 2, 0.025 ] }
      - { name: APP1, service: [G, 5, 0.001 ] }
      - { name: APP2, service: [G, 5, 0.001 ] }
      - { name: DB, service: [G, 2, 0.01] }
    """
    
    ln1_text = """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: TIER1
  queues: [ TIER1_0, TIER1_1 ]
queues:
- { name: INITIAL, service: [G, 20, 0.1] }
- { name: TIER1_0, service: [LN, 0.75, 0.5 ] }
- { name: TIER1_1, service: [LN, 0.75, 0.5 ] }
    """
#- { name: TIER1_1, service: [LN, 0.75, 0.5 ] }

    ln2_text = """
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
- { name: INITIAL, service: [G, 2, 0.01] }
- { name: TIER1_0, service: [LN, -3.75, 0.5 ] }
- { name: TIER1_1, service: [LN, -3.75, 0.5 ] }
- { name: TIER2_0, service: [LN, -3.75, 0.5 ] }
- { name: TIER2_1, service: [LN, -3.75, 0.5 ] }
- { name: TIER3_0, service: [LN, -3.75, 0.5 ] }
- { name: TIER3_1, service: [LN, -3.75, 0.5 ] }
    """
    
    oneq_text = """
    states:
      - name: INITIAL
        queues: [ INITIAL ]
        successors: [ STATE1 ]
        initial: TRUE 
      - name: STATE1
        queues: [ THE_Q ]
    queues:
      - { name: INITIAL, service: [M, 0.1 ] }
      - { name: THE_Q, service: M, type="GGk" }
    """

    ini_only_text = """
    states:
      - name: INITIAL
        queues: [ INITIAL ]
    queues:
      - { name: INITIAL, service: [M, 1.0 ] }
    """
    
    
    mix1_text = """
    states:
      - name: INITIAL
        queues: [ INITIAL ]
        successors: [ STATE1 ]
        initial: TRUE 
      - name: STATE1
        queues: [ THE_Q ]
    queues:
      - { name: INITIAL, service: [M, 2.0 ] }
      - { name: THE_Q, service: [MIX, 0.6, [M, 0.05], 0.4, [M, 2.0] ] }
    """
    
    lnmix1_text = """
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
  - { name: INITIAL, service: [G, 2, 10.0 ] }
  - { name: TIER1_0, service: [MIX, 0.5, [LN, -1, 1 ], 0.333, [LN, 2.5, 3], 0.167, [ LN, 5, 2 ] ] }
  - { name: TIER1_1, service: [MIX, 0.5, [LN, -1, 1 ], 0.333, [LN, 2.5, 3], 0.167, [ LN, 5, 2 ] ] }
  - { name: TIER2_0, service: [MIX, 0.5, [LN, -1, 1 ], 0.333, [LN, 2.5, 3], 0.167, [ LN, 5, 2 ] ] }
  - { name: TIER2_1, service: [MIX, 0.5, [LN, -1, 1 ], 0.333, [LN, 2.5, 3], 0.167, [ LN, 5, 2 ] ] }
  - { name: TIER3_0, service: [MIX, 0.5, [LN, -1, 1 ], 0.333, [LN, 2.5, 3], 0.167, [ LN, 5, 2 ] ] }
  - { name: TIER3_1, service: [MIX, 0.5, [LN, -1, 1 ], 0.333, [LN, 2.5, 3], 0.167, [ LN, 5, 2 ] ] }


"""
    poisson_text = """
    states:
      - name: INITIAL
        queues: [ INITIAL ]
        successors: [ STATE1 ]
        initial: TRUE 
      - name: STATE1
        queues: [ DUMMY_Q ]
    queues:
      - { name: INITIAL, service: [M, 0.5] }
      - { name: DUMMY_Q, service: [M, 0.00001 ] }
    """
    
    m2_text = """
    states:
        - name: INITIAL
          queues: [ INTIAL ]
          successors: [ STATE1 ]
          initial: TURE
        - name: STATE1
          queues: [ Q1 ]
          successors: [ STATE2 ]
        - name: STATE2
          queues: [ Q2 ]
    queues:
       - { name: INITIAL, service: [D, 10.0] }
       - { name: Q1, service: [M, 0.5  ] }
       - { name: Q2, service: [M, 0.25 ] }
    """
    
    
    m2b_text = """
    states:
        - name: INITIAL
          queues: [ INTIAL ]
          successors: [ STATE1 ]
          initial: TRUE
        - name: STATE1
          queues: [ Q1 ]
          successors: [ STATE2 ]
        - name: STATE2
          queues: [ Q2 ]
    queues:
       - { name: INITIAL, service: [D, 10.0] }
       - { name: Q1, service: [M, 0.5  ] }
       - { name: Q2, service: [M, 0.5 ] }
    """
    
    rails3_text = """
    states:
     - name: INITIAL
       queues: [INITIAL]
       successors: [NET]
       initial: TRUE
     - name: NET
       queues: [NET]
       successors: [RAILS]
     - name: RAILS
       queues: [r22_4100,r22_4101,r22_4102,r22_4103,r22_4104,r22_4105,r22_4106,r22_4107,r22_4108,r22_4109]
       successors: [DB]
     - name: DB
       queues: [DB]
    queues:
     - { name: INITIAL, service: M }
     - { name: NET, service: M }
     - { name: r22_4100, service: M }
     - { name: r22_4101, service: M }
     - { name: r22_4102, service: M }
     - { name: r22_4103, service: M }
     - { name: r22_4104, service: M }
     - { name: r22_4105, service: M }
     - { name: r22_4106, service: M }
     - { name: r22_4107, service: M }
     - { name: r22_4108, service: M }
     - { name: r22_4109, service: M }
     - { name: DB, service: M }
    """
    
    
    mm3_text = """
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
      - { name: WEB1, processors: 3, service: [M, 30.0] }
      - { name: APP1, processors: 4, service: [M, 30.0] }
    """

    rails3_arrv_text = """!Arrivals
    events:
    - !Event {arrival: 0.0, departure: 440.0, queue: INITIAL, state: 0, task: 16 }
    - !Event {arrival: 0.0, departure: 612.0, queue: INITIAL, state: 0, task: 15 }
    - !Event {arrival: 0.0, departure: 813.0, queue: INITIAL, state: 0, task: 7 }
    - !Event {arrival: 0.0, departure: 980.0, queue: INITIAL, state: 0, task: 13 }
    - !Event {arrival: 0.0, departure: 1040.0, queue: INITIAL, state: 0, task: 8 }
    - !Event {arrival: 0.0, departure: 1108.0, queue: INITIAL, state: 0, task: 4 }
    - !Event {arrival: 1108.0, departure: 1108.0, queue: NET, state: 1, task: 4}
    - !Event {arrival: 1108.0, departure: 1108.23395, queue: r22_4107, state: 1, task: 4}
    - !Event {arrival: 1108.23395, departure: 1108.3109999999999, queue: DB, state: 2, task: 4}
    - !Event {arrival: 1108.3109999999999, departure: 1108.3130000000001, queue: NET, state: 3,
      task: 4}
    - !Event {arrival: 813.0, departure: 813.0, queue: NET, state: 1, task: 7}
    - !Event {arrival: 813.0, departure: 813.01647000000003, queue: r22_4104, state: 1, task: 7}
    - !Event {arrival: 813.01647000000003, departure: 813.02099999999996, queue: DB, state: 2,
      task: 7}
    - !Event {arrival: 813.02099999999996, departure: 813.02200000000005, queue: NET, state: 3,
      task: 7} 
    - !Event {arrival: 1040.0, departure: 1040.0, queue: NET, state: 1, task: 8}
    - !Event {arrival: 1040.0, departure: 1040.0087599999999, queue: r22_4104, state: 1, task: 8}
    - !Event {arrival: 1040.0087599999999, departure: 1040.009, queue: DB, state: 2, task: 8}
    - !Event {arrival: 1040.009, departure: 1040.009, queue: NET, state: 3, task: 8}
    - !Event {arrival: 980.0, departure: 980.0, queue: NET, state: 1, task: 13}
    - !Event {arrival: 980.0, departure: 981.58758, queue: r22_4105, state: 1, task: 13}
    - !Event {arrival: 981.58758, departure: 981.63300000000004, queue: DB, state: 2, task: 13}
    - !Event {arrival: 981.63300000000004, departure: 981.75, queue: NET, state: 3, task: 13}
    - !Event {arrival: 612.0, departure: 612.0, queue: NET, state: 1, task: 15}
    - !Event {arrival: 612.0, departure: 612.10739999999998, queue: r22_4104, state: 1, task: 15}
    - !Event {arrival: 612.10739999999998, departure: 612.14700000000005, queue: DB, state: 2,
      task: 15}
    - !Event {arrival: 612.14700000000005, departure: 612.16700000000003, queue: NET, state: 3,
      task: 15}
    - !Event {arrival: 440.0, departure: 440.0, queue: NET, state: 1, task: 16}
    - !Event {arrival: 440.0, departure: 440.01567, queue: r22_4103, state: 1, task: 16}
    - !Event {arrival: 440.01567, departure: 440.02199999999999, queue: DB, state: 2, task: 16}
    """
    
    weird_arrv_text = """!Arrivals
    events:
    - !Event {arrival: 0.0, departure: 0.10452631782595653, obs: 1, queue: INITIAL, state: 0,
      task: 212}
    - !Event {arrival: 0.0, departure: 0.11646732583049313, obs: 1, queue: INITIAL, state: 0,
      task: 673}
    - !Event {arrival: 0.0, departure: 0.1384136306061286, obs: 1, queue: INITIAL, state: 0,
      task: 1729}
    - !Event {arrival: 0.0, departure: 0.16351009761370786, obs: 1, queue: INITIAL, state: 0,
      task: 1492}
    - !Event {arrival: 0.0, departure: 0.21410547922073025, obs: 1, queue: INITIAL, state: 0,
      task: 4}
    - !Event {arrival: 0.10452631782595653, departure: 0.17868226625405026, obs: 1, queue: WEB2,
      state: 1, task: 212}
    - !Event {arrival: 0.11646732583049313, departure: 0.21074134333605588, obs: 1, queue: WEB2,
      state: 1, task: 673}
    - !Event {arrival: 0.1384136306061286, departure: 0.16055713978366762, obs: 1, queue: WEB3,
      state: 1, task: 1729}
    - !Event {arrival: 0.16055713978366762, departure: 0.16302537577738552, obs: 1, queue: APP1,
      state: 2, task: 1729}
    - !Event {arrival: 0.16302537577738552, departure: 0.16397850142187836, obs: 1, queue: DB,
      state: 3, task: 1729}
    - !Event {arrival: 0.16351009761370786, departure: 0.19647926638923524, obs: 1, queue: WEB3,
      state: 1, task: 1492}
    - !Event {arrival: 0.17868226625405026, departure: 0.18375190168790723, obs: 1, queue: APP1,
      state: 2, task: 212}
    - !Event {arrival: 0.18375190168790723, departure: 0.21120474695589864, obs: 1, queue: DB,
      state: 3, task: 212}
    - !Event {arrival: 0.19647926638923524, departure: 0.20513038899218811, obs: 1, queue: APP2,
      state: 2, task: 1492}
    - !Event {arrival: 0.20513038899218811, departure: 0.22367865683238419, obs: 1, queue: DB,
      state: 3, task: 1492}
    - !Event {arrival: 0.21074134333605588, departure: 0.21558881350109046, obs: 1, queue: APP2,
      state: 2, task: 673}
    - !Event {arrival: 0.21410547922073025, departure: 0.23630202405888817, obs: 1, queue: WEB3,
      state: 1, task: 4}
    - !Event {arrival: 0.21558881350109046, departure: 0.25280315977450724, obs: 1, queue: DB,
      state: 3, task: 673}
    - !Event {arrival: 0.23630202405888817, departure: 0.23804664749304735, obs: 1, queue: APP2,
      state: 2, task: 4}
    - !Event {arrival: 0.23804664749304735, departure: 0.26553594352433596, obs: 1, queue: DB,
      state: 3, task: 4}

    """
    
    def setUp (self):
        self.qnet3 = qnetu.qnet_from_text (TestQnet.qnet3_text)
        self.gamma3 = qnetu.qnet_from_text (TestQnet.gamma3_text)
        self.oneq = qnetu.qnet_from_text (TestQnet.oneq_text)
        self.ini_only = qnetu.qnet_from_text (TestQnet.ini_only_text)
        self.mix1 = qnetu.qnet_from_text (TestQnet.mix1_text)
        self.lnmix1 = qnetu.qnet_from_text (TestQnet.lnmix1_text)
        self.m2 = qnetu.qnet_from_text (TestQnet.m2_text)
        self.m2b = qnetu.qnet_from_text (TestQnet.m2b_text)
        self.poisson = qnetu.qnet_from_text (TestQnet.poisson_text)
        self.mm3 = qnetu.qnet_from_text (TestQnet.mm3_text)
        self.rails3 = qnetu.qnet_from_text (TestQnet.rails3_text)
        self.ln1 = qnetu.qnet_from_text (TestQnet.ln1_text)
        self.rails3_arrv = qnetu.load_arrivals (self.rails3, TestQnet.rails3_arrv_text)
        
def zero_departure (e_old, e):
    e.a = e_old.a
    w = e_old.wait()
    e.s = 0
    e.d = e.a + w

def copy_evt (e_old, e):
    e.a = e_old.a
    e.d = e_old.d
    e.s = e_old.s
    e.proc = e_old.proc
    e.set_all_vars (e_old)
    e.copy_caches (e_old)

# assumes v sorted
def summarize (v):
    l = len(v)
    return "   Mean    Min     1Q    Med     3Q    Max\n%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f" % (numpy.mean(v), v[0], v[l/4], v[l/2], v[3*l/4], v[-1])
    
def main():
    if len(sys.argv) > 1:
        test_name = sys.argv[1]
        suite = unittest.TestLoader().loadTestsFromName("test_qnet.TestQnet.%s" % (test_name,))
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()

