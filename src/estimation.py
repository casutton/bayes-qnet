# Python file fo estimation procedures.  These tend to spend most of their time in functions
#  from qnet, so it's not as important that they be fast.

import qnet
import mytime
import numpy
import sampling
import qstats
import pwfun
from numpy import array
from numpy import random
import sys

# customization

MH = 0
SLICE = 1
_sampling_type = SLICE

def set_sampler (str):
    global _sampling_type
    if str == "SLICE": _sampling_type = SLICE
    elif str == "MH": _sampling_type = MH
    else: raise Exception ("Could not understand estimator %s" % str)
    print "StEM sampling type = ", str


# Stochastic EM
def sem (net, arrv, burn, num_iter, gibbs_iter=1, momentum=0, return_arrv=False, report_fn=None, burn_report_fn=None, debug=False):
   
    print "SEM: Starting..."

    tmr = mytime.timeit()
    
    total_iter = burn + num_iter
    mus = []
    all_sampled = []
    mu_old = None

    qnet.reset_stats()
    
    # initialize from the Gibbs initialization
    net.estimate (arrv)

    print "MEAN SERVICE (ALL) ", qstats.mean_service (arrv)
    print "MEAN SERVICE (OBS) ", qstats.mean_obs_service (arrv)

    print "MU 0",
    for x in net.parameters: print " %s" % x,
    print
        
    for i in xrange(total_iter):        

        if _sampling_type == MH:
            arrvl = net.gibbs_resample (arrv, burn=0, num_iter=gibbs_iter, return_arrv=False)
        elif _sampling_type == SLICE:
            arrvl,lp = net.slice_resample (arrv, gibbs_iter, return_arrv=False, return_lp=True)
        else: raise Exception ("huh?")

        assert 1 == len(arrvl)
        net.estimate (arrvl)
        arrv = arrvl[-1]

        if momentum > 0 and mu_old:
            params = momentum*array(mu_old) + (1-momentum)*array(net.parameters)
            net.parameters = list(params)

        if debug: print arrvl[-1]
        print "MU %d" % (i+1),
        for theta in net.parameters: print " %s" % theta,
        print
        sys.stdout.flush()
        
        if burn <= i:
            mus.append (net.parameters[:])
            if report_fn: report_fn (net, arrvl[-1], i, lp)
            if return_arrv: all_sampled.append (arrvl[-1])
        else:
            if burn_report_fn: burn_report_fn (net, arrvl[-1], i)
    
        mu_old = net.parameters[:]

    # always return the last arrivals, even in return_arrv false
    if not return_arrv:
        all_sampled.append (arrvl[-1])
        
    # diagnostics
    nobs = arrv.num_hidden_events ()
    nvars = nobs * total_iter
    
    time = tmr.total ("SEM time")
    stats = qnet.stats()
    if nvars > 0:
        print "Num events sampled = ", stats['nevents']
        print "Num M/H rejections = ", stats['nreject']
        print "Events per sec = ", (stats['nevents'] / time)
        print "Num diff lists = ", stats['num_diff']
        print "Avg diff list length = ", (stats['diff_len'] / float(stats['num_diff']))

    return (mus, all_sampled)
 


def bayes (net, arrv, num_iter, report_fn=None, return_arrv=False, reverse=False, sweeten=0, tag="MU"):
    
    tmr = mytime.timeit()
    
    mus = []
    all_sampled = []

    qnet.reset_stats()
    
    for i in range(sweeten):
        if _sampling_type == MH:
            arrvl = net.gibbs_resample (arrv, burn=0, num_iter=1, return_arrv=False)
        elif _sampling_type == SLICE:
            arrvl = net.slice_resample (arrv, 1, return_arrv=False, reverse=reverse)
        else: raise Exception ("Huh?")
        arrv = arrvl[-1]
        print "SWEETEN ", qstats.mean_service(arrv)
    
    # initialize from the Gibbs initialization
    net.estimate (arrv)

    print tag, 0,
    for x in net.parameters: print x,
    print
    sys.stdout.flush()

    for i in xrange(num_iter):        
        if _sampling_type == MH:
            arrvl = net.gibbs_resample (arrv, burn=0, num_iter=1, return_arrv=False)
        elif _sampling_type == SLICE:
            arrvl,lp = net.slice_resample (arrv, 1, return_arrv=False, reverse=reverse, return_lp = True)
        else: raise Exception ("Huh?")

        assert 1 == len(arrvl)
        arrv = arrvl[-1]
        arrv.validate()

        if report_fn: report_fn (net, arrvl[-1], i, lp)
        if return_arrv: all_sampled.append (arrvl[-1])

        net.sample_parameters (arrv)

        print tag, i+1,
        for theta in net.parameters: print theta,
        print
        sys.stdout.flush()
        
        mus.append (net.parameters[:])

    if not return_arrv:
        all_sampled = [arrv]

    # diagnostics
    nobs = arrv.num_hidden_events ()
    nvars = nobs * num_iter
    
    time = tmr.total ("Bayes time")
    stats = qnet.stats()
    if nvars > 0:
        print "Num events sampled = ", stats['nevents']
        print "Num M/H rejections = ", stats['nreject']
        print "Events per sec = ", (stats['nevents'] / time)
        print "Num diff lists = ", stats['num_diff']
        print "Avg diff list length = ", (stats['diff_len'] / float(stats['num_diff']))

    return mus,all_sampled

def sample_departures (net, arrv, num_iter, report_fn=None, debug=False):
   
    print "SEM: Starting..."

    tmr = mytime.timeit()
    
    all_sampled = []

    qnet.reset_stats()

    print "MEAN SERVICE (ALL) ", qstats.mean_service (arrv)
    print "MEAN SERVICE (OBS) ", qstats.mean_obs_service (arrv)

    print "MU 0",
    for x in net.parameters: print " %s" % x,        
    mus = []
    
    for i in xrange(num_iter):
        arrv = arrv.duplicate()
        
        if _sampling_type == MH:
            arrvl = net.gibbs_resample (arrv, burn=0, num_iter=1, return_arrv=False)
        elif _sampling_type == SLICE:
            arrvl,lp = net.slice_resample (arrv, 1, return_arrv=False, return_lp=True)
        else: raise Exception ("huh?")

        assert 1 == len(arrvl)
        arrv = arrvl[-1]

        if debug: print arrvl[-1]
        print "THETA %d" % (i+1),
        print " ".join (map (str, net.parameters))
        print "LP %d %.5f" % (i+1, net.log_prob (arrv))
        sys.stdout.flush()

        mus.append (net.parameters[:])
        if report_fn: report_fn (net, arrvl[-1], i, lp)
        all_sampled.append (arrv)

    # diagnostics
    nobs = arrv.num_hidden_events ()
    nvars = nobs * num_iter
    
    time = tmr.total ("Sampling time")
    stats = qnet.stats()
    if nvars > 0:
        print "Num events sampled = ", stats['nevents']
        print "Num M/H rejections = ", stats['nreject']
        print "Events per sec = ", (stats['nevents'] / time)
        print "Num diff lists = ", stats['num_diff']
        print "Avg diff list length = ", (stats['diff_len'] / float(stats['num_diff']))

    return (mus, all_sampled)
 

# Use report fns

def mean_service_fn (net, arrv, gi):
    mus = qstats.mean_service(arrv)
    print "MEAN_SERVICE     %d %s" % (gi, " ".join (map(str, mus)))
    mobs = qstats.mean_obs_service(arrv)
    print "MEAN_OBS_SERVICE %d %s" % (gi, " ".join (map(str, mobs)))
#    f = open ("arrv%d.txt" % gi, "w")
#    f.write (str(arrv))
#    f.close()

i_ggrf = 0
def gen_graphing_report_fn (pct, prefix="gibbs-dist"):
    def report_fn (net, arrv, gi):
        global i_ggrf
        for evt in arrv:
            if random.rand() < pct:
                # graph it
                f = open ("%s%d.txt" % (prefix, i_ggrf), "w")
                dfn = generate_true_dfn (net, arrv, evt)
                L,U = dfn_range (net, arrv, evt)
                eps = (U-L)/100
                x = L
                for i in range(100):
                    f.write ("%.5f %.5f\n" % (x, dfn(x)))
                    x += eps
                f.close()
                i_ggrf += 1
    return report_fn



# For model selection

def evidence_logProb_IS (net, arrv, niter=100, report_fn=None):
    allLp = []
    def mini_evidence_prob (net, arrv, iter, lp):
        arrv1 = net.slice_resample (arrv.duplicate(), 1)[-1]
        jp = net.log_prob (arrv1)
        kp = net.gibbs_kernel (arrv, arrv1)
        print "JP %.5f  KP %.5f" % (jp, kp) #GGG
        allLp.append (jp - kp)
    bayes (net, arrv, niter, report_fn=mini_evidence_prob)
    print allLp #GGG
    return pwfun.logsumexp0 (numpy.array (allLp)) - numpy.log(len(allLp))

# computes chib approximation to the evidence probability,
#  with the extension of Murray and Salakhutdinov (NIPS 2008)
def chib_evidence (net, obs, M_special=100, M=100, debug=False):

    tmr = mytime.timeit()

    # 1. find the special hidden state
    reporter = Reporter(M_special / 2)
    bayes (net, obs, M_special, report_fn=reporter.outputter, tag="CHIB1")
    
    arrv_star = reporter.best_arrv
    theta_star = reporter.best_theta
    
    tmr.tick ("CHIB step 1 :: High prob state")

    f = open ("chib_star.txt", "w")
    arrv_star.write_csv (f)
    f.close()

    # 2. forward sampling
    u = random.randint(M)
    cndProb = [None] * M
    def inc_prob (net, arrv, iter, lp):
        cndProb[u+iter] = net.gibbs_kernel (arrv, arrv_star, lp)
        cndProb[u+iter] += net.parameter_kernel (net.parameters, theta_star, arrv)
        print "TOT_S", qstats.total_service (arrv) # useful for exponential
        if debug:
            f = open ("chib2_dbg_%d.txt" % iter, "w")
            arrv.write_csv (f)
            f.close()
    bayes (net, arrv_star.duplicate(), (M-u), report_fn=inc_prob, tag="CHIB2")

    tmr.tick ("CHIB step 2 :: Forward state")

    # 3. backward sampling
    def inc_prob2 (net, arrv, iter, lp):
        cndProb[u-iter-1] = net.gibbs_kernel (arrv, arrv_star, lp)
        cndProb[u-iter-1] += net.parameter_kernel (net.parameters, theta_star, arrv)
        print "TOT_S", qstats.total_service (arrv) # useful for exponential
        if debug:
            f = open ("chib3_dbg_%d.txt" % iter, "w")
            arrv.write_csv (f)
            f.close()
    bayes (net, arrv_star.duplicate(), u, report_fn=inc_prob2, reverse=True, tag="CHIB3")

    tmr.tick ("CHIB step 3 :: Forward state")

    # 4. compute results
    #  trick : p(v) = p(h*, v) / p(h* | v)
    p_v_given_hstar = pwfun.logsumexp0 (numpy.array (cndProb)) - numpy.log (M)
    p_v_hstar = net.log_prob (arrv_star)
    # FIXME: assumes improper prior here

    tmr.tick ("CHIB total")
    
    print "CNDPROB", "\n".join (map (str, cndProb)), "/CNDPROB"
    print "CHIB log(p (v, h*)) = %.5f\nCHIB log(p (v|h*)) = %.5f" % (p_v_hstar, p_v_given_hstar)

    return p_v_hstar - p_v_given_hstar  


class Reporter:

   def __init__ (self, burn=0, allStats=[], train_stats_f=None, outputArrv=False, arrvIter=1):
       self.allStats = allStats
       self.train_stats_f = train_stats_f
       self.outputArrv = outputArrv
       self.arrvIter = arrvIter
       self.best_iter = -1
       self.max_lik = -numpy.inf
       self.burn=burn

   def outputter (self, net, arrv, iter, lp):
        if self.outputArrv:
            if iter % self.arrvIter == 0:
                arrv_out_f = open ("%s%d.txt" % (self.outputArrv, iter), 'w')
                arrv.write_csv (arrv_out_f)
                arrv_out_f.write ("\n")
                arrv_out_f.close ()

        # output statistics for training arrvials
        if self.train_stats_f:
            for i in range(len(self.allStats)):
                f = self.train_stats_f[i]
                fn = self.allStats[i]
                f.write ("%d " % iter)
                f.write (qstats.stringify(fn(arrv)))
                f.write ("\n")

        lik = net.log_prob (arrv)
        print "ITER %d LP %.10f" % (iter, lik)
        if (iter > self.burn) and  (lik > self.max_lik): 
            self.max_lik = lik
            self.best_iter = iter
            self.best_arrv = arrv.duplicate()
            self.best_theta = net.parameters[:]


