#!/usr/bin/python
#  Front end for running the Gibbs sampler


import getopt, sys
import numpy
import numpy.random
import yaml

import expt
import mytime
import estimation, qnetu
import sampling

import qstats

def usage():    
    print """Usage:
     gibbs.py <options>
    Front-end for running the queueing network Gibbs sampler.
     --network (-n)       YAML file describing network (required)
     --arrivals (-a)      YAML file describing arrivals (required)
     --multiarr (-A)      Three-file prefix describing arrivals (alternative to -a)
     --arrvtbl            Arrivals file in table format
     --burn <int>         Number of iterations to burn in (default 100)    
     --bayes <T|F>        If true, use Bayes instead of StEM
     --departures-only    If supplied, hold params constant; only sample departure times
     --gibbs-iter <int>   Number of iterations to sample (default 100)
     --output-mu <file>   If supplied, write series of average response times to file
     --output-arrv <prefix> If supplied, write all arrivals in YAML to files with given prefix 
     --arrv-iter <int>    Number of iterations at which to output arrivals
     --thin <int>         Number of iterations to thin in stochastic EM
     --statistics <list>  List of statistics to report after every iteration.
                          Example: --stat "mean_wait, mean_qsize_arriving"
     --do-init <int>      Whether to run initialization (True by default)
                          Use 0 if you're restarting from a previous sample
     --params FILE        Read in parameters from file
         Available statistics: %s
    """ % qstats.all_stats_string()
    
def main():
    
    nburn = 100
    niter = 100
    netf = None
    arrf = None
    multif = None
    arrvtbl = None
    samplePct = None
    outputMu = outputArrv = None
    arrvIter = 100
    configFile = None
    allStats = []
    statsIter = 0
    thin = 1
    doBayes=0
    doInit=1
    params=None
    departuresOnly=False
    sweeten=0

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hn:a:A:", ["help", "seed=", "network=", "arrivals=", "gibbs-iter=", "burn=", "output-mu=", "output-arrv=", "config=", "statistics=", "stats=", "thin=", "stats-iter=", "multiarr=", "arrvtbl=", "bayes=", "departures-only", "arrv-iter=", "do-init=", "sweeten=", "params=" ])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    output = None
    verbose = False
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-n", "--network"):
            netf = a
        elif o in ("-a", "--arrivals"):
            arrf = a
        elif o in ("-A", "--multiarr"):
            multif = a
        elif o == "--arrvtbl":
            arrvtbl = a
        elif o == "--do-init":
            doInit = int(a)
        elif o == "--burn":
            nburn = int(a)
        elif o == "--sweeten":
            sweeten = int(a)
        elif o == "--gibbs-iter":
            niter = int(a)
        elif o == "--output-mu":
            outputMu = a
        elif o == "--output-arrv":
            outputArrv = a
        elif o == "--arrv-iter":
            arrvIter = int(a)
        elif o == "--thin":
            thin = int(a)
        elif o == "--seed":
            sampling.set_seed (int(a))
        elif o == "--config":
            configFile = a
        elif o == "--departures-only":
            departuresOnly = True
        elif o == "--statistics" or o == "--stats":
            if a == "ALL":
                allStats = qstats.STATS[:]                
            else:
                statl = a.split (", ")
                allStats = [ qstats.__dict__[fn_name] for fn_name in statl ]
        elif o == "--stats-iter":
            statsIter = int(a)
        elif o == "--bayes":
            doBayes = int(a)
        elif o == "--params":
            params = a 
        else:
            assert False, "unhandled option"
    
    tmr = mytime.timeit()
    
    if configFile:
        f = open(configFile)
        qnetu.read_configuration (f)
        f.close()
        
    if netf is None: 
        print "ERROR: Specify --network"
        sys.exit(1)

    f = open(netf)
    net = qnetu.qnet_from_text (f)
    f.close()
    
    if params:
        f = open(params)
        line = f.readlines()[0]
        mu = map(float, line.split())
        net.parameters = mu

    if arrf is not None:
        f = open(arrf)
        arrv = qnetu.load_arrivals (net, f)
        f.close()
    elif multif is not None:
        statef = open ("%sstate.txt" % multif)
        af = open ("%sa.txt" % multif)
        df = open ("%sd.txt" % multif)
        arrv = qnetu.read_multifile_arrv (net, statef, af, df)
        statef.close()
        af.close()
        df.close()
    elif arrvtbl is not None:
        arrv = qnetu.read_from_table (net, arrvtbl)
        # hack
        print "WARNING: Hackily resampling mixture components"
        for evt in arrv:
            evt.queue().resample_auxiliaries(evt)
    else:
        raise Exception ("Must specify either -a or -A")
    
    tmr.tick ("Loading arrivals")
#       print arrv.dump()
    
    print net
    print "Number of events = ", arrv.num_events ()
    print "Number hidden events = ", arrv.num_hidden_events()
    print "Thin = ", thin

    if doInit:
        arrv = net.gibbs_initialize (arrv)
    if outputArrv:
        inif = open ("initial.txt", 'w')
        arrv.write_csv (inif)
        inif.write ("\n")
        inif.close ()
#    print "INITIALIZATION"
#    print arrv
    tmr.tick ("Initialization")

    arrv.validate()
    tmr.tick ("Validation")

    train_stats_f = []
    if outputMu: mu_f = open (outputMu, 'w')
    for fn in allStats: train_stats_f.append (open ("train_%s.txt" % fn.__name__, "w"))
    
    def output_train_stats (net, arrv, iter, lp):
        if outputMu:
            mu_f.write (" ".join (map (str, net.parameters)))
            mu_f.write ("\n")
        
        if outputArrv:
            if iter % arrvIter == 0:
                arrv_out_f = open ("%s%d.txt" % (outputArrv, iter), 'w')
                arrv.write_csv (arrv_out_f)
                arrv_out_f.write ("\n")
                arrv_out_f.close ()
                
        # output statistics for training arrvials
        for i in range(len(allStats)):
            f = train_stats_f[i]
            fn = allStats[i]
            f.write ("%d " % iter)
            f.write (qstats.stringify(fn(arrv)))
            f.write ("\n")
            f.flush ()
    
    if doBayes:
        mul,arrvl = estimation.bayes (net, arrv, niter, report_fn=output_train_stats, sweeten=sweeten)
    elif departuresOnly:
        mul,arrvl = estimation.sample_departures (net, arrv, niter, report_fn=output_train_stats)
    else:
        mul,arrvl = estimation.sem (net, arrv, nburn, niter, gibbs_iter=thin, report_fn=output_train_stats)

    tmr.tick ("Running EM")

    if outputMu: mu_f.close()
    for f in train_stats_f: f.close()
    
    # output statistics from a fresh sample at the ML soln, if desired
    if statsIter > 0:
        test_stats_f = []
        for fn in allStats: test_stats_f.append (open ("%s.txt" % fn.__name__, "w"))

        def output_test_stats (net, arrv, iter):
            for i in range(len(allStats)):
                f = test_stats_f[i]
                fn = allStats[i]
                f.write (qstats.stringify (fn (arrv)))
                f.write ("\n")

        net.parameters = mul[-1]
        arrv = arrvl[-1]
        arrvl2 = net.gibbs_resample(arrv, 0, statsIter, return_arrv= False, report_fn = output_test_stats)
        tmr.tick ("Sampling wait times")

        tmr.tick ("Outputing estimated statistics")

    lp = net.log_prob (arrvl[-1])
    n = numpy.sum (qstats.nobs_by_queue (arrvl[-1]))
    print "LOG_PROB ", lp
    print "AIC ", 2*len(net.parameters) - 2*lp
    print "BIC ", numpy.log(n)*len(net.parameters) -2 * lp

    tmr.total ("Total time")
    
if __name__ == "__main__":
    main()

