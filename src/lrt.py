#!/usr/bin/python
#  Front end for model selection


import getopt, sys
from optparse import OptionParser

import numpy
import numpy.random
from scipy.stats import distributions
import yaml

import expt
import mytime
import estimation, qnetu
import sampling

import qstats
    
def main():
    parser = OptionParser(description="Front-end for model selection using LRT.")
    parser.add_option("--network", "-n", dest="allNetf", action="append", default=[],
                      help="YAML file describing network for each hypothesis (required; use once for each network to compare)\nFirst network supplied is used as the null hypothesis", metavar="FILE")
    parser.add_option("--arrivals", dest="allArrvf", action="append", default=list(),
                        help="Prefix for multif files describing arrivals", metavar="FILE")
    parser.add_option("--config", dest="config", help="YAML file containing configuration")
    parser.add_option("--seed", dest="seed", type="int", default=12241, help="Integer random seed.")
    parser.add_option("--use-multif", dest="isMultif", action="store_true",
                      help="If supplied, input and output are in multif format.")
    parser.add_option ("--gibbs-iter", dest="niter", type="int", default=10,
                       help="Number of Gibbs iteration for each model")
    parser.add_option ("--stats-prefix", dest="statsPrefix",
                       help="Prefix of stats files output.")
    parser.add_option ("--output-arrv", dest="outputArrv", action="store_true",
                       help="If supplied, write arrivals after every interation")

    (options, args) = parser.parse_args()    
            
    sampling.set_seed (options.seed)
    
    tmr = mytime.timeit()
    
    if options.config:
        f = open(options.config)
        qnetu.read_configuration (f)
        f.close()

    if len(options.allNetf) == 0:
        parser.error ("No networks specified (use -n)")
        sys.exit(1)

    if len(options.allNetf) != len(options.allArrvf):
        parser.error ("Need to specify same # of arrivals as networks.\nNets: %s\nArrivals %s" % (options.allNetf, options.allArrvf))

    allNet = []
    for netf in options.allNetf:
        f = open(netf)
        allNet.append (qnetu.qnet_from_text (f))
        f.close()
    net0 = allNet[0]

    likf = []
    for i in range(len(allNet)):
        maxLik = computeLik (i, allNet[i], options.allArrvf[i], options, tmr)
        likf.append (maxLik)
    
    print "Model 0 Lik %.5f" % likf[0]
    for i in range(1, len(likf)):
        print "Model %d Lik %.10f" % (i, likf[i])
    print "MDL LOG.LR DF p"
    for i in range(1, len(likf)):
        llr = likf[i] - likf[0]
        df = len(allNet[i].parameters) - len(net0.parameters)
        if df > 0:
            p = 1 - distributions.chi2.cdf (-2*llr, df)
        else:
            p = -1
        print "%d %.5f %d %.5f" % (i, llr, df, p)

    tmr.total ("Total time")
    

def computeLik (mdlIx, net, arrvf, options, tmr):

    statef = open ("%sstate.txt" % arrvf)
    af = open ("%sa.txt" % arrvf)
    df = open ("%sd.txt" % arrvf)
    arrv = qnetu.read_multifile_arrv (net, statef, af, df)
    statef.close()
    af.close()
    df.close()
    
    print "Number of events = ", arrv.num_events ()
    print "Number hidden events = ", arrv.num_hidden_events()

    tmr.tick ("MDL %d : Loading arrivals" % mdlIx)

    arrv = net.gibbs_initialize (arrv.duplicate ())
    tmr.tick ("MDL %d : Initialization" % mdlIx)
    arrv.validate()
    tmr.tick ("MDL %d : Validation" % mdlIx)

    arrv.display()

    allStats = qstats.STATS[:]
    train_stats_f = [ open("train_%d_%s.txt" % (mdlIx, fn.__name__), "w") for fn in allStats ]
#    reporter = Reporter (allStats, train_stats_f, options)
      
    lp = estimation.chib_evidence (net, arrv, options.niter, options.niter)
    tmr.tick ("Running evidence log prob (Model %d)" % mdlIx)

    for f in train_stats_f: f.close()

    print "MDL %d: Evidence prob %.5f" % (mdlIx, lp)

    return lp
    

if __name__ == "__main__":
    main()

