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
import modelmgmt
import sampling
import qnetu
    
def main():
    parser = OptionParser(description="Front-end for greedy bottleneck finding.")
    parser.add_option("--network", "-n", dest="net0f",
                      help="YAML file describing network describing expected behavior", metavar="FILE")
    parser.add_option("--arrivals", dest="arrv0f",
                        help="Prefix for multif files describing arrivals", metavar="FILE_PREFIX")
    parser.add_option("--queue-text", dest="qtext", help="YAML description of the type of queue you want to try adding to the network.")
    parser.add_option("--state", dest="state", help="Name of the state you want to add bottlenecks to.")
    parser.add_option("--config", dest="config", help="YAML file containing configuration")
    parser.add_option("--seed", dest="seed", type="int", default=12241, help="Integer random seed.")
    parser.add_option ("--gibbs-iter", dest="niter", type="int", default=10,
                       help="Number of Gibbs iteration for each model")
    parser.add_option ("--stats-prefix", dest="statsPrefix",
                       help="Prefix of stats files output.")
    parser.add_option ("--output-arrv", dest="outputArrv", action="store_true",
                       help="If supplied, write arrivals after every interation")
    parser.add_option ("--mdlidx", dest="mdlIdx", default=-1, type="int",
                       help="If supplied, try only alternative model # MDLIDX (0 == base network)")

    (options, args) = parser.parse_args()    
            
    sampling.set_seed (options.seed)
    
    tmr = mytime.timeit()
    
    if options.config:
        f = open(options.config)
        qnetu.read_configuration (f)
        f.close()

    if not options.net0f:
        parser.error ("No networks specified (use -n)")
        sys.exit(1)

    if not options.arrv0f:
        parser.error ("No data specified (use --arrivals)")
        sys.exit(1)

    if not options.state:
        parser.error ("No state specified (use --state)")
        sys.exit(1)

    if not options.qtext:
        parser.error ("No queue text specified (use --queue-text)")
        sys.exit(1)

    f = open(options.net0f)
    net0 = qnetu.qnet_from_text (f)
    f.close()

    arrv0 = qnetu.read_multif_of_prefix (options.arrv0f, net0)
    arrv0.validate()
    tmr.tick ("Reading time")

    print "Evaluating model", options.mdlIdx
    modelmgmt.bottleneck_model_selection (
        net0, options.state, arrv0, options.qtext,
        gibbs_iter=options.niter,
        mdlIx=options.mdlIdx
        )

    tmr.total ("Total time")

if __name__ == "__main__":
    main()

