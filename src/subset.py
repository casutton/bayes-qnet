#!/usr/bin/python
#  Front end for running the Gibbs sampler


import sys
from optparse import OptionParser

import sampling
import qnetu
import yaml

    
def main():
    
    subsetPct = None
    retainAll = True
    
    parser = OptionParser(description="Take a random subset of a given set of arrivals.")
    parser.add_option("--network", dest="netf",
                      help="YAML file containing network description", metavar="FILE")
    parser.add_option("--arrivals", dest="arrvf",
                        help="YAML file describing arrivals", metavar="FILE")
    parser.add_option("--output-file", dest="outf", default="arrivals.yml", 
                         help="File to write YAML arrivals to")
    parser.add_option("--subset-pct", dest="subsetPct", type="float", default=0.5,
                        help="Percentage of tasks in subset.")
    parser.add_option("--max-time", dest="maxTime", default=-1, type="int")
    parser.add_option("--verbose", dest="verbose", action="store_true", default=False,
                        help="If true, ouptut lots of debugging information.")
    parser.add_option("--retain-all", dest="retain_all", action="store_true", default=True,
                        help="If supplied, keep all arrivals in the output file, marking the ones that aren't observed.")
    parser.add_option("--no-retain-all", dest="retain_all", action="store_false",
                        help="If true, keep all arrivals in the output file, marking the ones that aren't observed.")
    parser.add_option("--seed", dest="seed", type="int", default=12241, help="Integer random seed.")
    parser.add_option("--use-multif", dest="isMultif", action="store_true",
                      help="If supplied, input and output are in multif format.")
    parser.add_option("--input-multif", dest="inputMultif",
                      help="Prefix for input multif.")
    parser.add_option("--output-multif", dest="outputMultif",
                      help="Prefix for input multif.")

    (options, args) = parser.parse_args()    
            
    sampling.set_seed (options.seed)

    f = open(options.netf)
    net = qnetu.qnet_from_text (f)
    f.close()

    if options.isMultif:
        statef = open ("%sstate.txt" % options.inputMultif)
        af = open ("%sa.txt" % options.inputMultif)
        df = open ("%sd.txt" % options.inputMultif)
        arrv = qnetu.read_multifile_arrv (net, statef, af, df)
        statef.close()
        af.close()
        df.close()
    else:
        f = open (options.arrvf)
        arrv = qnetu.load_arrivals (net, f)
        f.close ()
    
    if options.verbose:
        print "Original arrivals", arrv
        arrv.validate()

    if options.maxTime > 0:
        if options.verbose: print "Subset up to time ", options.maxTime
        def task_in_time (arrv, tid):
            evtl = arrv.events_of_task (tid)
            task_arrival = min(e.d for e in evtl)
            ret = task_arrival < options.maxTime
            return ret
        subset = arrv.subset_by_task_fn (task_in_time)
    elif options.retain_all:
        subset = arrv.subset_by_task (options.subsetPct, adapt_fn=copy_evt)
    else:
        subset = arrv.subset_by_task (options.subsetPct)
        delete_unobserved_tasks (subset)

    if options.verbose:
        subset.validate ()
        print "Subsetted arrivals", subset
                  
    if options.isMultif:
        statef = open ("%sstate.txt" % options.outputMultif, "w")
        af = open ("%sa.txt" % options.outputMultif, "w")
        df = open ("%sd.txt" % options.outputMultif, "w")
        arrv = qnetu.write_multifile_arrv (subset, statef, af, df, obs_only=False)
        statef.close()
        af.close()
        df.close()
    else:
        f = open (options.outf, 'w')
        yaml.dump (subset, stream=f)
        f.close()

def copy_evt (e_old, e):
    e.a = e_old.a
    e.d = e_old.d
    e.s = e_old.s        
    
def delete_unobserved_tasks (subset):
    tids = dict()
    for e in subset: 
        if not e.obs:
            tids[e.tid] = True
    for tid in tids:
        subset.delete_task (tid)
    
if __name__ == "__main__":
    main()

