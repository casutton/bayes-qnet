#!/usr/bin/python
#  Front end for running the Gibbs sampler


import getopt, sys
import numpy.random

import sampling
import qnetu
import yaml

def usage():
    print """Usage:
     simulate.py <options>
    Front-end for simulating a general queueing network.
    Options:
      --network   YAML file containing network description (required)
      --num-tasks Number of tasks to simulate (default: 10)
      --output-file  File to write YAML arrivals to (default: arrv.txt)
      --subset-pct   If supplied, generate a file containing a sample of this percentage of tasks
      --subset-file  File to write task sumbset to (default: arrv-sampled.txt).  Must supply subset-pct as well.
      --seed         Integer random seed"""
    
    
def main():
    
    subsetPct = None
    outf = "arrv.yml"
    subf = "arrv-sampled.yml"
    numTasks = 10
    useMultif = False
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["help", "network=", "num-tasks=", "output-file=", "subset-pct=", "subset-file=", "seed=", "use-multif=" ])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    output = None
    verbose = False
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("--network"):
            netf = a
        elif o in ("--num-tasks"):
            numTasks = int(a)
        elif o in ("--output-file"):
            outf = a
        elif o in ("--use-multif"):
            useMultif = True
        elif o in ("--subset-file"):
            subf = a
        elif o in ("--subset-pct"):
            subsetPct = float(a)
            if subsetPct <= 0 or subsetPct > 1:
                raise Exception ("Invalid --subset-pct %s Must be in 0..1" % a)
        elif o in ("--seed"):
            sampling.set_seed (int(a))
        else:
            print "Error: Can't understand ", o
            sys.exit (1)
            
    f = open(netf)
    net = qnetu.qnet_from_text (f)
    f.close()

    arrv = net.sample (numTasks)
    
    write_arrivals (arrv, outf, useMultif)

    if subsetPct:
        subset = arrv.subset_by_task (subsetPct, adapt_fn=copy_evt)
        write_arrivals (subset, subf, useMultif)

def write_arrivals (arrv, outf, useMultif):
    if useMultif:
        statef = open ("%sstate.txt" % outf, "w")
        af = open ("%sa.txt" % outf, "w")
        df = open ("%sd.txt" % outf, "w")
        arrv = qnetu.write_multifile_arrv (arrv, statef, af, df, obs_only=False)
        statef.close()
        af.close()
        df.close()
    else:
        arrv = net.sample (numTasks)
        f = open (outf, 'w')
        yaml.dump (arrv, stream=f)
        f.close ()
        
def copy_evt (e_old, e):
    e.a = e_old.a
    e.d = e_old.d
    e.s = e_old.s        
 
if __name__ == "__main__":
    main()

